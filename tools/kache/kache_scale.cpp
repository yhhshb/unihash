// kache-hash side of the thread-scaling comparison (see tools/kache_compare.sh).
//
// Reads a dump directory written by `unihash-scale --dump DIR` (the same work
// items and query workloads unihash-scale measures, with their expected
// answers):
//   manifest.tsv   k, l, kmers, distinct, capacity, reps, threads, items,
//                  pos N HITS, neg N HITS, stream SEQS WINDOWS VALID HITS
//   items.txt      work items, one per line (64 kbp windows of long records
//                  overlapping by k - 1 bases, or whole reads)
//   stream.txt     sequences of the streaming query, one per line
//   pos.u64        positive k-mers (sampled from the input, shuffled)
//   neg.u64        uniform random k-mers
// and measures, as unihash-scale does for unihash:
//   1. single-threaded kache-hash (`Streaming_Kmer_Hash_Table<k, false>`),
//      growing from KACHE_INIT_CAP like the single-threaded comparison; then
//      positive, random and streaming queries on it from each thread count
//      (`find` is const and takes no lock; the table is not modified);
//   2. concurrent kache-hash (`Streaming_Kmer_Hash_Table<k, true>`, one token
//      per thread) built with each thread count, the work items handed out
//      dynamically in batches like unihash-scale. With `capacity 1` in the
//      manifest the table is created for the distinct k-mer count (as
//      unihash-scale --capacity pre-sizes the concurrent unihash table);
//      otherwise it grows from KACHE_INIT_CAP.
// Every answer is checked: table sizes against the distinct k-mer count, hit
// counts against the expected answers. Exit status 1 on any mismatch.
//
// Pre-parsed k-mers are decoded to ASCII outside the timed region (kache-hash
// only takes a k-mer through a `Kmer_Window` built from characters); `init`
// (hash and minimizer of the k-mer) is timed, matching unihash recomputing
// the minimizer of a packed k-mer.

#include "Streaming_Kmer_Hash_Table.hpp"

#include <malloc.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#ifndef KH_K
#define KH_K 31
#endif
#ifndef KH_L
#define KH_L 17
#endif

constexpr uint16_t k = KH_K;
constexpr uint16_t l = KH_L;
using TableST = kache_hash::Streaming_Kmer_Hash_Table<k, false, void, l>;
using TableMT = kache_hash::Streaming_Kmer_Hash_Table<k, true, void, l>;
using Window = kache_hash::Kmer_Window<k, l>;

// Work items handed to a thread at a time (as unihash-scale).
constexpr size_t BATCH = 32;

static double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// Live heap bytes according to glibc (all arenas, plus mmapped chunks).
static size_t heap_bytes() {
#if defined(__GLIBC__) && (__GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 33))
    struct mallinfo2 mi = mallinfo2();
    return mi.uordblks + mi.hblkhd;
#else
    struct mallinfo mi = mallinfo();
    return size_t(unsigned(mi.uordblks)) + size_t(unsigned(mi.hblkhd));
#endif
}

static bool is_acgt(char c) { return c == 'A' || c == 'C' || c == 'G' || c == 'T'; }

// Calls f(pointer, length) for every maximal ACGT run of length >= k of `s`
// (upper-cased beforehand).
template <typename F> static void for_each_segment(const std::string& s, F f) {
    size_t i = 0;
    while (i < s.size()) {
        while (i < s.size() && !is_acgt(s[i])) i++;
        size_t j = i;
        while (j < s.size() && is_acgt(s[j])) j++;
        if (j - i >= k) f(s.data() + i, j - i);
        i = j;
    }
}

static std::vector<std::string> read_lines(const std::string& path) {
    std::ifstream f(path);
    if (!f) { std::cerr << "cannot open " << path << "\n"; std::exit(2); }
    std::vector<std::string> v;
    std::string line;
    while (std::getline(f, line)) {
        for (auto& c : line) c = std::toupper(static_cast<unsigned char>(c));
        v.push_back(std::move(line));
    }
    return v;
}

// k-mers of a u64 file, decoded to k characters each.
static std::vector<char> read_decoded(const std::string& path, size_t& n) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) { std::cerr << "cannot open " << path << "\n"; std::exit(2); }
    n = size_t(f.tellg()) / 8;
    std::vector<uint64_t> xs(n);
    f.seekg(0);
    f.read(reinterpret_cast<char*>(xs.data()), n * 8);
    static const char D[4] = {'A', 'C', 'G', 'T'};
    std::vector<char> out(n * k);
    for (size_t i = 0; i < n; i++)
        for (int j = 0; j < k; j++) out[i * k + j] = D[(xs[i] >> (2 * (k - 1 - j))) & 3];
    return out;
}

// Speed of one thread on a fixed integer loop (the same loop as
// unihash-scale's CPU check; M iterations/s).
static double cpu_speed() {
    const uint64_t n = 4000000;
    double t0 = now();
    volatile uint64_t sink;
    uint64_t x = 0x9e3779b97f4a7c15ull;
    for (uint64_t i = 0; i < n; i++) {
        x ^= x >> 31;
        x = x * 0xbf58476d1ce4e5b9ull + i;
    }
    sink = x;
    (void)sink;
    return n / (now() - t0) / 1e6;
}

// One-thread CPU speed now (best of three), relative to the start of the run,
// taken after 0.3 s of spinning: well below 1 means the machine was throttled
// or busy when the row was measured.
static double cpu_ratio() {
    static double base = [] {
        double t0 = now();
        while (now() - t0 < 0.3) cpu_speed();
        double b = 0;
        for (int i = 0; i < 5; i++) b = std::max(b, cpu_speed());
        return b;
    }();
    double c = 0;
    for (int i = 0; i < 3; i++) c = std::max(c, cpu_speed());
    return c / base;
}

static double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

// Runs f(begin, end) over n parts in `threads` threads (contiguous ranges, as
// unihash-scale); returns (sum of results, seconds).
static std::pair<uint64_t, double> parallel_count_once(size_t n, size_t threads,
                                                       const std::function<uint64_t(size_t, size_t)>& f) {
    std::vector<uint64_t> out(threads, 0);
    std::vector<std::thread> ts;
    size_t per = (n + threads - 1) / threads;
    double t0 = now();
    for (size_t i = 0; i < threads; i++)
        ts.emplace_back([&, i] { out[i] = f(std::min(i * per, n), std::min((i + 1) * per, n)); });
    for (auto& t : ts) t.join();
    double s = now() - t0;
    uint64_t sum = 0;
    for (auto x : out) sum += x;
    return {sum, s};
}

// The same `reps` times: (hits of the first run, median seconds); `same` is
// cleared if the runs give different hits.
static std::pair<uint64_t, double> parallel_count(size_t reps, size_t n, size_t threads,
                                                  const std::function<uint64_t(size_t, size_t)>& f,
                                                  bool& same) {
    std::vector<double> ts;
    uint64_t h0 = 0;
    for (size_t r = 0; r < std::max<size_t>(reps, 1); r++) {
        auto [h, t] = parallel_count_once(n, threads, f);
        if (r == 0) h0 = h;
        same &= h == h0;
        ts.push_back(t);
    }
    return {h0, median(ts)};
}

// Inserts every k-mer of `item`.
template <typename Table, typename Token>
static void insert_item(Table& t, const std::string& item, const Token& token) {
    for_each_segment(item, [&](const char* p, size_t len) {
        Window w;
        w.init(p);
        for (size_t i = 0;; i++) {
            t.insert(w, token);
            if (i + k == len) break;
            w.advance(p[i + k]);
        }
    });
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: kache_scale DUMP_DIR [initial_capacity]\n";
        return 2;
    }
    const std::string dir = argv[1];
    const size_t init_cap = argc > 2 ? std::stoull(argv[2]) : (1u << 20);

    uint64_t kmers = 0, distinct = 0, reps = 1, n_items = 0;
    bool capacity = false;
    std::vector<size_t> threads;
    uint64_t pos_n = 0, pos_hits = 0, neg_n = 0, neg_hits = 0;
    uint64_t seq_n = 0, seq_windows = 0, seq_valid = 0, seq_hits = 0;
    {
        std::ifstream m(dir + "/manifest.tsv");
        if (!m) { std::cerr << "no manifest in " << dir << "\n"; return 2; }
        std::string line;
        while (std::getline(m, line)) {
            std::vector<std::string> f;
            std::stringstream ss(line);
            std::string x;
            while (std::getline(ss, x, '\t')) f.push_back(x);
            const auto u = [&](size_t i) { return uint64_t(std::stoull(f.at(i))); };
            if (f[0] == "k" && u(1) != k) {
                std::cerr << "dump has k = " << f[1] << " but this binary was built for k = " << k << "\n";
                return 2;
            }
            if (f[0] == "l" && u(1) != l) {
                std::cerr << "dump has l = " << f[1] << " but this binary was built for l = " << l << "\n";
                return 2;
            }
            if (f[0] == "kmers") kmers = u(1);
            if (f[0] == "distinct") distinct = u(1);
            if (f[0] == "capacity") capacity = u(1) != 0;
            if (f[0] == "reps") reps = std::max<uint64_t>(1, u(1));
            if (f[0] == "items") n_items = u(1);
            if (f[0] == "pos") pos_n = u(1), pos_hits = u(2);
            if (f[0] == "neg") neg_n = u(1), neg_hits = u(2);
            if (f[0] == "stream") seq_n = u(1), seq_windows = u(2), seq_valid = u(3), seq_hits = u(4);
            if (f[0] == "threads") {
                std::stringstream ts(f[1]);
                std::string t;
                while (std::getline(ts, t, ',')) threads.push_back(std::stoull(t));
            }
        }
    }
    const auto items = read_lines(dir + "/items.txt");
    const auto seqs = read_lines(dir + "/stream.txt");
    size_t np = 0, nn = 0;
    const auto pos = read_decoded(dir + "/pos.u64", np);
    const auto neg = read_decoded(dir + "/neg.u64", nn);
    int bad = 0;
    const auto check = [&](bool ok, const std::string& what) {
        if (!ok) {
            std::cout << "  MISMATCH: " << what << "\n";
            bad++;
        }
    };
    check(items.size() == n_items, "number of work items");
    check(seqs.size() == seq_n, "number of streaming sequences");
    check(np == pos_n && nn == neg_n, "number of query k-mers");
    {
        uint64_t n = 0;
        for (auto& r : items) for_each_segment(r, [&](const char*, size_t len) { n += len - k + 1; });
        check(n == kmers, "k-mers in the work items");
    }
    const unsigned cores = std::max(1u, std::thread::hardware_concurrency());

    cpu_ratio();  // take the CPU-speed reference now
    printf("== kache-hash scaling (k = %u, l = %u), %u cores ==\n", k, l, cores);
    printf("input: %lu items, %lu k-mers, %lu distinct\n", (unsigned long)items.size(), (unsigned long)kmers,
           (unsigned long)distinct);

    // ------------------------------------------- single-threaded build
    std::vector<double> times;
    TableST* st = nullptr;
    size_t st_heap = 0;
    for (uint64_t r = 0; r < reps; r++) {
        delete st;
        size_t h0 = heap_bytes();
        double t0 = now();
        st = new TableST(init_cap, 1);
        const auto token = st->register_user();
        for (auto& item : items) insert_item(*st, item, token);
        times.push_back(now() - t0);
        st_heap = heap_bytes() - h0;
    }
    double st_time = median(times);
    check(st->size() == distinct, "single-threaded kache-hash size");
    printf("single-threaded kache-hash: %lu k-mers in %.3f s, %.2f M k-mers/s, %lu distinct k-mers, %.3f B/k-mer "
           "(heap; capacity %lu keys, overflow %lu)\n",
           (unsigned long)kmers, st_time, kmers / st_time / 1e6, (unsigned long)st->size(),
           double(st_heap) / distinct, (unsigned long)st->capacity(), (unsigned long)st->overflow_size());

    // ------------------------------------------- queries on it
    // Streaming: the segments of each sequence, so threads split by sequence
    // like unihash-scale.
    std::vector<std::vector<std::pair<const char*, size_t>>> segs(seqs.size());
    for (size_t i = 0; i < seqs.size(); i++)
        for_each_segment(seqs[i], [&](const char* p, size_t len) { segs[i].push_back({p, len}); });
    printf("queries on the single-threaded table, shared read-only (M queries/s):\n");
    printf("threads | kache positive | random | streaming | CPU check\n");
    const auto lookups = [&](const std::vector<char>& xs) {
        return [&](size_t s, size_t e) {
            uint64_t h = 0;
            Window w;
            for (size_t i = s; i < e; i++) {
                w.init(xs.data() + i * k);
                h += st->find(w) ? 1 : 0;
            }
            return h;
        };
    };
    for (size_t th : threads) {
        if (th > cores) {
            printf("%7lu | skipped: more threads than cores\n", (unsigned long)th);
            continue;
        }
        bool same = true;
        auto [hp, dp] = parallel_count(reps, np, th, lookups(pos), same);
        auto [hn, dn] = parallel_count(reps, nn, th, lookups(neg), same);
        std::atomic<uint64_t> valid{0};
        auto [hs, ds] = parallel_count(reps, seqs.size(), th, [&](size_t s, size_t e) {
            uint64_t h = 0, n = 0;
            for (size_t i = s; i < e; i++)
                for (auto [p, len] : segs[i]) {
                    Window w;
                    w.init(p);
                    for (size_t j = 0;; j++) {
                        h += st->find(w) ? 1 : 0;
                        n++;
                        if (j + k == len) break;
                        w.advance(p[j + k]);
                    }
                }
            valid += n;
            return h;
        }, same);
        check(same, "answers differ between repetitions, " + std::to_string(th) + " threads");
        check(hp == pos_hits, "positive hits, " + std::to_string(th) + " threads");
        check(hn == neg_hits, "random-k-mer hits, " + std::to_string(th) + " threads");
        check(hs == seq_hits && valid == seq_valid * std::max<uint64_t>(reps, 1),
              "streaming hits, " + std::to_string(th) + " threads");
        printf("%7lu | %14.2f | %6.2f | %9.2f | %9.2f\n", (unsigned long)th, np / dp / 1e6, nn / dn / 1e6,
               seq_windows / ds / 1e6, cpu_ratio());
    }
    delete st;

    // ------------------------------------------- concurrent builds
    const size_t mt_cap = capacity ? size_t(distinct) : init_cap;
    printf("concurrent kache-hash (Streaming_Kmer_Hash_Table<k, true>, one token per thread, %s %lu k-mers):\n",
           capacity ? "created for" : "growing from", (unsigned long)mt_cap);
    printf("threads | insert s | M k-mers/s | speed-up vs 1 thread | B/k-mer (heap) | CPU check\n");
    double r1 = 0;
    for (size_t th : threads) {
        if (th > cores) {
            printf("%7lu | skipped: more threads than cores (kache-hash keeps per-core state)\n", (unsigned long)th);
            continue;
        }
        std::vector<double> ts;
        size_t heap = 0, size = 0;
        for (uint64_t r = 0; r < reps; r++) {
            size_t h0 = heap_bytes();
            double t0 = now();
            auto* t = new TableMT(mt_cap, th);
            std::vector<TableMT::Token> tokens;
            for (size_t i = 0; i < th; i++) tokens.push_back(t->register_user());
            std::atomic<size_t> next{0};
            const size_t batch = std::clamp<size_t>(items.size() / (th * 16), 1, BATCH);
            std::vector<std::thread> ws;
            for (size_t i = 0; i < th; i++)
                ws.emplace_back([&, i] {
                    for (;;) {
                        size_t s = next.fetch_add(batch, std::memory_order_relaxed);
                        if (s >= items.size()) break;
                        for (size_t j = s; j < std::min(s + batch, items.size()); j++)
                            insert_item(*t, items[j], tokens[i]);
                    }
                });
            for (auto& w : ws) w.join();
            ts.push_back(now() - t0);
            heap = heap_bytes() - h0;
            size = t->size();
            delete t;
        }
        check(size == distinct, "concurrent kache-hash size, " + std::to_string(th) + " threads");
        double time = median(ts);
        double r = kmers / time / 1e6;
        if (r1 == 0) r1 = r;
        printf("%7lu | %8.3f | %10.2f | %20.2f | %14.3f | %9.2f\n", (unsigned long)th, time, r, r / r1,
               double(heap) / distinct, cpu_ratio());
    }

    printf("\ncheck: %s\n", bad ? "MISMATCH" : "all kache-hash sizes and answers equal unihash's");
    return bad ? 1 : 0;
}
