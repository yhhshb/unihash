// kache-hash side of the unihash comparison (see tools/kache_compare.sh).
//
// Reads a dump directory written by `unihash-bench --dump DIR`:
//   manifest.tsv           k, insertion inputs, query workloads, negatives
//   q<i>.fa                windows of each query workload (one record each)
//   q<i>.inorder.u64       the same k-mers, canonical 2-bit, in order
//   q<i>.shuffled.u64      the same k-mers, shuffled
//   neg.u64                random absent k-mers
// and runs the same measurements on kache-hash (single-threaded, hash set):
//   insert                 every k-mer of the inputs, streamed per segment
//   kmers-in-order         pre-parsed k-mers in sequence order
//   kmers-random-order     the same k-mers shuffled
//   seq-streaming          the window sequences, streamed (kache's native mode)
//   negatives              random absent k-mers
// Answers are checked against the hit counts recorded by unihash-bench.
//
// kache-hash only accepts a k-mer through a `Kmer_Window` built from
// characters, so pre-parsed k-mers are decoded to ASCII *outside* the timed
// region (in batches); `init` (hash + minimizer of the k-mer) is timed, which
// matches unihash recomputing the minimizer of a packed k-mer.

#include "Streaming_Kmer_Hash_Table.hpp"
#include "kseq/kseq.h"

#include <zlib.h>
#include <malloc.h>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

KSEQ_INIT(gzFile, gzread)

#ifndef KH_K
#define KH_K 31
#endif
#ifndef KH_L
#define KH_L 17
#endif

constexpr uint16_t k = KH_K;
constexpr uint16_t l = KH_L;
using Table = kache_hash::Streaming_Kmer_Hash_Table<k, false, void, l>;
using Window = kache_hash::Kmer_Window<k, l>;

static double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// Live heap bytes according to glibc (small + mmapped chunks).
static size_t heap_bytes() {
#if defined(__GLIBC__) && (__GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 33))
    struct mallinfo2 mi = mallinfo2();
    return mi.uordblks + mi.hblkhd;
#else
    struct mallinfo mi = mallinfo();
    return size_t(unsigned(mi.uordblks)) + size_t(unsigned(mi.hblkhd));
#endif
}

// (VmRSS, VmHWM) in bytes.
static std::pair<size_t, size_t> rss() {
    std::ifstream f("/proc/self/status");
    std::string line;
    size_t r = 0, h = 0;
    while (std::getline(f, line)) {
        if (line.rfind("VmRSS:", 0) == 0) r = std::stoull(line.substr(6)) * 1024;
        if (line.rfind("VmHWM:", 0) == 0) h = std::stoull(line.substr(6)) * 1024;
    }
    return {r, h};
}

static std::string rate(uint64_t n, double s) {
    char b[64];
    double r = n / s;
    if (r >= 1e6) snprintf(b, sizeof b, "%.2f M/s", r / 1e6);
    else snprintf(b, sizeof b, "%.1f k/s", r / 1e3);
    return b;
}

static bool is_acgt(char c) { return c == 'A' || c == 'C' || c == 'G' || c == 'T'; }

// Calls f(segment pointer, length) for every maximal ACGT run of length >= k
// (sequence is upper-cased in place).
template <typename F> static void for_each_segment(std::string& s, F f) {
    for (auto& c : s) c = std::toupper(static_cast<unsigned char>(c));
    size_t i = 0;
    while (i < s.size()) {
        while (i < s.size() && !is_acgt(s[i])) i++;
        size_t j = i;
        while (j < s.size() && is_acgt(s[j])) j++;
        if (j - i >= k) f(s.data() + i, j - i);
        i = j;
    }
}

static std::vector<uint64_t> read_u64(const std::string& path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) { std::cerr << "cannot open " << path << "\n"; std::exit(2); }
    size_t n = size_t(f.tellg()) / 8;
    std::vector<uint64_t> v(n);
    f.seekg(0);
    f.read(reinterpret_cast<char*>(v.data()), n * 8);
    return v;
}

static void decode(uint64_t x, char* out) {
    static const char D[4] = {'A', 'C', 'G', 'T'};
    for (int i = 0; i < k; i++) out[i] = D[(x >> (2 * (k - 1 - i))) & 3];
}

// Queries pre-parsed k-mers; decoding to ASCII is done per batch outside the
// timer. Returns (seconds, hits).
static std::pair<double, uint64_t> query_kmers(const Table& t, const std::vector<uint64_t>& xs) {
    constexpr size_t BATCH = 1 << 20;
    std::vector<char> buf(BATCH * k);
    double secs = 0;
    uint64_t hits = 0;
    for (size_t b = 0; b < xs.size(); b += BATCH) {
        size_t n = std::min(BATCH, xs.size() - b);
        for (size_t i = 0; i < n; i++) decode(xs[b + i], buf.data() + i * k);
        double t0 = now();
        Window w;
        for (size_t i = 0; i < n; i++) {
            w.init(buf.data() + i * k);
            hits += t.find(w) ? 1 : 0;
        }
        secs += now() - t0;
    }
    return {secs, hits};
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: kache_bench DUMP_DIR [initial_capacity]\n";
        return 2;
    }
    std::string dir = argv[1];
    size_t init_cap = argc > 2 ? std::stoull(argv[2]) : (1u << 20);

    struct Insert { std::string path; uint64_t skip, max; };
    struct Query { int id; std::string label; uint64_t n, hits; };
    std::vector<Insert> inserts;
    std::vector<Query> queries;
    uint64_t neg_n = 0;
    {
        std::ifstream m(dir + "/manifest.tsv");
        if (!m) { std::cerr << "no manifest in " << dir << "\n"; return 2; }
        std::string line;
        while (std::getline(m, line)) {
            std::vector<std::string> f;
            std::stringstream ss(line);
            std::string x;
            while (std::getline(ss, x, '\t')) f.push_back(x);
            if (f[0] == "k" && std::stoul(f[1]) != k) {
                std::cerr << "dump has k = " << f[1] << " but this binary was built for k = " << k << "\n";
                return 2;
            }
            if (f[0] == "insert") inserts.push_back({f[1], std::stoull(f[2]), std::stoull(f[3])});
            if (f[0] == "query") queries.push_back({std::stoi(f[1]), f[2], std::stoull(f[3]), std::stoull(f[4])});
            if (f[0] == "neg") neg_n = std::stoull(f[1]);
        }
    }

    std::cout << "== kache-hash (k = " << k << ", l = " << l << ", single-threaded hash set, initial capacity "
              << init_cap << ") ==\n";

    // ------------------------------------------------------------ insertion
    // Records are buffered in ~4 MiB batches and each batch inserted under
    // one timer, like unihash-bench.
    auto [rss0, hwm0] = rss();
    size_t heap0 = heap_bytes();
    auto* t = new Table(init_cap, 1);
    const auto token = t->register_user();
    double t_ins = 0;
    uint64_t kmers = 0, records = 0;
    std::vector<std::string> batch;
    size_t batch_bytes = 0;
    auto flush = [&]() {
        double t0 = now();
        for (auto& s : batch) {
            for_each_segment(s, [&](const char* p, size_t len) {
                Window w;
                w.init(p);
                for (size_t i = 0;; i++) {
                    t->insert(w, token);
                    if (i + k == len) break;
                    w.advance(p[i + k]);
                }
                kmers += len - k + 1;
            });
        }
        t_ins += now() - t0;
        batch.clear();
        batch_bytes = 0;
    };
    for (auto& in : inserts) {
        gzFile fp = gzopen(in.path.c_str(), "r");
        if (!fp) { std::cerr << "cannot open " << in.path << "\n"; return 2; }
        kseq_t* seq = kseq_init(fp);
        uint64_t r = 0;
        while (r < in.skip + in.max && kseq_read(seq) >= 0) {
            if (r++ < in.skip) continue;
            batch.emplace_back(seq->seq.s, seq->seq.l);
            batch_bytes += seq->seq.l;
            records++;
            if (batch_bytes >= (4u << 20)) flush();
        }
        kseq_destroy(seq);
        gzclose(fp);
    }
    flush();
    batch.shrink_to_fit();
    size_t heap_table = heap_bytes() - heap0;
    auto [rss1, hwm1] = rss();
    size_t distinct = t->size();

    std::cout << "\n-- insertion --\n";
    std::cout << "records " << records << ", k-mers " << kmers << ", distinct " << distinct << "\n";
    printf("kache   %9.3f s  %s\n", t_ins, rate(kmers, t_ins).c_str());
    std::cout << "\n-- memory --\n";
    printf("kache heap (glibc)         %12zu B  %6.3f B/k-mer  (capacity %zu keys, overflow %zu)\n", heap_table,
           double(heap_table) / distinct, size_t(t->capacity()), t->overflow_size());
    printf("RSS: before table %.1f MiB, after table %.1f MiB (+%.1f MiB), peak %.1f MiB\n", rss0 / 1048576.0,
           rss1 / 1048576.0, (rss1 - rss0) / 1048576.0, hwm1 / 1048576.0);

    // ------------------------------------------------------------ queries
    int bad = 0;
    for (auto& q : queries) {
        std::string base = dir + "/q" + std::to_string(q.id);
        auto in_order = read_u64(base + ".inorder.u64");
        auto shuffled = read_u64(base + ".shuffled.u64");
        auto [s_in, h_in] = query_kmers(*t, in_order);
        auto [s_rnd, h_rnd] = query_kmers(*t, shuffled);

        // Window sequences, streamed: init at each ACGT segment, then advance.
        std::vector<std::string> wins;
        {
            gzFile fp = gzopen((base + ".fa").c_str(), "r");
            kseq_t* seq = kseq_init(fp);
            while (kseq_read(seq) >= 0) wins.emplace_back(seq->seq.s, seq->seq.l);
            kseq_destroy(seq);
            gzclose(fp);
        }
        std::vector<std::pair<const char*, size_t>> segs;
        for (auto& w : wins) for_each_segment(w, [&](const char* p, size_t len) { segs.push_back({p, len}); });
        uint64_t n_str = 0, h_str = 0;
        double t0 = now();
        for (auto [p, len] : segs) {
            Window w;
            w.init(p);
            for (size_t i = 0;; i++) {
                h_str += t->find(w) ? 1 : 0;
                n_str++;
                if (i + k == len) break;
                w.advance(p[i + k]);
            }
        }
        double s_str = now() - t0;

        std::cout << "\n-- queries: " << q.label << " --\n";
        printf("workload: %lu k-mers, expected hits %lu\n", (unsigned long)q.n, (unsigned long)q.hits);
        printf("%-26s %8.3f s  %11s   hits %lu\n", "kmers-in-order", s_in, rate(q.n, s_in).c_str(), (unsigned long)h_in);
        printf("%-26s %8.3f s  %11s   hits %lu\n", "kmers-random-order", s_rnd, rate(q.n, s_rnd).c_str(), (unsigned long)h_rnd);
        printf("%-26s %8.3f s  %11s   hits %lu\n", "seq-streaming", s_str, rate(n_str, s_str).c_str(), (unsigned long)h_str);
        uint64_t mm = (h_in != q.hits) + (h_rnd != q.hits) + (h_str != q.hits) + (n_str != q.n);
        printf("  MISMATCHES vs unihash oracle %lu\n", (unsigned long)mm);
        bad += mm != 0;
    }
    if (neg_n) {
        auto negs = read_u64(dir + "/neg.u64");
        auto [s, h] = query_kmers(*t, negs);
        std::cout << "\n-- negative queries --\n";
        printf("kache   %9.3f s  %s\n", s, rate(negs.size(), s).c_str());
        printf("  FALSE POSITIVES              %lu\n", (unsigned long)h);
        bad += h != 0;
    }
    delete t;
    return bad ? 1 : 0;
}
