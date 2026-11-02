// Stress test of the concurrent kache-hash table (see kache-hash-lock.patch):
// THREADS threads insert the same k-mers of a random sequence of LEN bases,
// in different orders, into a table created for the distinct k-mer count,
// ROUNDS times; each round checks size() and that every k-mer is found.
//
//   bash tools/kache_compare.sh --setup        # builds external/bin/kache_lock_stress
//   external/bin/kache_lock_stress 100000 4000 2
//   KACHE_PATCH=0 bash tools/kache_compare.sh --setup   # upstream code, to see the bug
//
// Exit status 1 if any round was wrong. THREADS must not exceed the
// hardware threads (kache-hash keeps one counter per hardware thread).
#include "Streaming_Kmer_Hash_Table.hpp"
#include <atomic>
#include <cstdio>
#include <random>
#include <string>
#include <thread>
#include <vector>
constexpr uint16_t k = 31, l = 17;
using TableST = kache_hash::Streaming_Kmer_Hash_Table<k, false, void, l>;
using TableMT = kache_hash::Streaming_Kmer_Hash_Table<k, true, void, l>;
using Window = kache_hash::Kmer_Window<k, l>;
template <typename T, typename Tok> void ins(T& t, const std::string& s, const Tok& tok) {
    Window w; w.init(s.data());
    for (size_t i = 0;; i++) { t.insert(w, tok); if (i + k == s.size()) break; w.advance(s[i + k]); }
}
int main(int argc, char** argv) {
    size_t len = std::stoull(argv[1]), rounds = std::stoull(argv[2]), th = std::stoull(argv[3]);
    std::mt19937_64 rng(7);
    std::string g(len, 'A');
    for (auto& c : g) c = "ACGT"[rng() & 3];
    std::vector<std::string> items;
    for (size_t p = 0; p + k < len; p += 500) items.push_back(g.substr(p, std::min<size_t>(500 + k - 1, len - p)));
    TableST st(1 << 16); auto t0 = st.register_user();
    for (auto& it : items) ins(st, it, t0);
    const size_t distinct = st.size();
    size_t bad = 0, lost = 0;
    for (size_t r = 0; r < rounds; r++) {
        TableMT mt(distinct, th);
        std::vector<TableMT::Token> tok; for (size_t i = 0; i < th; i++) tok.push_back(mt.register_user());
        std::vector<std::thread> ws;
        for (size_t i = 0; i < th; i++) ws.emplace_back([&, i] {
            std::vector<size_t> ord(items.size()); for (size_t j = 0; j < ord.size(); j++) ord[j] = j;
            std::shuffle(ord.begin(), ord.end(), std::mt19937_64(r * 31 + i));
            for (size_t j : ord) ins(mt, items[j], tok[i]);
        });
        for (auto& w : ws) w.join();
        size_t miss = 0;
        for (auto& it : items) { Window w; w.init(it.data()); for (size_t i = 0;; i++) { if (!mt.find(w)) miss++; if (i + k == it.size()) break; w.advance(it[i + k]); } }
        if (mt.size() != distinct || miss) { bad++; lost += miss; if (bad <= 5) printf("round %zu: size %zu (expected %zu), %zu k-mer lookups missed\n", r, mt.size(), distinct, miss); }
    }
    printf("%zu of %zu rounds wrong (%zu missed lookups)\n", bad, rounds, lost);
    return bad != 0;
}
