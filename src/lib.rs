#![forbid(unsafe_code)]
//! # unihash — exact k-mer presence table over piece storage
//!
//! Stores a set of canonical k-mers as sequence *pieces* (maximal runs of
//! previously absent k-mers taken from inserted reads) and answers exact
//! membership queries through a minimizer index.
//!
//! ```
//! use unihash::{Config, KmerTable};
//!
//! let mut t = KmerTable::new(Config { k: 5, l: 3, ..Config::default() }).unwrap();
//! t.insert_read(b"ACGTACGGTCA");
//! assert!(t.contains_seq(b"CGTAC").unwrap());
//! assert!(t.contains_seq(b"GTACG").unwrap());   // reverse complement of CGTAC
//! assert!(!t.contains_seq(b"AAAAA").unwrap());
//! ```
//!
//! ## Correctness contract
//! After any sequence of `insert_read` and `delete_read` / `delete` calls,
//! `contains(x)` returns `true` iff the canonical form of `x` was inserted as
//! a k-mer (made only of `ACGTacgt` characters) and not deleted since. There
//! are no false positives and no false negatives; answers do not depend on
//! the order of operations, nor on `defrag` calls.
//!
//! ## Several threads
//! [`ConcurrentKmerTable`] accepts insertions from many threads (one
//! [`Writer`] each) while others delete and query, with the same guarantees
//! for every insertion and deletion that has completed, also through
//! streaming queries ([`ConcurrentKmerTable::query_stream`]).
//! Defragmentation takes `&mut` and is done in place. A `&KmerTable` can be
//! queried from many threads without any locking.
//!
//! ```
//! use unihash::{Config, ConcurrentKmerTable};
//!
//! let t = ConcurrentKmerTable::new(Config { k: 5, l: 3, ..Config::default() }).unwrap();
//! std::thread::scope(|s| {
//!     for read in [&b"ACGTACGGTCA"[..], b"TTGACCGTACGT"] {
//!         let t = &t;
//!         s.spawn(move || t.writer().insert_read(read));
//!     }
//! });
//! assert!(t.contains(0b01_10_11_00_01)); // CGTAC
//! let table = t.into_table();
//! assert_eq!(table.redundant_kmers(), 0);
//! ```

pub mod codec;
pub mod concurrent;
pub mod edit;
pub mod fastx;
pub mod index;
pub mod minimizer;
pub mod oracle;
pub mod piece_store;
pub mod rng;
mod seg;
pub mod stream_query;
pub mod table;

pub use concurrent::{
    ConcurrentKmerTable, ConcurrentStats, ConcurrentStreamingQuery, ContentionStats, ReadSnapshot,
    Writer,
};
pub use edit::{EditStats, Fragmentation};
pub use oracle::Oracle;
pub use stream_query::{StreamStats, StreamingQuery};
pub use table::{
    CandidateScan, Config, ConfigError, Hit, InsertStats, KmerError, KmerTable, QueryStats,
    QueryTrace, RunIndexing, TableStats,
};
