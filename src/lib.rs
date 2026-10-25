#![forbid(unsafe_code)]
//! # unihash — exact k-mer presence table over piece storage (MVP)
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

pub mod codec;
pub mod edit;
pub mod fastx;
pub mod index;
pub mod minimizer;
pub mod oracle;
pub mod piece_store;
pub mod rng;
pub mod stream_query;
pub mod table;

pub use edit::{EditStats, Fragmentation};
pub use oracle::Oracle;
pub use stream_query::{StreamStats, StreamingQuery};
pub use table::{
    CandidateScan, Config, ConfigError, Hit, InsertStats, KmerError, KmerTable, QueryStats,
    QueryTrace, TableStats,
};
