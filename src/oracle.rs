//! Hash-set oracle used by tests and benchmarks.

use crate::codec::{canonicalize, kmer_mask, KmerStream};
use std::collections::HashSet;

#[derive(Clone, Debug)]
pub struct Oracle {
    k: usize,
    set: HashSet<u64>,
    stream: KmerStream,
}

impl Oracle {
    pub fn new(k: usize) -> Self {
        Self {
            k,
            set: HashSet::new(),
            stream: KmerStream::new(k),
        }
    }

    /// Insert every canonical k-mer of `read` (split at non-ACGT characters,
    /// exactly like [`crate::KmerTable::insert_read`]).
    pub fn insert_read(&mut self, read: &[u8]) {
        self.insert_chunk(read);
        self.end_read();
    }

    /// Streaming counterpart of [`crate::KmerTable::insert_chunk`].
    pub fn insert_chunk(&mut self, chunk: &[u8]) {
        let set = &mut self.set;
        self.stream.push(chunk, |x| {
            set.insert(x);
        });
    }

    pub fn end_read(&mut self) {
        self.stream.end();
    }

    /// Remove every canonical k-mer of `read`; returns how many were present.
    pub fn delete_read(&mut self, read: &[u8]) -> u64 {
        let mut s = KmerStream::new(self.k);
        let mut n = 0;
        let set = &mut self.set;
        s.push(read, |x| n += set.remove(&x) as u64);
        n
    }

    /// Remove one k-mer (either orientation); returns whether it was present.
    pub fn delete(&mut self, kmer: u64) -> bool {
        kmer & !kmer_mask(self.k) == 0 && self.set.remove(&canonicalize(kmer, self.k))
    }

    pub fn contains(&self, kmer: u64) -> bool {
        kmer & !kmer_mask(self.k) == 0 && self.set.contains(&canonicalize(kmer, self.k))
    }

    pub fn len(&self) -> usize {
        self.set.len()
    }

    pub fn is_empty(&self) -> bool {
        self.set.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &u64> + '_ {
        self.set.iter()
    }

    /// Approximate heap bytes of the std `HashSet<u64>` (SwissTable): one
    /// 8-byte slot plus one control byte per bucket, plus 16 trailing
    /// control bytes. Buckets are the power of two ≥ capacity·8/7.
    pub fn heap_bytes(&self) -> usize {
        let cap = self.set.capacity();
        if cap == 0 {
            return 0;
        }
        let buckets = if cap < 8 {
            if cap < 4 {
                4
            } else {
                8
            }
        } else {
            (cap * 8 / 7).next_power_of_two()
        };
        buckets * 9 + 16
    }

    /// Bytes per k-mer of an ideal, perfectly sized open-addressing table of
    /// `u64` at load factor 7/8 — a lower bound for any plain hash set.
    pub fn ideal_bytes(&self) -> usize {
        (self.set.len() as f64 * 8.0 * 8.0 / 7.0) as usize
    }

    pub fn shrink_to_fit(&mut self) {
        self.set.shrink_to_fit();
    }
}
