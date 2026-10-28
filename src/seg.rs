//! Segmented array of `AtomicU64` with stable addresses, for the concurrent
//! table.
//!
//! Elements live in segments that are allocated (zeroed) on first write and
//! never move, so readers need no lock; reading an element that was never
//! written returns 0. Segment sizes grow geometrically in groups: group `g`
//! has `2^group_bits` segments of `2^(seg_bits + g)` elements, so the
//! directory stays small (a few hundred cells), an index is located with a
//! few shifts, and the unused part of the last segment is at most
//! `2^-group_bits` of the allocated size.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

type Segment = Box<[AtomicU64]>;

/// `SEG`: log2 of the first segments' length; `GROUP`: log2 of the number
/// of segments of each length.
pub(crate) struct SegArray<const SEG: u32, const GROUP: u32> {
    dir: Box<[OnceLock<Segment>]>,
    capacity: u64,
}

impl<const SEG: u32, const GROUP: u32> SegArray<SEG, GROUP> {
    /// An array of at least `max` elements.
    pub(crate) fn new(max: u64) -> Self {
        let (seg_bits, group_bits) = (SEG, GROUP);
        let unit = 1u64 << (seg_bits + group_bits);
        let mut groups = 0u32;
        while unit * ((1u64 << groups) - 1) < max {
            groups += 1;
        }
        Self {
            dir: (0..(groups as usize) << group_bits)
                .map(|_| OnceLock::new())
                .collect(),
            capacity: unit * ((1u64 << groups) - 1),
        }
    }

    /// Elements the array can address.
    pub(crate) fn capacity(&self) -> u64 {
        self.capacity
    }

    /// (directory cell, offset in segment, segment length) of element `i`.
    #[inline(always)]
    fn locate(&self, i: u64) -> (usize, usize, usize) {
        let u = i + (1u64 << (SEG + GROUP));
        let top = 63 - u.leading_zeros();
        let shift = top - GROUP; // = SEG + group
        let g = shift - SEG;
        let s = (u >> shift) - (1 << GROUP);
        (
            ((g as usize) << GROUP) + s as usize,
            (u & ((1u64 << shift) - 1)) as usize,
            1usize << shift,
        )
    }

    /// Element `i` if its segment exists.
    #[inline(always)]
    pub(crate) fn get(&self, i: u64) -> Option<&AtomicU64> {
        let (d, o, _) = self.locate(i);
        self.dir.get(d)?.get().map(|s| &s[o])
    }

    #[inline(always)]
    pub(crate) fn load(&self, i: u64, order: Ordering) -> u64 {
        self.get(i).map_or(0, |a| a.load(order))
    }

    /// Elements `i` and `i + 1` (relaxed), locating the segment once when
    /// both are in it.
    #[inline(always)]
    pub(crate) fn load2(&self, i: u64) -> (u64, u64) {
        let (d, o, _) = self.locate(i);
        match self.dir.get(d).and_then(OnceLock::get) {
            Some(seg) => match seg.get(o + 1) {
                Some(b) => (seg[o].load(Ordering::Relaxed), b.load(Ordering::Relaxed)),
                None => (
                    seg[o].load(Ordering::Relaxed),
                    self.load(i + 1, Ordering::Relaxed),
                ),
            },
            None => (0, 0),
        }
    }

    /// The segment holding element `i` (allocated if needed) and the index
    /// of its first element.
    pub(crate) fn segment(&self, i: u64) -> (u64, &[AtomicU64]) {
        let (d, o, n) = self.locate(i);
        assert!(d < self.dir.len(), "segmented array full");
        let seg = self.dir[d].get_or_init(|| (0..n).map(|_| AtomicU64::new(0)).collect());
        (i - o as u64, seg)
    }

    /// Element `i`, allocating its segment if needed.
    #[inline]
    pub(crate) fn slot(&self, i: u64) -> &AtomicU64 {
        let (first, seg) = self.segment(i);
        &seg[(i - first) as usize]
    }

    /// Copy elements `0..n` into a vector (unwritten ones are 0).
    pub(crate) fn to_vec(&self, n: u64) -> Vec<u64> {
        let mut v = Vec::with_capacity(n as usize);
        let mut i = 0;
        while i < n {
            let (d, o, len) = self.locate(i);
            let take = ((len - o) as u64).min(n - i) as usize;
            match self.dir.get(d).and_then(OnceLock::get) {
                Some(seg) => v.extend(seg[o..o + take].iter().map(|a| a.load(Ordering::Relaxed))),
                None => v.resize(v.len() + take, 0),
            }
            i += take as u64;
        }
        v
    }

    /// Fill elements `0..src.len()` from a slice (single-threaded setup).
    pub(crate) fn fill_from(&self, src: &[u64]) {
        let mut i = 0;
        while i < src.len() {
            let (first, seg) = self.segment(i as u64);
            let o = i - first as usize;
            let take = (seg.len() - o).min(src.len() - i);
            for (a, &x) in seg[o..o + take].iter().zip(&src[i..i + take]) {
                a.store(x, Ordering::Relaxed);
            }
            i += take;
        }
    }

    /// Heap bytes: directory plus allocated segments.
    pub(crate) fn heap_bytes(&self) -> usize {
        self.dir.len() * std::mem::size_of::<OnceLock<Segment>>()
            + self
                .dir
                .iter()
                .filter_map(OnceLock::get)
                .map(|s| s.len() * 8)
                .sum::<usize>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn geometric_segments() {
        let a = SegArray::<2, 1>::new(1000); // segments 4,4 | 8,8 | 16,16 | ...
        assert!(a.capacity() >= 1000);
        let mut start = 0u64;
        let mut expect = Vec::new();
        for g in 0..a.dir.len() / 2 {
            for _ in 0..2 {
                expect.push((start, 4usize << g));
                start += 4 << g;
            }
        }
        for (d, &(s, n)) in expect.iter().enumerate() {
            for i in [s, s + n as u64 - 1] {
                assert_eq!(a.locate(i), (d, (i - s) as usize, n), "i={i}");
            }
        }
    }

    #[test]
    fn lazily_allocated_and_stable() {
        let a = SegArray::<3, 2>::new(128);
        assert!(a.capacity() >= 128);
        assert_eq!(a.load(77, Ordering::Relaxed), 0);
        assert!(a.get(77).is_none());
        a.slot(77).store(5, Ordering::Relaxed);
        let p = a.get(77).unwrap() as *const AtomicU64;
        a.slot(127).store(9, Ordering::Relaxed);
        a.slot(0).store(1, Ordering::Relaxed);
        assert_eq!(a.get(77).unwrap() as *const AtomicU64, p);
        let v = a.to_vec(128);
        assert_eq!((v[0], v[77], v[127], v[1]), (1, 5, 9, 0));
        assert_eq!(v.iter().filter(|&&x| x != 0).count(), 3);
        let b = SegArray::<3, 2>::new(128);
        b.fill_from(&v);
        assert_eq!(b.to_vec(128), v);
        let (first, seg) = b.segment(77);
        assert!(first <= 77 && 77 < first + seg.len() as u64);
        assert_eq!(seg[(77 - first) as usize].load(Ordering::Relaxed), 5);
        assert_eq!(b.load2(76), (0, 5));
        // across a segment boundary
        let (f2, s2) = b.segment(first - 1);
        assert_eq!(f2 + s2.len() as u64, first);
        b.slot(first - 1).store(3, Ordering::Relaxed);
        b.slot(first).store(4, Ordering::Relaxed);
        assert_eq!(b.load2(first - 1), (3, 4));
    }
}
