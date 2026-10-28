//! Counting global allocator: exact live and peak heap bytes of the
//! process (shared by the benchmark binaries).

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};

pub struct Counting;
static CURRENT: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);

#[inline]
fn add(n: usize) {
    let c = CURRENT.fetch_add(n, Relaxed) + n;
    PEAK.fetch_max(c, Relaxed);
}

// SAFETY: every call is forwarded unchanged to the system allocator;
// only byte counters are updated around it.
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        let p = unsafe { System.alloc(l) };
        if !p.is_null() {
            add(l.size());
        }
        p
    }
    unsafe fn alloc_zeroed(&self, l: Layout) -> *mut u8 {
        let p = unsafe { System.alloc_zeroed(l) };
        if !p.is_null() {
            add(l.size());
        }
        p
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        unsafe { System.dealloc(p, l) };
        CURRENT.fetch_sub(l.size(), Relaxed);
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, new: usize) -> *mut u8 {
        let q = unsafe { System.realloc(p, l, new) };
        if !q.is_null() {
            if new >= l.size() {
                add(new - l.size());
            } else {
                CURRENT.fetch_sub(l.size() - new, Relaxed);
            }
        }
        q
    }
}

/// Live heap bytes (sum of requested allocation sizes).
pub fn current() -> usize {
    CURRENT.load(Relaxed)
}
/// Highest live heap since the last `reset_peak`.
pub fn peak() -> usize {
    PEAK.load(Relaxed)
}
pub fn reset_peak() {
    PEAK.store(CURRENT.load(Relaxed), Relaxed);
}
