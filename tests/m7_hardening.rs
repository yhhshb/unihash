//! Input validation and error messages.

use unihash::{Config, ConfigError, KmerError, KmerTable};

fn err(c: Config) -> ConfigError {
    KmerTable::new(c).unwrap_err()
}

#[test]
fn config_bounds() {
    let d = Config::default();
    assert_eq!(
        err(Config {
            k: 1,
            l: 1,
            ..d.clone()
        }),
        ConfigError::KOutOfRange(1)
    );
    assert_eq!(
        err(Config { k: 33, ..d.clone() }),
        ConfigError::KOutOfRange(33)
    );
    assert_eq!(
        err(Config {
            k: 31,
            l: 0,
            ..d.clone()
        }),
        ConfigError::LOutOfRange { l: 0, k: 31 }
    );
    assert_eq!(
        err(Config {
            k: 21,
            l: 22,
            ..d.clone()
        }),
        ConfigError::LOutOfRange { l: 22, k: 21 }
    );
    assert_eq!(
        err(Config {
            bucket_bits: 31,
            ..d.clone()
        }),
        ConfigError::BucketBitsTooLarge(31)
    );
    assert_eq!(
        err(Config {
            max_piece_len: 30,
            ..d.clone()
        }),
        ConfigError::MaxPieceLenTooSmall {
            max_piece_len: 30,
            k: 31
        }
    );
    // boundary values are accepted
    for (k, l) in [(2, 1), (2, 2), (32, 32), (32, 1)] {
        KmerTable::new(Config {
            k,
            l,
            bucket_bits: 0,
            max_piece_len: k as u32,
            ..d.clone()
        })
        .unwrap();
    }
    KmerTable::new(Config {
        max_piece_len: unihash::piece_store::MAX_PIECE_LEN,
        ..d.clone()
    })
    .unwrap();
    assert_eq!(
        err(Config {
            max_piece_len: unihash::piece_store::MAX_PIECE_LEN + 1,
            ..d
        }),
        ConfigError::MaxPieceLenTooLarge(unihash::piece_store::MAX_PIECE_LEN + 1)
    );
}

#[test]
fn config_error_messages_are_descriptive() {
    let m = ConfigError::KOutOfRange(40).to_string();
    assert!(m.contains("k = 40") && m.contains("2 ≤ k ≤ 32"), "{m}");
    let m = ConfigError::LOutOfRange { l: 25, k: 21 }.to_string();
    assert!(m.contains("ℓ = 25") && m.contains("k = 21"), "{m}");
    let m = ConfigError::MaxPieceLenTooSmall {
        max_piece_len: 10,
        k: 31,
    }
    .to_string();
    assert!(m.contains("max_piece_len = 10"), "{m}");
}

#[test]
fn string_query_validation() {
    let mut t = KmerTable::new(Config {
        k: 5,
        l: 3,
        bucket_bits: 4,
        ..Config::default()
    })
    .unwrap();
    t.insert_read(b"ACGTACGT");
    assert_eq!(
        t.contains_seq(b"ACGT"),
        Err(KmerError::WrongLength {
            expected: 5,
            got: 4
        })
    );
    assert_eq!(
        t.contains_seq(b"ACNTA"),
        Err(KmerError::InvalidBase { pos: 2, byte: b'N' })
    );
    assert_eq!(t.contains_seq(b"acgta"), Ok(true));
    assert_eq!(t.contains_seq(b"TACGT"), Ok(true)); // rc of ACGTA
}

#[test]
fn packed_query_with_high_bits_is_rejected() {
    let mut t = KmerTable::new(Config {
        k: 5,
        l: 3,
        bucket_bits: 4,
        ..Config::default()
    })
    .unwrap();
    t.insert_read(b"AAAAA");
    assert!(t.contains(0));
    assert!(!t.contains(1 << 10)); // bit above 2k
    assert!(!t.contains(u64::MAX));
}

#[test]
fn fallback_layer_configuration() {
    let d = Config::default();
    // automatic ℓ2 = min(k, ℓ + 8)
    assert_eq!(d.fallback_l(), Some(25));
    assert_eq!(
        Config {
            k: 21,
            l: 17,
            ..d.clone()
        }
        .fallback_l(),
        Some(21)
    );
    // routing disabled by threshold 0, or when no ℓ2 > ℓ fits in k
    assert_eq!(
        Config {
            hot_threshold: 0,
            ..d.clone()
        }
        .fallback_l(),
        None
    );
    assert_eq!(
        Config {
            k: 17,
            l: 17,
            ..d.clone()
        }
        .fallback_l(),
        None
    );
    // explicit ℓ2 must satisfy ℓ < ℓ2 ≤ k
    for l2 in [17usize, 10, 32] {
        assert_eq!(
            err(Config { l2, ..d.clone() }),
            ConfigError::L2OutOfRange { l2, l: 17, k: 31 }
        );
    }
    let m = ConfigError::L2OutOfRange {
        l2: 10,
        l: 17,
        k: 31,
    }
    .to_string();
    assert!(m.contains("ℓ2 = 10") && m.contains("ℓ < ℓ2 ≤ k"), "{m}");
    KmerTable::new(Config { l2: 31, ..d }).unwrap();
}
