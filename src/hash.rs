// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! The Fx hash (Firefox's, as in rustc): one rotate, xor and multiply per word. Fast for the
//! small integer and short byte-string keys this tool counts, and not flooding-resistant, which
//! is fine: keys come from a bounded sample of one input and every count is order-independent
//! (each map's result is sorted, summed or tie-broken on the key), so hash quality only affects
//! speed, never output.

use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hasher};

pub type FxHashMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;
pub type FxHashSet<K> = HashSet<K, BuildHasherDefault<FxHasher>>;

/// A map sized for `capacity` entries up front, so filling it does not reallocate.
pub fn map_with_capacity<K, V>(capacity: usize) -> FxHashMap<K, V> {
    let map = FxHashMap::with_capacity_and_hasher(capacity, BuildHasherDefault::default());
    assert!(map.capacity() >= capacity, "a map asked for {capacity} slots holds fewer");
    assert!(map.is_empty(), "a new map holds nothing");
    map
}

/// Rustc's multiplier: odd, with well-spread bits.
const SEED: u64 = 0x51_7c_c1_b7_27_22_0a_95;
/// Rustc's rotation per word.
const ROTATE: u32 = 5;
/// A short tail's length goes in the top byte, where no tail byte lands.
const TAIL_LENGTH_SHIFT: u32 = 56;
const _: () = assert!(SEED % 2 == 1, "the multiplier must be odd to be invertible");

const _: () = assert!(usize::BITS <= u64::BITS, "a usize widens to a u64 without loss");

#[derive(Default, Clone, Copy)]
pub struct FxHasher {
    hash: u64,
}

impl FxHasher {
    fn add(&mut self, word: u64) {
        self.hash = (self.hash.rotate_left(ROTATE) ^ word).wrapping_mul(SEED);
    }
}

impl Hasher for FxHasher {
    fn write(&mut self, bytes: &[u8]) {
        let mut chunks = bytes.chunks_exact(8);
        for chunk in &mut chunks {
            let Ok(word) = <[u8; 8]>::try_from(chunk) else {
                unreachable!("chunks_exact(8) yields 8 bytes")
            };
            self.add(u64::from_le_bytes(word));
        }
        let tail = chunks.remainder();
        assert!(tail.len() < 8, "chunks_exact leaves under 8 bytes");
        let mut word = [0u8; 8];
        let Some(head) = word.get_mut(..tail.len()) else { unreachable!("tail fits a word") };
        head.copy_from_slice(tail);
        // The length goes in too, so "A" and "A\0" hash apart.
        self.add(u64::from_le_bytes(word) ^ (tail.len() as u64) << TAIL_LENGTH_SHIFT);
    }

    fn write_u32(&mut self, value: u32) {
        self.add(u64::from(value));
    }

    fn write_u64(&mut self, value: u64) {
        self.add(value);
    }

    fn write_usize(&mut self, value: usize) {
        self.add(value as u64);
    }

    fn finish(&self) -> u64 {
        self.hash
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::hash::{BuildHasher, BuildHasherDefault};

    fn hash_of(value: impl std::hash::Hash) -> u64 {
        BuildHasherDefault::<FxHasher>::default().hash_one(value)
    }

    #[test]
    fn distinct_keys_hash_apart_and_equal_keys_alike() {
        assert_eq!(hash_of(42u32), hash_of(42u32));
        assert_ne!(hash_of(42u32), hash_of(43u32));
        assert_ne!(hash_of(&b"ACGT"[..]), hash_of(&b"ACGA"[..]));
        assert_ne!(hash_of(&b"A"[..]), hash_of(&b"A\0"[..]));
        let keys: FxHashSet<u64> = (0..100_000u64).map(hash_of).collect();
        assert_eq!(keys.len(), 100_000, "consecutive integers collide");
    }

    #[test]
    fn maps_count_like_the_standard_map() {
        let mut fx: FxHashMap<&[u8], usize> = map_with_capacity(16);
        let mut std_map: HashMap<&[u8], usize> = HashMap::new();
        for word in [&b"ACGT"[..], b"TTTT", b"ACGT", b"", b"ACGTACGTACGT", b"TTTT", b"ACGT"] {
            *fx.entry(word).or_insert(0) += 1;
            *std_map.entry(word).or_insert(0) += 1;
        }
        assert_eq!(fx.len(), std_map.len());
        assert!(std_map.iter().all(|(k, v)| fx.get(k) == Some(v)));
    }
}
