// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Per-position base composition of a read file, and the one-line sketch that shows a layout at a
//! glance.

use crate::constants::{
    FIXED_BASE_SHARE_MIN, PROFILE_POSITIONS, SKETCH_LEAD_SHARE_MIN, SKETCH_RANDOM_BITS_MIN,
};
use crate::hash::FxHashMap;
use crate::input::fastq::Reads;

pub struct Profile {
    pub reads: usize,
    pub len_min: usize,
    pub len_mode: usize,
    pub len_max: usize,
    counts: Vec<[u32; BASE_CLASSES]>,
}

/// Counts below this convert to f64 exactly (a 53-bit mantissa).
const F64_EXACT_INTEGERS: usize = 1 << f64::MANTISSA_DIGITS;

/// Base classes counted per position: A, C, G, T, and anything else (N).
pub const BASE_CLASSES: usize = 5;

/// A, C, G, T map to 0..=3; anything else (N, IUPAC) to 4.
pub fn base_index(base: u8) -> usize {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

/// Adds one observation of `base` to a position's A/C/G/T/other counts.
pub fn count_base(slot: &mut [u32; BASE_CLASSES], base: u8) {
    let index = base_index(base);
    #[cfg(debug_assertions)]
    let before = *slot;
    let Some(count) = slot.get_mut(index) else {
        unreachable!("base_index maps every byte into 0..=4, not {index}")
    };
    *count += 1;
    // Hot: checked build only. Exactly the base's own count moved, by one.
    #[cfg(debug_assertions)]
    {
        debug_assert_eq!(slot.get(index), before.get(index).map(|n| n + 1).as_ref());
        debug_assert_eq!(
            slot.iter().map(|&n| u64::from(n)).sum::<u64>(),
            before.iter().map(|&n| u64::from(n)).sum::<u64>() + 1,
            "one observation adds one count"
        );
    }
}

/// `n / d` as a fraction. Counts here are read counts, at most `SAMPLE_READS_MAX` (10^6), far
/// inside f64's exact integer range, so the conversion loses nothing.
#[allow(clippy::cast_precision_loss)]
pub fn ratio(n: usize, d: usize) -> f64 {
    // Exact below 2^53: both counts convert to f64 without rounding.
    assert!(n < F64_EXACT_INTEGERS, "numerator {n} past f64's exact integers");
    assert!(d < F64_EXACT_INTEGERS, "denominator {d} past f64's exact integers");
    let fraction = n as f64 / d.max(1) as f64;
    assert!(fraction >= 0.0, "a ratio of counts is never negative: {fraction}");
    if n <= d {
        assert!(fraction <= 1.0, "a part over its whole is at most 1: {n}/{d} = {fraction}");
    }
    fraction
}

/// Shannon entropy over A/C/G/T at a position, in bits; None when no A/C/G/T was seen.
pub fn entropy(c: &[u32; BASE_CLASSES]) -> Option<f64> {
    let [n_a, n_c, n_g, n_t, _] = *c;
    let total = n_a + n_c + n_g + n_t;
    if total == 0 {
        return None;
    }
    let h: f64 = [n_a, n_c, n_g, n_t]
        .into_iter()
        .filter(|&n| n > 0)
        .map(|n| {
            let p = f64::from(n) / f64::from(total);
            -p * p.log2()
        })
        .sum();
    assert!((0.0..=2.0 + 1e-9).contains(&h), "entropy {h} outside 0..=2 bits for {c:?}");
    let bases_seen = [n_a, n_c, n_g, n_t].iter().filter(|&&n| n > 0).count();
    if bases_seen == 1 {
        assert!(h == 0.0, "a column of one base carries no entropy, not {h} for {c:?}");
    }
    Some(h)
}

/// The most frequent of A/C/G/T and its share. On a tie the later base wins, as `max_by_key` does.
pub fn dominant(c: &[u32; BASE_CLASSES]) -> (u8, f64) {
    let [n_a, n_c, n_g, n_t, _] = *c;
    let (base, top) = [(b'A', n_a), (b'C', n_c), (b'G', n_g), (b'T', n_t)]
        .into_iter()
        .fold((b'A', n_a), |best, (base, n)| if n >= best.1 { (base, n) } else { best });
    assert!(
        [n_a, n_c, n_g, n_t].iter().all(|&n| n <= top),
        "the dominant base {} has the largest count in {c:?}",
        base as char
    );
    let share = f64::from(top) / f64::from((n_a + n_c + n_g + n_t).max(1));
    assert!((0.0..=1.0).contains(&share), "dominant share {share} outside 0..=1");
    (base, share)
}

impl Profile {
    pub fn of(sequences: &Reads) -> Self {
        assert!(!sequences.is_empty(), "profiling an empty sample");
        let mut counts = vec![[0u32; BASE_CLASSES]; PROFILE_POSITIONS];
        let mut lens: FxHashMap<usize, usize> = FxHashMap::default();
        let (mut len_min, mut len_max) = (usize::MAX, 0);
        for sequence in sequences.iter() {
            *lens.entry(sequence.len()).or_insert(0) += 1;
            len_min = len_min.min(sequence.len());
            len_max = len_max.max(sequence.len());
            for (slot, &base) in counts.iter_mut().zip(sequence) {
                count_base(slot, base);
            }
        }
        // The most common length; on a tie the longer one.
        let Some((&len_mode, _)) = lens.iter().max_by_key(|&(len, n)| (*n, *len)) else {
            unreachable!("a non-empty sample has at least one read length")
        };
        assert!(len_min <= len_mode, "shortest read {len_min} longer than the mode {len_mode}");
        assert!(len_mode <= len_max, "the mode {len_mode} longer than the longest read {len_max}");
        assert_eq!(
            lens.values().sum::<usize>(),
            sequences.len(),
            "every read counted once by length"
        );
        Self { reads: sequences.len(), len_min, len_mode, len_max, counts }
    }

    /// Columns where one base exceeds 70% among at least half the reads: the fixed positions a
    /// linker, barcode set or tail leaves. Near zero for cDNA, where reads start anywhere.
    pub fn fixed_columns(&self) -> usize {
        assert_eq!(self.counts.len(), PROFILE_POSITIONS, "one count slot per profiled position");
        let fixed = self
            .counts
            .iter()
            .filter(|c| {
                let [n_a, n_c, n_g, n_t, n_other] = **c;
                let total = usize::try_from(n_a + n_c + n_g + n_t).unwrap_or(usize::MAX);
                let observed = total + usize::try_from(n_other).unwrap_or(usize::MAX);
                assert!(observed <= self.reads, "a column saw {observed} bases from fewer reads");
                total * 2 >= self.reads && dominant(c).1 > FIXED_BASE_SHARE_MIN
            })
            .count();
        assert!(fixed <= PROFILE_POSITIONS, "{fixed} fixed columns out of {PROFILE_POSITIONS}");
        fixed
    }

    /// One character per position: '.' random (> 1.8 bits), an upper-case base when one base is
    /// over 70%, lower-case when it merely leads, '-' in between. Shows the layout at a glance.
    pub fn sketch(&self) -> String {
        assert!(self.reads > 0, "a profile covers at least one read");
        let sketch: String = self
            .counts
            .iter()
            .map_while(|c| {
                let h = entropy(c)?;
                let (base, frac) = dominant(c);
                Some(if h > SKETCH_RANDOM_BITS_MIN {
                    '.'
                } else if frac > FIXED_BASE_SHARE_MIN {
                    base as char
                } else if frac > SKETCH_LEAD_SHARE_MIN {
                    base.to_ascii_lowercase() as char
                } else {
                    '-'
                })
            })
            .collect();
        assert!(sketch.len() <= PROFILE_POSITIONS, "one character per profiled position at most");
        assert!(
            sketch
                .chars()
                .all(|c| matches!(c, '.' | '-' | 'A' | 'C' | 'G' | 'T' | 'a' | 'c' | 'g' | 't')),
            "a sketch is '.', '-' or a base: {sketch}"
        );
        sketch
    }
}
