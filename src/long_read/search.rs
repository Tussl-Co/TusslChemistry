// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Approximate search: where a short pattern (an adapter, a primer) best matches inside a read,
//! counting substitutions, insertions and deletions alike. Long reads carry indel errors, so an
//! exact or substitution-only match would miss most adapters.
//!
//! Myers' bit-parallel algorithm (J. ACM 1999): one column of the edit-distance matrix per text
//! base, held as bit vectors over the pattern, so a pattern of up to 64 bases costs a handful of
//! word operations per base. Semi-global: the pattern must match whole, anywhere in the text, so
//! the text's start and end are free. Tested against the textbook dynamic programme below.

/// Longest pattern a search handles: one bit per pattern base in a u64.
pub const PATTERN_BASES_MAX: usize = 64;

/// The best place a pattern matches: its edit distance, and the text index just past the match.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Hit {
    pub distance: usize,
    pub end: usize,
}

/// A pattern prepared for searching: for each base A, C, G and T, a bit per pattern position where
/// it occurs. A text base other than A/C/G/T matches no pattern position.
pub struct Pattern {
    equal: [u64; 4],
    length: usize,
}

fn base_code(base: u8) -> Option<usize> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

impl Pattern {
    /// Prepares `bases` (A/C/G/T only, 1 to 64 of them) for searching.
    pub fn new(bases: &[u8]) -> Self {
        assert!(!bases.is_empty(), "an empty pattern matches everywhere");
        assert!(bases.len() <= PATTERN_BASES_MAX, "a pattern of {} bases", bases.len());
        let mut equal = [0u64; 4];
        for (index, &base) in bases.iter().enumerate() {
            let Some(code) = base_code(base) else {
                unreachable!("pattern base {base:#04x} is not A, C, G or T")
            };
            let Some(slot) = equal.get_mut(code) else { unreachable!("base codes are 0..4") };
            *slot |= 1 << index;
        }
        let total: u32 = equal.iter().map(|bits| bits.count_ones()).sum();
        assert_eq!(total as usize, bases.len(), "every pattern position has exactly one base");
        Self { equal, length: bases.len() }
    }

    /// The best match in `text`: the lowest edit distance, ties going to the earliest end. None
    /// when even the best is over `distance_max` edits.
    pub fn search(&self, text: &[u8], distance_max: usize) -> Option<Hit> {
        assert!(
            distance_max < self.length,
            "a pattern within {distance_max} edits matches anything"
        );
        let high = 1u64 << (self.length - 1);
        let mask =
            if self.length == PATTERN_BASES_MAX { u64::MAX } else { (1u64 << self.length) - 1 };
        // Column 0: distance i to match the first i pattern bases against nothing.
        let (mut positive, mut negative, mut score) = (mask, 0u64, self.length);
        let mut best: Option<Hit> = None;
        for (index, &base) in text.iter().enumerate() {
            let equal =
                base_code(base).map_or(0, |code| self.equal.get(code).copied().unwrap_or(0));
            let vertical = equal | negative;
            let horizontal = (((equal & positive).wrapping_add(positive)) ^ positive) | equal;
            let mut plus = negative | !(horizontal | positive);
            let mut minus = positive & horizontal;
            if plus & high != 0 {
                score += 1;
            } else if minus & high != 0 {
                score -= 1;
            }
            // Semi-global: the top row is all zeros (the match may start anywhere), so nothing
            // shifts in at the bottom.
            plus = (plus << 1) & mask;
            minus = (minus << 1) & mask;
            positive = (minus | !(vertical | plus)) & mask;
            negative = (plus & vertical) & mask;
            // Hot: checked build only.
            debug_assert!(score <= self.length + index + 1, "a column score past its bound");
            if score <= distance_max && best.is_none_or(|hit| score < hit.distance) {
                best = Some(Hit { distance: score, end: index + 1 });
            }
        }
        if let Some(hit) = best {
            assert!(hit.distance <= distance_max, "a reported match is within the bound");
            assert!(hit.end <= text.len(), "a match ends inside the text");
        }
        best
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The textbook semi-global dynamic programme: rows over the pattern, columns over the text,
    /// a free start anywhere in the text.
    fn reference(pattern: &[u8], text: &[u8], distance_max: usize) -> Option<Hit> {
        let mut column: Vec<usize> = (0..=pattern.len()).collect();
        let mut best: Option<Hit> = None;
        for (index, &t) in text.iter().enumerate() {
            let mut previous_diagonal = column[0];
            column[0] = 0;
            for row in 1..=pattern.len() {
                let same = base_code(t).is_some() && base_code(t) == base_code(pattern[row - 1]);
                let substitute = previous_diagonal + usize::from(!same);
                let value = substitute.min(column[row] + 1).min(column[row - 1] + 1);
                previous_diagonal = column[row];
                column[row] = value;
            }
            let score = column[pattern.len()];
            if score <= distance_max && best.is_none_or(|hit| score < hit.distance) {
                best = Some(Hit { distance: score, end: index + 1 });
            }
        }
        best
    }

    fn bases(state: &mut u64, length: usize, alphabet: &[u8]) -> Vec<u8> {
        (0..length)
            .map(|_| {
                *state = state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                alphabet[usize::try_from(*state >> 33).unwrap() % alphabet.len()]
            })
            .collect()
    }

    #[test]
    fn matches_the_textbook_programme() {
        let mut state = 17;
        for pattern_length in [1, 2, 7, 16, 22, 31, 32, 33, 63, 64] {
            for _ in 0..300 {
                let pattern = bases(&mut state, pattern_length, b"ACGT");
                let text_length = usize::try_from(state >> 58).unwrap() * 3;
                let mut text = bases(&mut state, text_length, b"ACGTN");
                // Plant a mutated copy of the pattern half the time, so close matches occur.
                if state & 1 == 1 && text.len() > pattern.len() {
                    let at = usize::try_from(state >> 40).unwrap() % (text.len() - pattern.len());
                    text[at..at + pattern.len()].copy_from_slice(&pattern);
                    let flip = at + usize::try_from(state >> 20).unwrap() % pattern.len();
                    text[flip] = b'N';
                }
                let distance_max = pattern_length.saturating_sub(1).min(5);
                let prepared = Pattern::new(&pattern);
                assert_eq!(
                    prepared.search(&text, distance_max),
                    reference(&pattern, &text, distance_max),
                    "pattern {} text {}",
                    String::from_utf8_lossy(&pattern),
                    String::from_utf8_lossy(&text)
                );
            }
        }
    }

    #[test]
    fn finds_an_adapter_through_an_insertion_and_a_deletion() {
        let adapter = b"CTACACGACGCTCTTCCGATCT";
        let prepared = Pattern::new(adapter);
        // One base inserted, one deleted, inside a read.
        let read = b"GGGGGCTACACGACGTCTCTTCGATCTAAACCCGGGTTT";
        let hit = prepared.search(read, 3).unwrap();
        assert_eq!(hit.distance, 2);
        assert_eq!(&read[hit.end - 5..hit.end], b"GATCT");
        assert_eq!(prepared.search(b"ACGTACGTACGTACGT", 3), None);
    }
}
