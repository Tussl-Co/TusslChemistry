// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! The barcode list registry: every list the layouts name, with its file in the whitelist
//! directory, its barcode length and count, the sha256 of the file as fetched, where it comes from
//! and on whose terms. This table is the one source: `data/lists.tsv` is rendered from it (a test
//! refuses a stale copy) and `scripts/fetch_whitelists.sh` reads that file, so the fetch script
//! checks every list against the same count and checksum the code was built with.
//!
//! The checksum is checked when a list is fetched; loading checks each barcode's length, since a
//! library user may load lists of their own.

mod bd;
mod mgi;
mod parse;
pub mod registry;
mod scalebio;
mod seekgene;
mod singleron;
mod split_seq;
mod tenx;

use registry::{Host, ListEntry, REGISTRY};

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum ListId {
    TenxV2,
    TenxV3,
    TenxV4ThreePrime,
    TenxV3FivePrime,
    BdCls1,
    BdCls2,
    BdCls3,
    SplitRound1,
    SplitRound2,
    SplitRound3,
    SingleronV1,
    SingleronBc1,
    SingleronBc2,
    SingleronBc3,
    TenxArcGex,
    TenxFlexV1,
    ParseV3Rounds23,
    ParseV3Round1,
    MgiCSeries,
    BdV2Cls1,
    BdV2Cls2,
    BdV2Cls3,
    ScaleQuantumRt,
    SeekgeneP3cb2,
}

impl ListId {
    pub const ALL: [ListId; 24] = [
        ListId::TenxV2,
        ListId::TenxV3,
        ListId::TenxV4ThreePrime,
        ListId::TenxV3FivePrime,
        ListId::BdCls1,
        ListId::BdCls2,
        ListId::BdCls3,
        ListId::SplitRound1,
        ListId::SplitRound2,
        ListId::SplitRound3,
        ListId::SingleronV1,
        ListId::SingleronBc1,
        ListId::SingleronBc2,
        ListId::SingleronBc3,
        ListId::TenxArcGex,
        ListId::TenxFlexV1,
        ListId::ParseV3Rounds23,
        ListId::ParseV3Round1,
        ListId::MgiCSeries,
        ListId::BdV2Cls1,
        ListId::BdV2Cls2,
        ListId::BdV2Cls3,
        ListId::ScaleQuantumRt,
        ListId::SeekgeneP3cb2,
    ];

    /// How many lists there are; sizes per-list tables.
    pub const COUNT: usize = ListId::ALL.len();

    /// This list's position in `ListId::ALL`, the registry and every per-list table. On the
    /// per-read lookup path, so both properties are proved for every id at compile time
    /// (`ids_index_in_order`) and checked here in checked builds only.
    pub const fn index(self) -> usize {
        let index = self as usize;
        // Hot: checked build only.
        debug_assert!(index < ListId::COUNT, "a list id indexes the per-list tables");
        #[allow(clippy::indexing_slicing)] // earned by ids_index_in_order at compile time
        let listed = ListId::ALL[index] as usize == index;
        // Hot: checked build only.
        debug_assert!(listed, "ListId::ALL is in declaration order");
        index
    }

    /// This list's registry entry.
    pub const fn entry(self) -> &'static ListEntry {
        let index = self.index();
        assert!(index < REGISTRY.len(), "every list id has an entry");
        #[allow(clippy::indexing_slicing)] // earned by the assertion above
        let entry = &REGISTRY[index];
        assert!(entry.id.index() == index, "the registry is in list id order");
        entry
    }

    /// File name in the whitelist directory.
    pub const fn file_name(self) -> &'static str {
        self.entry().file
    }

    /// The length of every barcode on this list: the layouts' barcode segments must match it
    /// (checked at compile time) and so must the list file (checked on loading).
    pub const fn barcode_length(self) -> usize {
        self.entry().length
    }

    pub const fn label(self) -> &'static str {
        self.entry().label
    }

    /// Whether loading may go on without this list. Lists that only a Cell Ranger install holds
    /// are optional: without one, the layouts that need them are reported as not checked.
    pub const fn optional(self) -> bool {
        let entry = self.entry();
        let optional = matches!(entry.source.0, Host::CellRanger);
        // A list with a public source anywhere would not need to be optional.
        assert!(!optional || entry.fallback.is_none(), "a Cell Ranger list has no public fallback");
        optional
    }
}

/// A set of barcode lists, one bit per `ListId`: fixed size, copied by value, built at compile
/// time for each layout.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub struct ListSet(u32);

const _: () = assert!(ListId::COUNT <= u32::BITS as usize, "every list has a bit in a ListSet");

/// Every id indexes the per-list tables, and `ListId::ALL` lists the ids in declaration order,
/// so `index` needs no run-time check in release builds.
const fn ids_index_in_order() -> bool {
    let mut i = 0;
    while i < ListId::COUNT {
        #[allow(clippy::indexing_slicing)] // earned by the loop bound
        let id = ListId::ALL[i];
        assert!((id as usize) < ListId::COUNT, "a list id indexes the per-list tables");
        assert!(id as usize == i, "ListId::ALL is in declaration order");
        i += 1;
    }
    true
}

const _: () = assert!(ids_index_in_order(), "list ids do not index in order");

impl ListSet {
    pub const EMPTY: ListSet = ListSet(0);

    #[must_use]
    pub const fn with(self, id: ListId) -> ListSet {
        let bit = 1 << id.index();
        let set = ListSet(self.0 | bit);
        assert!(set.0 & bit != 0, "the list is in the set");
        assert!(set.0 & self.0 == self.0, "adding a list keeps the others");
        set
    }

    pub const fn contains(self, id: ListId) -> bool {
        let bit = 1 << id.index();
        assert!(self.0 < 1 << ListId::COUNT, "a set holds only list ids");
        assert!(bit < 1 << ListId::COUNT, "a list's bit is inside the set's range");
        self.0 & bit != 0
    }

    /// Whether every list in `self` is also in `other`.
    pub const fn within(self, other: ListSet) -> bool {
        assert!(self.0 < 1 << ListId::COUNT, "a set holds only list ids");
        assert!(other.0 < 1 << ListId::COUNT, "a set holds only list ids");
        self.0 & !other.0 == 0
    }

    pub const fn is_empty(self) -> bool {
        assert!(self.0 < 1 << ListId::COUNT, "a set holds only list ids");
        let empty = self.0 == 0;
        assert!(empty == (self.len() == 0), "an empty set has no lists");
        empty
    }

    pub const fn len(self) -> usize {
        assert!(self.0 < 1 << ListId::COUNT, "a set holds only list ids");
        let len = self.0.count_ones() as usize;
        assert!(len <= ListId::COUNT, "a set holds at most every list");
        len
    }

    /// The lists in the set, in `ListId` order.
    pub fn iter(self) -> impl Iterator<Item = ListId> {
        assert!(self.0 < 1 << ListId::COUNT, "a set holds only list ids");
        assert!(self.len() <= ListId::COUNT, "a set holds at most every list");
        ListId::ALL.into_iter().filter(move |&id| self.contains(id))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The compile-time check that lets `index` skip release assertions, run again here so it
    /// still guards the ids if the `const` assertion that calls it is ever lost.
    #[test]
    fn ids_index_in_order_at_run_time() {
        assert!(ids_index_in_order());
        assert!(ListId::ALL.iter().enumerate().all(|(i, id)| id.index() == i));
    }
}
