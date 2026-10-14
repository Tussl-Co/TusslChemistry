// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! `simpleaf quant`: the layout's own geometry string, its permit list (or knee filtering when
//! the kit has no list), and the expected orientation from the measured cDNA strand. `-1` is the
//! file holding the layout's read 1 and `-2` the other, whichever way they were given.

use super::{
    Arg, EmitRefusal, Emitted, Matched, Pipeline, Placeholder, cdna, emitted, path_text,
    single_barcode, strand, unsplit, unsupported,
};
use crate::Strand;
use crate::layouts::Mate;

const INDEX: Placeholder =
    Placeholder { name: "SIMPLEAF_INDEX", meaning: "a simpleaf index directory (simpleaf index)" };

pub(super) fn emit(m: &Matched<'_>) -> Result<Emitted, EmitRefusal> {
    let Some((_, listed, _)) = single_barcode(m.structure) else {
        return Err(unsupported(
            Pipeline::Simpleaf,
            m,
            "tchem writes simpleaf commands for one barcode and a UMI from the read's first \
             base so far",
        ));
    };
    let cdna = cdna(m)?;
    let (read1, read2) = match m.structure.read {
        Mate::R1 => (m.barcode_file, cdna),
        Mate::R2 => (cdna, m.barcode_file),
    };
    let orientation = match strand(m)? {
        Strand::Forward => "fw",
        Strand::Reverse => "rc",
    };
    let text = |s: &str| Arg::Text(s.to_owned());
    let mut argv = vec![
        text("simpleaf"),
        text("quant"),
        text("--chemistry"),
        text(m.layout.geometry),
        text("--index"),
        Arg::Var(INDEX),
        text("--reads1"),
        Arg::Text(unsplit(path_text(read1)?, "--reads1")?),
        text("--reads2"),
        Arg::Text(unsplit(path_text(read2)?, "--reads2")?),
        text("--expected-ori"),
        text(orientation),
    ];
    match (listed, m.layout.barcode_lists.as_slice()) {
        (true, [list]) => argv.extend([text("--unfiltered-pl"), Arg::Text(path_text(list)?)]),
        (false, []) => argv.push(text("--knee")),
        _ => unreachable!("{}: one listed barcode has one list", m.layout.name),
    }
    argv.extend([text("--resolution"), text("cr-like"), text("--output"), text("simpleaf_quant")]);
    let mut e = emitted(Pipeline::Simpleaf, m, argv);
    if crate::layouts::nuclei(m.structure) {
        e.notes.push(
            "nuclei (Multiome): with a spliceu or splici index (simpleaf index's own), sum the \
             spliced, unspliced and ambiguous counts, as Cell Ranger ARC counts introns"
                .into(),
        );
    }
    assert!(e.placeholders == [INDEX], "simpleaf needs the index alone");
    assert!(read1 != read2, "two different files");
    Ok(e)
}

#[cfg(test)]
mod tests {
    use crate::layouts::{LAYOUTS, Layout, Mate, Segment};

    /// One piece of a geometry string: a kind (`b`, `u`, `f`, `x`) and what its brackets hold.
    fn pieces(read: &str) -> Vec<(char, String)> {
        let mut out = Vec::new();
        let mut rest = read;
        while let Some(kind) = rest.chars().next() {
            let open = rest.find('[').expect("a piece opens its brackets");
            let close = rest.find(']').expect("a piece closes its brackets");
            assert_eq!(open, 1, "one letter names a piece: {rest}");
            out.push((kind, rest[open + 1..close].to_string()));
            rest = &rest[close + 1..];
        }
        out
    }

    /// The geometry of the barcode read and of the other read, as simpleaf reads them.
    fn reads(layout: &Layout) -> (String, String) {
        let g = layout.geometry;
        let (one, two) = g.split_once("}2{").expect("a geometry describes read 1 and read 2");
        let one = one.strip_prefix("1{").expect("read 1 first").to_string();
        let two = two.strip_suffix('}').expect("read 2 closes the geometry").to_string();
        match layout.read {
            Mate::R1 => (one, two),
            Mate::R2 => (two, one),
        }
    }

    /// Whether the geometry's pieces describe the layout's segments, base for base. Accepted
    /// encodings: start offsets as a leading `x[0-n]` or as a range on the first barcode
    /// (`b[8-11]`); adjacent barcodes merged into one `b`; a linker as `f[its sequence]` when it
    /// has one alternative, or as `x[its length]`. The tail is open-ended, so it has no piece.
    fn agrees(layout: &Layout) -> Result<(), String> {
        let (barcode_read, other_read) = reads(layout);
        if other_read != "r:" {
            return Err(format!("the cDNA read is {other_read:?}, not r:"));
        }
        let body = barcode_read.strip_suffix("x:").ok_or("the barcode read ends open (x:)")?;
        let mut got = pieces(body);
        let offset_max = layout.offsets.iter().copied().max().unwrap_or(0);
        let mut first_range = None;
        if offset_max > 0 {
            if got.first().is_some_and(|(k, v)| *k == 'x' && *v == format!("0-{offset_max}")) {
                got.remove(0);
            } else {
                first_range = Some(offset_max);
            }
        }
        let mut want = layout.segments.iter().filter(|s| !matches!(s, Segment::Tail(_))).peekable();
        for (index, (kind, value)) in got.iter().enumerate() {
            match *kind {
                'b' => {
                    let (low, high) = value.split_once('-').unwrap_or((value, value));
                    let (low, high): (usize, usize) = (low.parse().unwrap(), high.parse().unwrap());
                    let mut len = 0;
                    while len < low {
                        match want.next() {
                            Some(
                                Segment::Barcode { len: l, .. } | Segment::Unlisted { len: l, .. },
                            ) => {
                                len += l;
                            }
                            other => return Err(format!("b[{value}] meets {other:?}")),
                        }
                    }
                    let range = if index == 0 { first_range.take().unwrap_or(0) } else { 0 };
                    if len != low || high != low + range {
                        return Err(format!("b[{value}]: the segments give {len} (+{range})"));
                    }
                }
                'u' => match want.next() {
                    Some(Segment::Umi(l)) if value.parse() == Ok(*l) => {}
                    other => return Err(format!("u[{value}] meets {other:?}")),
                },
                'f' => match want.next() {
                    Some(Segment::Fixed([only])) if only == &value.as_bytes() => {}
                    other => return Err(format!("f[{value}] meets {other:?}")),
                },
                'x' => match want.next() {
                    Some(s @ Segment::Fixed(_)) if value.parse() == Ok(s.len()) => {}
                    other => return Err(format!("x[{value}] meets {other:?}")),
                },
                _ => return Err(format!("unknown piece {kind}")),
            }
        }
        if let Some(left) = want.next() {
            return Err(format!("{left:?} has no piece in the geometry"));
        }
        if first_range.is_some() {
            return Err(format!("offsets up to {offset_max} are not in the geometry"));
        }
        Ok(())
    }

    /// simpleaf receives the geometry string as it is, so a geometry that disagrees with the
    /// segments tchem detects would print a wrong command with no error.
    #[test]
    fn every_geometry_describes_its_segments() {
        let wrong: Vec<String> = LAYOUTS
            .iter()
            .filter_map(|l| agrees(l).err().map(|e| format!("{}: {} ({e})", l.name, l.geometry)))
            .collect();
        assert!(LAYOUTS.len() >= 25, "every layout is tested");
        assert!(
            wrong.is_empty(),
            "geometries that disagree with their segments:\n{}",
            wrong.join("\n")
        );
    }
}
