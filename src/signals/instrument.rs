// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Which instrument read the sample, from what the source records: a BAM read group's platform
//! tag when there is one, else the shape of the first read name. Archives often keep the original
//! name after their own accession (`@SRR123.1 A00836:478:...`), and often strip it (`@SRR123.1`,
//! `1/1`): then the answer is unknown, never a guess. Detection never reads any of this.
//!
//! Each shape was checked on real names (8 Oct 2026, ENA runs and heads cached for the benches);
//! archive metadata disagreed with the names often enough (an "Element AVITI" run with `NextSeq`
//! names, "UG 100" runs with `NovaSeq` names) that the names are reported as what they are.

use crate::input::fastq::Origin;

/// The sequencing platform.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum Platform {
    /// Illumina.
    Illumina,
    /// Element Biosciences (AVITI).
    Element,
    /// MGI / Complete Genomics DNBSEQ (BGISEQ, G400, T7, T20).
    Mgi,
    /// Ultima Genomics.
    Ultima,
    /// Pacific Biosciences.
    PacBio,
    /// Oxford Nanopore Technologies.
    OxfordNanopore,
}

impl Platform {
    /// The platform as the report names it.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Platform::Illumina => "Illumina",
            Platform::Element => "Element",
            Platform::Mgi => "MGI DNBSEQ",
            Platform::Ultima => "Ultima",
            Platform::PacBio => "PacBio",
            Platform::OxfordNanopore => "Oxford Nanopore",
        }
    }
}

/// Where the answer came from.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum Evidence {
    /// The first read's name.
    ReadName,
    /// A BAM `@RG` line's `PL`/`PM` tags.
    ReadGroup,
}

/// The instrument, as far as the source says.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
#[non_exhaustive]
pub struct Instrument {
    /// The platform, or None when nothing the source records says.
    pub platform: Option<Platform>,
    /// The model, when the name or tags fix it (e.g. `NovaSeq 6000`).
    pub model: Option<String>,
    /// The instrument, run or flowcell identifier the name carries.
    pub id: Option<String>,
    /// What the answer rests on; None when unknown.
    pub evidence: Option<Evidence>,
}

/// Illumina instrument-ID prefixes and their models; no two match one ID. An ID is
/// a prefix then digits only; a custom name (`HISEQ2500-10`) gives Illumina with no model.
const ILLUMINA_MODELS: [(&str, &str); 13] = [
    ("LH", "NovaSeq X"),
    ("VH", "NextSeq 2000"),
    ("VL", "NextSeq 1000"),
    ("NB", "NextSeq 500/550"),
    ("NS", "NextSeq 500/550"),
    ("MN", "MiniSeq"),
    ("FS", "iSeq 100"),
    ("SN", "HiSeq 2000"),
    ("A", "NovaSeq 6000"),
    ("M", "MiSeq"),
    ("K", "HiSeq 3000/4000"),
    ("E", "HiSeq X"),
    ("D", "HiSeq 2500"),
];

/// Fields of an Illumina read name: instrument, run, flowcell, lane, tile, x, y.
const ILLUMINA_NAME_FIELDS: usize = 7;
/// Digits of an Ultima run id (`415630_1-…`).
const ULTIMA_RUN_DIGITS: usize = 6;
/// Digits of a PacBio movie's date (YYMMDD).
const PACBIO_DATE_DIGITS: usize = 6;
/// Lengths of a UUID's hex groups (8-4-4-4-12), as Oxford Nanopore read ids are.
const UUID_GROUPS: [usize; 5] = [8, 4, 4, 4, 12];

/// The instrument the origin records: the read group's tags first, else the read name.
#[must_use]
pub fn of(origin: &Origin) -> Instrument {
    let instrument = match &origin.platform_tag {
        Some(tag) => from_read_group(tag, origin.model_tag.as_deref()),
        None => None,
    }
    .or_else(|| origin.first_name.as_deref().and_then(from_name))
    .unwrap_or_default();
    assert!(instrument.platform.is_some() == instrument.evidence.is_some(), "evidence iff known");
    assert!(instrument.platform.is_some() || instrument.model.is_none(), "no model without one");
    instrument
}

/// A BAM read group's `PL` tag (SAM spec values) and `PM` model.
fn from_read_group(platform: &str, model: Option<&str>) -> Option<Instrument> {
    let platform = match platform.to_ascii_uppercase().as_str() {
        "ILLUMINA" => Platform::Illumina,
        "ELEMENT" => Platform::Element,
        "DNBSEQ" | "BGI" | "MGI" => Platform::Mgi,
        "ULTIMA" => Platform::Ultima,
        "PACBIO" => Platform::PacBio,
        "ONT" => Platform::OxfordNanopore,
        _ => return None,
    };
    let model = model.map(str::to_owned);
    let found = Instrument {
        platform: Some(platform),
        model,
        id: None,
        evidence: Some(Evidence::ReadGroup),
    };
    Some(found)
}

/// The original name inside an archive's: `SRR123.1 A00836:...` gives `A00836:...`; a bare
/// accession (`SRR123.1`) gives None, its original name stripped.
fn original(name: &str) -> Option<&str> {
    let mut tokens = name.split_whitespace();
    let first = tokens.next()?;
    let archived = ["SRR", "ERR", "DRR"].iter().any(|p| first.starts_with(p))
        && first.split_once('.').is_some_and(|(_, n)| n.bytes().all(|b| b.is_ascii_digit()));
    if archived { tokens.next() } else { Some(first) }
}

/// The platform a read name's shape shows, if any.
fn from_name(name: &str) -> Option<Instrument> {
    let token = original(name)?;
    let shape = illumina_like(token)
        .or_else(|| mgi(token))
        .or_else(|| ultima(token))
        .or_else(|| pacbio(token))
        .or_else(|| nanopore(name));
    assert!(shape.as_ref().is_none_or(|i| i.evidence == Some(Evidence::ReadName)), "from the name");
    assert!(
        shape.as_ref().is_none_or(|i| i.id.as_deref().is_some_and(|id| name.contains(id))),
        "the instrument id is part of the name"
    );
    shape
}

fn named(platform: Platform, model: Option<&str>, id: &str) -> Instrument {
    assert!(!id.is_empty(), "a name carries an identifier");
    let found = Instrument {
        platform: Some(platform),
        model: model.map(str::to_owned),
        id: Some(id.to_owned()),
        evidence: Some(Evidence::ReadName),
    };
    assert!(
        found.id.as_ref().is_some_and(|i| i.len() <= crate::input::fastq::KEPT_NAME_BYTES_MAX),
        "a kept name is bounded"
    );
    found
}

/// `instrument:run:flowcell:lane:tile:x:y` (and a `/1`), Illumina's and Element's shape: lane,
/// tile and coordinates are numbers. Illumina's run is a number; Element's instrument is `AV…`.
/// bcl2fastq and BCL Convert may add the UMI as an eighth field (`…:y:ACGTACGT`). A DNBSEQ
/// flowcell written in this shape (`V350033550:…`) is not Illumina: unknown, not a guess.
fn illumina_like(token: &str) -> Option<Instrument> {
    let token = token.split('/').next().unwrap_or(token);
    let fields: Vec<&str> = token.split(':').collect();
    let number = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    let bases = |s: &str| !s.is_empty() && s.bytes().all(|b| b"ACGTN+".contains(&b));
    let shape = match fields.as_slice() {
        [first @ .., umi] if first.len() == ILLUMINA_NAME_FIELDS && bases(umi) => first,
        all => all,
    };
    let [instrument, run, flowcell, lane, tile, x, y] = shape else { return None };
    if instrument.is_empty() || flowcell.is_empty() || ![lane, tile, x, y].iter().all(|f| number(f))
    {
        return None;
    }
    if dnbseq_flowcell(instrument) {
        return None;
    }
    if instrument.starts_with("AV") {
        return Some(named(Platform::Element, Some("AVITI"), instrument));
    }
    if !number(run) {
        return None;
    }
    let upper = instrument.to_ascii_uppercase();
    let model = ILLUMINA_MODELS.iter().find_map(|(prefix, model)| {
        let digits = upper.strip_prefix(prefix)?;
        (!digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit())).then_some(*model)
    });
    // The table's prefixes never overlap: an id names at most one model.
    let matching = ILLUMINA_MODELS
        .iter()
        .filter(|(prefix, _)| {
            upper
                .strip_prefix(prefix)
                .is_some_and(|d| !d.is_empty() && d.bytes().all(|b| b.is_ascii_digit()))
        })
        .count();
    assert!(matching <= 1, "an instrument id matches one model prefix at most: {instrument}");
    assert_eq!(model.is_some(), matching == 1, "the model is the one prefix that matched");
    Some(named(Platform::Illumina, model, instrument))
}

/// A DNBSEQ flowcell id: one or two capitals and at least nine digits (`V350033550`,
/// `E100012345`, `FP200001234`). Illumina instrument ids hold at most eight (`FS10000123`).
fn dnbseq_flowcell(id: &str) -> bool {
    let letters = id.bytes().take_while(u8::is_ascii_uppercase).count();
    let digits = id.get(letters..).unwrap_or("");
    (1..=DNBSEQ_FLOWCELL_LETTERS_MAX).contains(&letters)
        && digits.len() >= DNBSEQ_FLOWCELL_DIGITS_MIN
        && digits.bytes().all(|b| b.is_ascii_digit())
}

/// Capitals that open a DNBSEQ flowcell id at most (`V…`, `FP…`).
const DNBSEQ_FLOWCELL_LETTERS_MAX: usize = 2;
/// Digits in a DNBSEQ flowcell id at least; Illumina's longest instrument number is 8 (iSeq).
const DNBSEQ_FLOWCELL_DIGITS_MIN: usize = 9;

/// DNBSEQ: flowcell, then `L<lane>C<column>R<row>` (`V350033550L1C001R001…`, `E250067724L1…`).
/// The flowcell may hold an `L` itself (`DP8400021809BLL1C001R005…`): the lane marker is the `L`
/// that the rest of the shape follows.
fn mgi(token: &str) -> Option<Instrument> {
    let digits = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    let lane_at = token.match_indices('L').map(|(at, _)| at).find(|&at| {
        let rest = token.get(at + 1..).unwrap_or("");
        let Some(c) = rest.find('C') else { return false };
        let column = rest.get(c + 1..).unwrap_or("");
        (1..=2).contains(&c)
            && rest.get(..c).is_some_and(digits)
            && column.get(..3).is_some_and(digits)
            && column.get(3..4) == Some("R")
            && column.get(4..7).is_some_and(digits)
    })?;
    let flowcell = token.get(..lane_at)?;
    let shaped = flowcell.bytes().next().is_some_and(|b| b.is_ascii_uppercase())
        && flowcell.bytes().any(|b| b.is_ascii_digit());
    shaped.then(|| named(Platform::Mgi, None, flowcell))
}

/// Ultima: a six-digit run, `_` and a number, then `-` (`415630_1-…`, `017464_1-UGAv3-…`).
fn ultima(token: &str) -> Option<Instrument> {
    let (run, rest) = token.split_once('_')?;
    let (part, _) = rest.split_once('-')?;
    let digits = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    let shaped = run.len() == ULTIMA_RUN_DIGITS && digits(run) && digits(part);
    shaped.then(|| named(Platform::Ultima, None, run))
}

/// PacBio: `m<model><serial>[e]_<date>_<time>[_s<cell>]/<hole>/…`; m84 is Revio, m64 Sequel
/// II/IIe, m54 Sequel. Revio adds the SMRT Cell slot (`m84047_230404_192312_s1/…`).
fn pacbio(token: &str) -> Option<Instrument> {
    let movie = token.split('/').next()?;
    let parts: Vec<&str> = movie.split('_').collect();
    let digits = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    let (instrument, date, time) = match parts.as_slice() {
        [instrument, date, time] => (instrument, date, time),
        [instrument, date, time, cell] if cell.strip_prefix('s').is_some_and(digits) => {
            (instrument, date, time)
        }
        _ => return None,
    };
    if !instrument.starts_with('m')
        || date.len() != PACBIO_DATE_DIGITS
        || !digits(date)
        || !digits(time)
    {
        return None;
    }
    let model = match instrument.get(1..3) {
        Some("84") => Some("Revio"),
        Some("64") => Some("Sequel II/IIe"),
        Some("54") => Some("Sequel"),
        _ => None,
    };
    Some(named(Platform::PacBio, model, instrument))
}

/// Oxford Nanopore: a UUID read id (8-4-4-4-12 hex) anywhere in the name.
fn nanopore(name: &str) -> Option<Instrument> {
    let uuid = name.split(|c: char| !(c.is_ascii_hexdigit() || c == '-')).find(|t| {
        let groups: Vec<usize> = t.split('-').map(str::len).collect();
        groups == UUID_GROUPS
    })?;
    Some(named(Platform::OxfordNanopore, None, uuid))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn by_name(name: &str) -> (Option<Platform>, Option<String>) {
        let i = of(&Origin { first_name: Some(name.to_owned()), ..Origin::default() });
        (i.platform, i.model)
    }

    /// Every shape on a real name (ENA runs and bench heads, 8 Oct 2026).
    #[test]
    fn real_names() {
        use Platform::*;
        let some = |s: &str| Some(s.to_owned());
        let cases: [(&str, Option<Platform>, Option<String>); 19] = [
            (
                "A00836:478:HGTYWDSXY:1:1101:10004:1000 1:N:0:ACGT",
                Some(Illumina),
                some("NovaSeq 6000"),
            ),
            (
                "SRR34855738.1 LH00746:77:22NHWHLT4:1:1101:1000:1000",
                Some(Illumina),
                some("NovaSeq X"),
            ),
            (
                "ERR10026091.1 VH00574:5:AAAL2LYHV:1:1101:37621:1000",
                Some(Illumina),
                some("NextSeq 2000"),
            ),
            (
                "ERR13764620.1 NB551556:794:HNMT7BGXV:2:11101:12562:1137",
                Some(Illumina),
                some("NextSeq 500/550"),
            ),
            (
                "ERR10038570.1 M01332:100:000000000-JH7RJ:1:1101:10009:10813/1",
                Some(Illumina),
                some("MiSeq"),
            ),
            ("ERR10015002.1 HISEQ2500-10:576:H3VGGBCX2:1:1101:1349:2166", Some(Illumina), None),
            (
                "ERR14064734.1 AV240401:AVT0016:2345645132:1:10102:0001:0213/1",
                Some(Element),
                some("AVITI"),
            ),
            ("SRR25532325.1 DP8400021809BLL1C001R00500000030/1", Some(Mgi), None),
            ("ERR10030051.1 V350033550L1C001R00100000882:0:0:0:0/1", Some(Mgi), None),
            ("ERR10038832.1 CL100076516L2C001R001_100024/1", Some(Mgi), None),
            ("SRR36319711.1 400675_1-Z0101-1661031103/1", Some(Ultima), None),
            ("SRR28496332.1 017464_1-UGAv3-207-0000000100/1", Some(Ultima), None),
            (
                "ERR10824956.1 m64270e_220523_002302/69/27170_29083",
                Some(PacBio),
                some("Sequel II/IIe"),
            ),
            ("m84012_230101_101010/1/ccs", Some(PacBio), some("Revio")),
            // A real Revio name carries the cell slot.
            ("m84047_230404_192312_s1/12345/ccs", Some(PacBio), some("Revio")),
            ("A00123:45:HLKXXDSXY:1:1101:1000:1000:ACGTACGT", Some(Illumina), some("NovaSeq 6000")),
            ("a00123:45:hlkxxdsxy:1:1101:1000:1000", Some(Illumina), some("NovaSeq 6000")),
            (
                "ERR16927611.1 67fd17f7-2c4a-4b6e-9b1a-0123456789ab runid=x",
                Some(OxfordNanopore),
                None,
            ),
            (
                "ERR10092105.1 62:393|235552a6-6cec-48ba-8315-3fb59b2d8046",
                Some(OxfordNanopore),
                None,
            ),
        ];
        for (name, platform, model) in cases {
            assert_eq!(by_name(name), (platform, model), "{name}");
        }
    }

    /// Names that say nothing stay unknown: stripped by an archive, processed, or numbered.
    #[test]
    fn stripped_names_are_unknown() {
        for name in [
            "SRR36374443.1 1/1",
            "SRR12345.1",
            "1/1",
            "OTU426/1",
            "KNX_DTOL14815334_transcript/691061",
            "V143:418786:418786:418786:2:2:1:4529:9:1:1780:N:0.71:CATCTCTGGCTG",
            // A DNBSEQ flowcell in Illumina's shape is not Illumina.
            "V350033550:1:V350033550:1:1101:1000:1000",
            "E100012345:1:E100012345:1:1101:1000:1000",
            "",
        ] {
            assert_eq!(by_name(name), (None, None), "{name}");
        }
        let tagged = Origin {
            first_name: Some("SRR1.1".into()),
            platform_tag: Some("ONT".into()),
            model_tag: Some("PromethION".into()),
        };
        assert_eq!(of(&tagged).platform, Some(Platform::OxfordNanopore), "the read group wins");
        assert_eq!(of(&tagged).model.as_deref(), Some("PromethION"));
    }
}
