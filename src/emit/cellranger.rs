// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! `cellranger count`. Cell Ranger finds a sample's reads by file name in a directory, not by the
//! order given, so the files must be Illumina-named and the barcode read must be the one named
//! R1: a swapped pair would be counted with the cDNA as barcodes, so it is refused.

use super::{
    Arg, EmitRefusal, Emitted, Matched, Pipeline, Placeholder, cdna, emitted, path_text, unsplit,
    unsupported,
};
use std::path::Path;

const REFERENCE: Placeholder = Placeholder {
    name: "CELLRANGER_REF",
    meaning: "a 10x-compatible transcriptome reference directory",
};

/// Digits of an Illumina file name's lane (`_L001`).
const LANE_DIGITS: usize = 3;

/// The directory and sample of an Illumina-named FASTQ for `mate` (`R1` or `R2`):
/// `<sample>_S<n>_L<lane>_<mate>_001.fastq[.gz]`, which is what `--fastqs`/`--sample` select.
fn illumina(file: &Path, mate: &str) -> Result<(String, String), EmitRefusal> {
    assert!(mate == "R1" || mate == "R2", "a read mate, not {mate}");
    let refuse = || EmitRefusal::NotIlluminaNamed(file.display().to_string());
    let name = file.file_name().and_then(|n| n.to_str()).ok_or_else(refuse)?;
    let stem = [".fastq.gz", ".fastq"]
        .iter()
        .find_map(|ext| name.strip_suffix(&format!("_{mate}_001{ext}")))
        .ok_or_else(refuse)?;
    // `_L` and three digits, then `_S` and a number, from the end.
    let (rest, lane) = stem.rsplit_once("_L").ok_or_else(refuse)?;
    let (sample, number) = rest.rsplit_once("_S").ok_or_else(refuse)?;
    let digits = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    if lane.len() != LANE_DIGITS || !digits(lane) || !digits(number) || sample.is_empty() {
        return Err(refuse());
    }
    let parent = file.parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or(Path::new("."));
    let directory = unsplit(path_text(parent)?, "--fastqs")?;
    let sample = unsplit(sample.to_owned(), "--sample")?;
    assert!(!sample.is_empty() && name.starts_with(&sample), "the sample prefixes the name");
    Ok((directory, sample))
}

/// Cell Ranger's run id: the sample with anything outside `[A-Za-z0-9_-]` made `_`.
fn run_id(sample: &str) -> String {
    let id: String = sample
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() || c == '_' || c == '-' { c } else { '_' })
        .collect();
    assert_eq!(id.chars().count(), sample.chars().count(), "one character for each");
    assert!(id.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'-'), "{id}");
    id
}

pub(super) fn emit(m: &Matched<'_>) -> Result<Emitted, EmitRefusal> {
    let Some(chemistry) = m.structure.cellranger_chemistry else {
        return Err(unsupported(
            Pipeline::CellRanger,
            m,
            "Cell Ranger has no count chemistry for it",
        ));
    };
    if m.swapped {
        return Err(EmitRefusal::SwappedForCellRanger);
    }
    let (directory, sample) = illumina(m.barcode_file, "R1")?;
    let (cdna_directory, cdna_sample) = illumina(cdna(m)?, "R2")?;
    if (cdna_directory.as_str(), cdna_sample.as_str()) != (directory.as_str(), sample.as_str()) {
        return Err(EmitRefusal::InputFormat(
            "R1 and R2 are not one sample in one directory, as Cell Ranger needs".into(),
        ));
    }
    let text = |s: &str| Arg::Text(s.to_owned());
    let argv = vec![
        text("cellranger"),
        text("count"),
        text("--id"),
        Arg::Text(run_id(&sample)),
        text("--transcriptome"),
        Arg::Var(REFERENCE),
        text("--fastqs"),
        Arg::Text(directory.clone()),
        text("--sample"),
        Arg::Text(sample.clone()),
        text("--chemistry"),
        text(chemistry),
        text("--create-bam"),
        text("true"),
    ];
    let mut e = emitted(Pipeline::CellRanger, m, argv);
    e.notes.push(format!("Cell Ranger reads every lane of sample {sample} in {directory}"));
    assert!(e.placeholders == [REFERENCE], "Cell Ranger needs the reference alone");
    assert!(!m.swapped, "a swapped pair never reaches here");
    Ok(e)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn illumina_names_give_the_directory_and_sample() {
        let ok = illumina(Path::new("/runs/a b/pbmc_1k_v3_S1_L001_R1_001.fastq.gz"), "R1").unwrap();
        assert_eq!(ok, ("/runs/a b".to_owned(), "pbmc_1k_v3".to_owned()));
        assert_eq!(illumina(Path::new("x_S12_L004_R2_001.fastq"), "R2").unwrap().0, ".");
        for bad in [
            "SRR1_1.fastq.gz",
            "x_S1_L01_R1_001.fastq.gz",
            "x_L001_R1_001.fastq.gz",
            "x_S1_L001_R2_001.fastq.gz",
            "_S1_L001_R1_001.fastq.gz",
        ] {
            assert!(illumina(Path::new(bad), "R1").is_err(), "{bad}");
        }
        assert_eq!(run_id("pbmc 1k.v3"), "pbmc_1k_v3");
    }
}
