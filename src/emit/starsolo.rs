// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! `STARsolo`, `--soloType CB_UMI_Simple`: one barcode then a UMI from the barcode read's first
//! base. The cDNA file comes first in `--readFilesIn` and the barcode file last, as `STARsolo`
//! requires; positions are 1-based.

use super::{
    Arg, EmitRefusal, Emitted, Matched, Pipeline, Placeholder, Stored, cdna, emitted, path_text,
    single_barcode, strand, unsplit, unsupported,
};
use crate::Strand;

const INDEX: Placeholder = Placeholder {
    name: "STAR_INDEX",
    meaning: "a STAR genome index built with a gene annotation (GTF)",
};

/// The `--readFilesCommand` for the two files, if they need one. One command serves both files,
/// so they must be stored alike. `gunzip -c`, not `zcat`: macOS's zcat reads only `.Z` files.
fn read_command(stored: (Stored, Option<Stored>)) -> Result<Option<&'static str>, EmitRefusal> {
    let command = match stored {
        (Stored::Gzip, Some(Stored::Gzip)) => Some("gunzip -c"),
        (Stored::Plain, Some(Stored::Plain)) => None,
        _ => {
            let why = "one file is gzipped and the other is not; STARsolo reads both one way";
            return Err(EmitRefusal::InputFormat(why.into()));
        }
    };
    Ok(command)
}

pub(super) fn emit(m: &Matched<'_>) -> Result<Emitted, EmitRefusal> {
    let Some((barcode, listed, umi)) = single_barcode(m.structure) else {
        return Err(unsupported(
            Pipeline::Starsolo,
            m,
            "CB_UMI_Simple takes one barcode and a UMI from the read's first base, and tchem \
             does not write CB_UMI_Complex yet",
        ));
    };
    let cdna = cdna(m)?;
    let strand = match strand(m)? {
        Strand::Forward => "Forward",
        Strand::Reverse => "Reverse",
    };
    let list = match (listed, m.layout.barcode_lists.as_slice()) {
        (true, [list]) => path_text(list)?,
        (false, []) => "None".to_owned(),
        _ => unreachable!("{}: one listed barcode has one list", m.layout.name),
    };
    let nuclei = crate::layouts::nuclei(m.structure);
    let text = |s: &str| Arg::Text(s.to_owned());
    let mut argv = vec![
        text("STAR"),
        text("--genomeDir"),
        Arg::Var(INDEX),
        text("--readFilesIn"),
        Arg::Text(unsplit(path_text(cdna)?, "--readFilesIn")?),
        Arg::Text(unsplit(path_text(m.barcode_file)?, "--readFilesIn")?),
    ];
    if let Some(command) = read_command(m.stored)? {
        argv.extend([text("--readFilesCommand"), text(command)]);
    }
    argv.extend([
        text("--soloType"),
        text("CB_UMI_Simple"),
        text("--soloCBwhitelist"),
        Arg::Text(list),
        text("--soloCBstart"),
        text("1"),
        text("--soloCBlen"),
        Arg::Text(barcode.to_string()),
        text("--soloUMIstart"),
        Arg::Text((barcode + 1).to_string()),
        text("--soloUMIlen"),
        Arg::Text(umi.to_string()),
        // Reads longer than barcode + UMI (often R1 is 50+ bases) are fine: do not check length.
        text("--soloBarcodeReadLength"),
        text("0"),
        text("--soloStrand"),
        text(strand),
        text("--soloFeatures"),
        // Nuclei: GeneFull counts intronic reads too, as Cell Ranger ARC does; Gene would drop
        // the unspliced majority.
        text(if nuclei { "GeneFull" } else { "Gene" }),
    ]);
    let mut e = emitted(Pipeline::Starsolo, m, argv);
    if nuclei {
        e.notes.push(
            "nuclei (Multiome): GeneFull counts exons and introns; read Solo.out/GeneFull".into(),
        );
    }
    if !listed {
        e.notes.push("no barcode list exists for this kit: STARsolo counts every barcode".into());
    }
    assert!(e.placeholders == [INDEX], "STARsolo needs the index alone");
    assert!(barcode >= 1 && umi >= 1, "positions are positive");
    Ok(e)
}
