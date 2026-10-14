// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! `kb count` with a kallisto technology string derived from the layout, `barcode:umi:cDNA`,
//! each `file,start,stop` (0-based, stop exclusive, 0 for "to the end"). The barcode file is
//! given first (file 0) and the cDNA second (file 1).

use super::{
    Arg, EmitRefusal, Emitted, Matched, Pipeline, Placeholder, cdna, emitted, path_text,
    single_barcode, strand, unsupported,
};
use crate::Strand;

const INDEX: Placeholder =
    Placeholder { name: "KB_INDEX", meaning: "a kallisto index (kb ref -i)" };
const T2G: Placeholder =
    Placeholder { name: "KB_T2G", meaning: "the transcript-to-gene map (kb ref -g)" };
/// Nuclei count with kb's nac workflow, from an index built by `kb ref --workflow nac`.
const NAC_INDEX: Placeholder =
    Placeholder { name: "KB_INDEX", meaning: "a kallisto nac index (kb ref --workflow nac -i)" };
const MATURE: Placeholder = Placeholder {
    name: "KB_C1",
    meaning: "the mature transcripts list (kb ref --workflow nac -c1)",
};
const NASCENT: Placeholder = Placeholder {
    name: "KB_C2",
    meaning: "the nascent transcripts list (kb ref --workflow nac -c2)",
};

/// The technology string for one barcode of `barcode` bases then a UMI of `umi`.
fn technology(barcode: usize, umi: usize) -> String {
    assert!(barcode >= 1 && umi >= 1, "a barcode and a UMI occupy bases");
    let tech = format!("0,0,{barcode}:0,{barcode},{}:1,0,0", barcode + umi);
    assert_eq!(tech.split(':').count(), 3, "barcode, UMI and cDNA");
    tech
}

pub(super) fn emit(m: &Matched<'_>) -> Result<Emitted, EmitRefusal> {
    let Some((barcode, listed, umi)) = single_barcode(m.structure) else {
        return Err(unsupported(
            Pipeline::Kb,
            m,
            "tchem writes kb technology strings for one barcode and a UMI from the read's first \
             base so far",
        ));
    };
    let cdna = cdna(m)?;
    let strand = match strand(m)? {
        Strand::Forward => "forward",
        Strand::Reverse => "reverse",
    };
    let nuclei = crate::layouts::nuclei(m.structure);
    let text = |s: &str| Arg::Text(s.to_owned());
    let mut argv = vec![text("kb"), text("count")];
    if nuclei {
        // Nuclei: count nascent (intronic) and mature reads together, as Cell Ranger ARC does.
        argv.extend([text("--workflow"), text("nac"), text("-i"), Arg::Var(NAC_INDEX)]);
        argv.extend([text("-g"), Arg::Var(T2G), text("-c1"), Arg::Var(MATURE)]);
        argv.extend([text("-c2"), Arg::Var(NASCENT), text("--sum=total")]);
    } else {
        argv.extend([text("-i"), Arg::Var(INDEX), text("-g"), Arg::Var(T2G)]);
    }
    argv.extend([text("-x"), Arg::Text(technology(barcode, umi))]);
    match (listed, m.layout.barcode_lists.as_slice()) {
        (true, [list]) => argv.extend([text("-w"), Arg::Text(path_text(list)?)]),
        (false, []) => {}
        _ => unreachable!("{}: one listed barcode has one list", m.layout.name),
    }
    argv.extend([
        text("--strand"),
        text(strand),
        text("-o"),
        text("kb_count"),
        Arg::Text(path_text(m.barcode_file)?),
        Arg::Text(path_text(cdna)?),
    ]);
    let mut e = emitted(Pipeline::Kb, m, argv);
    if !listed {
        e.notes.push("no barcode list exists for this kit: kb builds one from the reads".into());
    }
    if nuclei {
        e.notes.push("nuclei (Multiome): the nac workflow sums mature and nascent counts".into());
        assert!(e.placeholders == [NAC_INDEX, T2G, MATURE, NASCENT], "nac needs its four files");
    } else {
        assert!(e.placeholders == [INDEX, T2G], "kb needs the index and the map");
    }
    assert!(m.barcode_file != cdna, "two different files");
    Ok(e)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn technology_strings_match_kb_builtins() {
        // kb --list (0.30.2): 10XV2 is 0,0,16 / 0,16,26 / 1; 10XV3 is 0,0,16 / 0,16,28 / 1;
        // DROPSEQ is 0,0,12 / 0,12,20.
        assert_eq!(technology(16, 10), "0,0,16:0,16,26:1,0,0");
        assert_eq!(technology(16, 12), "0,0,16:0,16,28:1,0,0");
        assert_eq!(technology(12, 8), "0,0,12:0,12,20:1,0,0");
    }
}
