// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Why a command line was refused, in words the user can act on.

use super::{ARGUMENTS_MAX, FILES_MAX};
use std::fmt;

/// Why the command line was refused.
#[derive(Debug, PartialEq, Eq)]
pub enum ArgumentRefusal {
    Unknown(String),
    Repeated(&'static str),
    MissingValue(&'static str),
    BadValue {
        flag: &'static str,
        value: String,
        expected: &'static str,
    },
    Missing(&'static str),
    PackWithReads,
    TooMany,
    /// Files given both as `--r1`/`--r2` and by position.
    MixedInputs,
    /// More files than one lane's R1, R2, I1 and I2.
    TooManyFiles,
    /// `--emit` with a run by accession: pipelines need the files on disk.
    EmitWithAccession,
}

impl fmt::Display for ArgumentRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unknown(arg) => write!(f, "unexpected argument '{arg}'"),
            Self::Repeated(flag) => write!(f, "{flag} given more than once"),
            Self::MissingValue(flag) => write!(f, "{flag} needs a value"),
            Self::BadValue { flag, value, expected } => {
                write!(f, "{flag} '{value}': expected {expected}")
            }
            Self::Missing(flag) => write!(f, "{flag} is required"),
            Self::PackWithReads => write!(f, "--pack takes no --r1 or --r2"),
            Self::TooMany => write!(f, "more than {ARGUMENTS_MAX} arguments"),
            Self::MixedInputs => {
                write!(f, "give files either as --r1/--r2 or by position, not both")
            }
            Self::EmitWithAccession => write!(
                f,
                "--emit needs the files on disk: download the run, then emit for its files"
            ),
            Self::TooManyFiles => {
                write!(f, "more than {FILES_MAX} files: give one lane's R1, R2, I1 and I2")
            }
        }
    }
}

impl std::error::Error for ArgumentRefusal {}
