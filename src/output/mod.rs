// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! The two ways a result leaves the tool: a report for people (`text`) and one JSON object for
//! pipelines (`json`). Both render the same `detect::Report`.

pub mod json;
pub mod text;
