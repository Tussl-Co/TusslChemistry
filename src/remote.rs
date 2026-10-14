// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! A public run by accession (`tchem SRR12345678`): ENA's file report names the run's FASTQ files,
//! and each is streamed through `curl` until the sample is full, then the download is dropped.
//! A 100,000-read sample costs tens of megabytes of a run that may be 100 GB.
//!
//! Everything here crosses the network, so all of it is untrusted: the accession is checked
//! against its pattern before it goes into a URL, ENA's answer is capped and parsed strictly,
//! every file URL must be on ENA's FASTQ host with a plain path, and each stream is capped in
//! bytes. `curl` runs with its arguments passed directly (no shell), HTTPS only, with timeouts.
//! It is a separate program so the crate keeps its one dependency; a missing `curl` is a refusal.

use crate::cli::FILES_MAX; // a run has at most a lane's R1, R2, I1 and I2
use std::fmt;
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use tussl_chemistry::Reads;

/// Run and experiment accessions of the three INSDC archives (SRA, ENA, DDBJ).
const PREFIXES: [&str; 6] = ["SRR", "ERR", "DRR", "SRX", "ERX", "DRX"];
/// Digits after the prefix: SRA's oldest runs have 6, and 10 leaves room for growth.
const DIGITS_MIN: usize = 6;
const DIGITS_MAX: usize = 10;
/// ENA's file report for a run or experiment: one row per run.
const FILE_REPORT: &str = "https://www.ebi.ac.uk/ena/portal/api/filereport?result=read_run\
                           &fields=run_accession,fastq_ftp,instrument_model&format=tsv&accession=";
/// The only host a FASTQ URL may name (ENA serves its FASTQ from here, over HTTPS too).
const FASTQ_HOST: &str = "ftp.sra.ebi.ac.uk/vol1/";
/// ENA's answer is a few hundred bytes per run; an experiment has at most a few hundred runs.
const REPORT_BYTES_MAX: u64 = 1 << 20;
const REPORT_ROWS_MAX: usize = 10_000;
/// Compressed bytes read from one file: 100,000 reads of 150 bases are about 20 MB gzipped.
const STREAM_BYTES_MAX: u64 = 256 << 20;
/// Bytes of an archive's instrument model kept (`Illumina NovaSeq 6000` is 22): longer is not a
/// model name, and it goes into a report line.
const INSTRUMENT_MODEL_BYTES_MAX: usize = 64;
/// Bytes to megabytes (MiB), for messages.
pub const MEGABYTE_SHIFT: u32 = 20;
/// Redirects followed: ENA serves reports and files directly, so a chain is never needed.
const REDIRECTS_MAX: &str = "3";
/// curl's error output kept for a refusal: its first line is the message.
const STDERR_BYTES_MAX: u64 = 4096;
/// Seconds before `curl` gives up on the report and on a stream.
const REPORT_SECONDS_MAX: &str = "120";
const STREAM_SECONDS_MAX: &str = "900";

/// Why a run could not be read by accession.
#[derive(Debug)]
pub enum RemoteRefusal {
    /// `curl` is not installed (or not on PATH).
    CurlMissing,
    /// `curl` failed: the network, a timeout, or an HTTP error.
    Download { what: String, reason: String },
    /// ENA knows no run by this accession.
    NotFound(String),
    /// ENA lists no FASTQ for the run (some runs are SRA archives only).
    NoFastq(String),
    /// ENA's answer was not what its file report promises.
    BadReport(String),
    /// The run has more files than a lane's R1, R2, I1 and I2.
    TooManyFiles { run: String, files: usize },
    /// Reading a stream failed.
    Read { url: String, source: tussl_chemistry::ReadsRefusal },
    /// The sample asked for did not fit in `STREAM_BYTES_MAX` bytes of the file.
    StreamCap { url: String, bytes: u64 },
}

impl fmt::Display for RemoteRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::CurlMissing => write!(f, "reading a run by accession needs curl, not found"),
            Self::Download { what, reason } => write!(f, "could not download {what}: {reason}"),
            Self::NotFound(accession) => write!(f, "ENA knows no public run {accession}"),
            Self::NoFastq(run) => write!(
                f,
                "ENA lists no FASTQ for {run}; fetch the SRA copy with `fasterq-dump \
                 --include-technical {run}` and give the files"
            ),
            Self::BadReport(why) => write!(f, "ENA's file report was not as expected: {why}"),
            Self::TooManyFiles { run, files } => {
                write!(f, "{run} has {files} files; expected a lane's R1, R2, I1 and I2 at most")
            }
            Self::Read { url, source } => write!(f, "{url}: {source}"),
            Self::StreamCap { url, bytes } => write!(
                f,
                "{url}: the reads asked for did not fit in the first {} MB; ask for fewer \
                 (--reads)",
                bytes >> MEGABYTE_SHIFT
            ),
        }
    }
}

impl std::error::Error for RemoteRefusal {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Read { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// Whether `text` is a run or experiment accession: a known prefix, then 6 to 10 digits.
pub fn is_accession(text: &str) -> bool {
    let digits = PREFIXES.iter().find_map(|p| text.strip_prefix(p)).unwrap_or("");
    let accession = (DIGITS_MIN..=DIGITS_MAX).contains(&digits.len())
        && digits.bytes().all(|b| b.is_ascii_digit());
    // Only ASCII letters and digits can follow from a match: safe to put into a URL as is.
    assert!(!accession || text.bytes().all(|b| b.is_ascii_alphanumeric()), "{text}");
    assert!(!accession || text.len() <= 3 + DIGITS_MAX, "{text} is longer than an accession");
    accession
}

/// `curl`, HTTPS only (redirects too, at most `REDIRECTS_MAX`), failing on HTTP errors, quiet but
/// for errors, with a time limit. `--disable` (first, as curl requires) skips the user's
/// `~/.curlrc`, which could otherwise add options such as an output file or insecure TLS.
fn curl(url: &str, seconds: &str) -> Command {
    assert!(url.starts_with("https://"), "only HTTPS URLs reach curl: {url}");
    assert!(!url.bytes().any(|b| b.is_ascii_whitespace() || b.is_ascii_control()), "{url}");
    let mut command = Command::new("curl");
    command
        .args(["--disable", "--silent", "--show-error", "--fail", "--proto", "=https"])
        .args(["--location", "--proto-redir", "=https", "--max-redirs", REDIRECTS_MAX])
        .args(["--max-time", seconds, "--", url])
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    command
}

fn spawn(mut command: Command) -> Result<Child, RemoteRefusal> {
    assert_eq!(command.get_program(), "curl", "only curl is run");
    let child = command.spawn().map_err(|error| match error.kind() {
        io::ErrorKind::NotFound => RemoteRefusal::CurlMissing,
        _ => RemoteRefusal::Download { what: "curl".into(), reason: error.to_string() },
    })?;
    assert!(child.stdout.is_some() && child.stderr.is_some(), "curl's output is piped");
    Ok(child)
}

/// How a curl process ended: stopped by us (no exit code: a signal), done with the whole file, or
/// failed on its own, with its first error line.
enum Ended {
    Stopped,
    Finished,
    Failed(String),
}

/// Stops curl if it is still running and waits for it. Its output pipe must still be open: had it
/// closed first, curl could fail writing before the stop and look like a download error.
fn stop(child: &mut Child) -> Ended {
    assert!(child.stdout.is_none(), "the caller holds curl's output pipe while it stops curl");
    let _ = child.kill(); // fails only when curl has already ended, which `wait` then reports
    let status = match child.wait() {
        Ok(status) => status,
        Err(error) => return Ended::Failed(error.to_string()),
    };
    let ended = match status.code() {
        None => Ended::Stopped,
        Some(0) => Ended::Finished,
        Some(code) => {
            let mut stderr = String::new();
            if let Some(pipe) = child.stderr.take() {
                let _ = pipe.take(STDERR_BYTES_MAX).read_to_string(&mut stderr);
            }
            let first = stderr.lines().next().unwrap_or("").trim();
            Ended::Failed(format!("{first} (curl exit {code})").trim_start().to_owned())
        }
    };
    assert!(!matches!(&ended, Ended::Failed(r) if r.is_empty()), "a failure says why");
    ended
}

/// A reader that passes at most `cap` bytes through, counting them, so a caller can tell a stream
/// that ended from one cut at the cap.
struct Capped<R> {
    inner: R,
    cap: u64,
    passed: u64,
}

impl<R: Read> Capped<R> {
    fn new(inner: R, cap: u64) -> Self {
        assert!(cap > 0, "a cap lets something through");
        assert!(cap <= STREAM_BYTES_MAX, "a cap is one of the named caps, not {cap}");
        Capped { inner, cap, passed: 0 }
    }

    /// Whether the cap was reached: the stream may go on past it.
    fn reached(&self) -> bool {
        assert!(self.passed <= self.cap, "{} bytes passed a {} byte cap", self.passed, self.cap);
        let reached = self.passed == self.cap;
        assert!(!reached || self.passed > 0, "a reached cap let bytes through");
        reached
    }
}

impl<R: Read> Read for Capped<R> {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let room = usize::try_from(self.cap - self.passed).unwrap_or(usize::MAX).min(buffer.len());
        let Some(window) = buffer.get_mut(..room) else { unreachable!("room fits the buffer") };
        let n = if room == 0 { 0 } else { self.inner.read(window)? };
        self.passed += n as u64;
        assert!(n <= room && self.passed <= self.cap, "read {n} past the cap {}", self.cap);
        Ok(n)
    }
}

/// The run's FASTQ URLs (HTTPS, on ENA's host), and the run they belong to. An experiment with
/// several runs gives its first run by accession: one library, one chemistry.
pub fn fastq_urls(accession: &str) -> Result<Run, RemoteRefusal> {
    assert!(is_accession(accession), "{accession} was checked before the URL is built");
    let mut child = spawn(curl(&format!("{FILE_REPORT}{accession}"), REPORT_SECONDS_MAX))?;
    let Some(stdout) = child.stdout.take() else { unreachable!("curl's stdout is piped") };
    // One byte past the cap tells an over-size answer from one that fits exactly.
    let mut stdout = Capped::new(stdout, REPORT_BYTES_MAX + 1);
    let mut report = String::new();
    let read = stdout.read_to_string(&mut report);
    let ended = stop(&mut child);
    if let Ended::Failed(reason) = ended {
        let what = format!("ENA's report for {accession}");
        return Err(RemoteRefusal::Download { what, reason });
    }
    read.map_err(|e| RemoteRefusal::BadReport(e.to_string()))?;
    if stdout.reached() {
        return Err(RemoteRefusal::BadReport(format!("larger than {REPORT_BYTES_MAX} bytes")));
    }
    let run = parse_report(&report, accession)?;
    assert!(run.urls.iter().all(|u| u.starts_with("https://")), "every URL is HTTPS");
    assert!(run.urls.len() <= FILES_MAX, "a run has at most a lane's files");
    Ok(run)
}

/// A run as ENA's file report describes it.
pub struct Run {
    /// The run accession.
    pub run: String,
    /// Its FASTQ files, as checked HTTPS URLs.
    pub urls: Vec<String>,
    /// The instrument model ENA records for it (`instrument_model`), when it records one.
    pub instrument_model: Option<String>,
}

/// ENA's TSV: a header naming `run_accession` and `fastq_ftp`, then one row per run.
fn parse_report(report: &str, accession: &str) -> Result<Run, RemoteRefusal> {
    let bad = |why: &str| RemoteRefusal::BadReport(why.to_owned());
    let mut lines = report.lines();
    let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    let column = |name: &str| header.iter().position(|h| *h == name);
    let (Some(run_at), Some(fastq_at)) = (column("run_accession"), column("fastq_ftp")) else {
        return Err(if report.trim().is_empty() {
            RemoteRefusal::NotFound(accession.to_owned())
        } else {
            bad("no run_accession and fastq_ftp columns")
        });
    };
    let model_at = column("instrument_model");
    let mut rows: Vec<(&str, &str, &str)> = Vec::new();
    for line in lines.take(REPORT_ROWS_MAX + 1) {
        let fields: Vec<&str> = line.split('\t').collect();
        let (Some(run), fastq) = (fields.get(run_at), fields.get(fastq_at)) else {
            return Err(bad("a row without a run"));
        };
        let model = model_at.and_then(|at| fields.get(at)).copied().unwrap_or("");
        rows.push((run, fastq.copied().unwrap_or(""), model));
    }
    if rows.len() > REPORT_ROWS_MAX {
        return Err(bad("more runs than an experiment has"));
    }
    rows.sort_unstable();
    let Some((run, fastq, model)) = rows.first().copied() else {
        return Err(RemoteRefusal::NotFound(accession.to_owned()));
    };
    if !is_accession(run) {
        return Err(bad("a run accession that is not one"));
    }
    let urls =
        fastq.split(';').filter(|f| !f.is_empty()).map(url).collect::<Result<Vec<_>, _>>()?;
    let urls = numbered_only(urls);
    if urls.is_empty() {
        return Err(RemoteRefusal::NoFastq(run.to_owned()));
    }
    if urls.len() > FILES_MAX {
        return Err(RemoteRefusal::TooManyFiles { run: run.to_owned(), files: urls.len() });
    }
    assert!(rows.len() <= REPORT_ROWS_MAX, "rows are capped");
    assert!(!urls.is_empty() && urls.len() <= FILES_MAX, "one to a lane's files");
    // A model is short printable text, or it is not kept: it goes into a report line.
    let plain = model.chars().all(|c| c.is_ascii_graphic() || c == ' ')
        && model.len() <= INSTRUMENT_MODEL_BYTES_MAX;
    let instrument_model = Some(model.trim().to_owned()).filter(|m| plain && !m.is_empty());
    Ok(Run { run: run.to_owned(), urls, instrument_model })
}

/// ENA lists a paired run as `_1`, `_2` (and `_3`, `_4` for index reads), and sometimes also an
/// unnumbered file of reads whose mate was lost; that file is not one of the lane's reads.
fn numbered_only(urls: Vec<String>) -> Vec<String> {
    let count = urls.len();
    let numbered = |u: &String| {
        let stem = u.trim_end_matches(".gz").trim_end_matches(".fastq").trim_end_matches(".fq");
        stem.rsplit_once('_')
            .is_some_and(|(_, n)| n.len() == 1 && n.bytes().all(|b| b.is_ascii_digit()))
    };
    if !urls.iter().any(numbered) {
        return urls;
    }
    let kept: Vec<String> = urls.into_iter().filter(numbered).collect();
    assert!(!kept.is_empty(), "a numbered file was listed, so one is kept");
    assert!(kept.len() <= count, "the filter only removes files");
    kept
}

/// A `fastq_ftp` entry as an HTTPS URL, only on ENA's host and only with a plain path.
fn url(entry: &str) -> Result<String, RemoteRefusal> {
    let plain = entry.bytes().all(|b| b.is_ascii_alphanumeric() || b"._/-".contains(&b));
    if !entry.starts_with(FASTQ_HOST) || !plain || entry.contains("..") {
        return Err(RemoteRefusal::BadReport(format!("a file outside ENA's FASTQ host: {entry}")));
    }
    Ok(format!("https://{entry}"))
}

/// The first `limit` reads of the file at `url`, streamed, and the bytes it took; the download
/// stops once they are in.
pub fn stream(url: &str, limit: usize) -> Result<(PathBuf, Reads, u64), RemoteRefusal> {
    assert!(url.starts_with("https://") && url.contains(FASTQ_HOST), "{url} was checked");
    let mut child = spawn(curl(url, STREAM_SECONDS_MAX))?;
    let Some(stdout) = child.stdout.take() else { unreachable!("curl's stdout is piped") };
    let name = PathBuf::from(url);
    let mut stdout = Capped::new(stdout, STREAM_BYTES_MAX);
    let result = Reads::from_reader(&mut stdout, &name, limit);
    // Enough read: stop the download (the pipe still open) rather than wait for the whole file.
    let ended = stop(&mut child);
    let (bytes, capped) = (stdout.passed, stdout.reached());
    let what = url.to_owned();
    match (result, ended) {
        // curl failed before it was stopped: what came through may be cut short.
        (_, Ended::Failed(reason)) => Err(RemoteRefusal::Download { what, reason }),
        (Ok(reads), _) if !reads.is_empty() => {
            assert!(reads.len() <= limit, "a sample is capped at the limit");
            assert!(bytes > 0, "reads came from bytes that passed the cap");
            Ok((name, reads, bytes))
        }
        (Err(_), _) if capped => Err(RemoteRefusal::StreamCap { url: what, bytes }),
        (Ok(_), _) => Err(RemoteRefusal::Download { what, reason: "no reads".into() }),
        (Err(source), _) => Err(RemoteRefusal::Read { url: what, source }),
    }
}

/// A short name for a streamed file in the report: the file name of its URL.
pub fn short_name(url: &Path) -> &str {
    let name = url.file_name().and_then(|n| n.to_str()).unwrap_or("");
    assert!(!name.is_empty(), "a FASTQ URL ends in a file name");
    assert!(!name.contains('/'), "a file name, not a path");
    name
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A `fastq_ftp` entry becomes an HTTPS URL on ENA's FASTQ host, and nothing else does.
    #[test]
    fn urls_stay_on_the_host() {
        let good = format!("{FASTQ_HOST}SRR123/001/SRR123001/SRR123001_1.fastq.gz");
        assert_eq!(url(&good).ok(), Some(format!("https://{good}")));
        for bad in [
            "evil.example.org/vol1/x.fastq.gz".to_string(),
            format!("{FASTQ_HOST}../../etc/passwd"),
            format!("{FASTQ_HOST}a b.fastq.gz"),
            format!("{FASTQ_HOST}x.fastq.gz?q=1"),
            format!("http://{FASTQ_HOST}x.fastq.gz"),
        ] {
            assert!(url(&bad).is_err(), "{bad} is refused");
        }
    }

    /// Numbered files are kept when any are numbered; otherwise every file is.
    #[test]
    fn unnumbered_files_go_when_numbered_ones_exist() {
        let files = |names: &[&str]| names.iter().map(|n| (*n).to_string()).collect::<Vec<_>>();
        let kept = numbered_only(files(&["a/R.fastq.gz", "a/R_1.fastq.gz", "a/R_2.fastq.gz"]));
        assert_eq!(kept, files(&["a/R_1.fastq.gz", "a/R_2.fastq.gz"]));
        assert_eq!(numbered_only(files(&["a/R.fastq.gz"])), files(&["a/R.fastq.gz"]));
    }

    /// A ~/.curlrc with `output = <file>` would send the ENA report to that file, and the run
    /// would be refused as unknown, so curl reads no config file.
    #[test]
    fn curl_skips_the_users_config_and_bounds_redirects() {
        let command = curl("https://www.ebi.ac.uk/x", "5");
        let args: Vec<_> = command.get_args().filter_map(|a| a.to_str()).collect();
        assert_eq!(args.first(), Some(&"--disable"), "curl takes --disable only first: {args:?}");
        assert!(args.windows(2).any(|w| w == ["--max-redirs", REDIRECTS_MAX]), "{args:?}");
        assert!(args.windows(2).any(|w| w == ["--proto-redir", "=https"]), "{args:?}");
    }

    /// The cap: an exact fit is told from an over-size stream by one byte.
    #[test]
    fn capped_reader_stops_and_says_so() {
        let mut fits = Capped::new(&b"abcd"[..], 5);
        let mut out = Vec::new();
        assert_eq!(fits.read_to_end(&mut out).ok(), Some(4));
        assert!(!fits.reached(), "4 bytes under a 5 byte cap");
        let mut over = Capped::new(&b"abcdefgh"[..], 5);
        out.clear();
        assert_eq!(over.read_to_end(&mut out).ok(), Some(5));
        assert!(over.reached(), "cut at the cap");
        assert_eq!(out, b"abcde");
    }

    #[test]
    fn accessions_and_not_accessions() {
        for yes in ["SRR123456", "ERR1234567890", "DRR983898", "SRX15432977"] {
            assert!(is_accession(yes), "{yes}");
        }
        for no in [
            "SRR12345",
            "SRR12345678901",
            "SRR12a456",
            "srr123456",
            "GSM123456",
            "SRR123456;rm",
            "R1.fastq.gz",
            "",
        ] {
            assert!(!is_accession(no), "{no}");
        }
    }

    #[test]
    fn reports_parse_strictly() {
        let ok = "run_accession\tfastq_ftp\n\
                  SRR1000002\tftp.sra.ebi.ac.uk/vol1/fastq/a_1.fastq.gz;\
                  ftp.sra.ebi.ac.uk/vol1/fastq/a_2.fastq.gz\n\
                  SRR1000001\tftp.sra.ebi.ac.uk/vol1/x.fastq.gz\n";
        let found = parse_report(ok, "SRX1000001").unwrap_or_else(|e| unreachable!("{e}"));
        assert_eq!(found.run, "SRR1000001", "the first run by accession");
        assert_eq!(found.urls, ["https://ftp.sra.ebi.ac.uk/vol1/x.fastq.gz"]);
        assert_eq!(found.instrument_model, None, "no instrument_model column, none kept");
        let modelled = "run_accession\tfastq_ftp\tinstrument_model\n\
                        SRR1000001\tftp.sra.ebi.ac.uk/vol1/x.fastq.gz\tIllumina NovaSeq X\n";
        let found = parse_report(modelled, "SRR1000001").unwrap_or_else(|e| unreachable!("{e}"));
        assert_eq!(found.instrument_model.as_deref(), Some("Illumina NovaSeq X"));
        let evil = "run_accession\tfastq_ftp\nSRR111111\tevil.example/vol1/x.fastq.gz\n";
        assert!(matches!(parse_report(evil, "SRR111111"), Err(RemoteRefusal::BadReport(_))));
        let climb = "run_accession\tfastq_ftp\nSRR111111\tftp.sra.ebi.ac.uk/vol1/../x\n";
        assert!(matches!(parse_report(climb, "SRR111111"), Err(RemoteRefusal::BadReport(_))));
        let unpaired = "run_accession\tfastq_ftp\nSRR111111\tftp.sra.ebi.ac.uk/vol1/r.fastq.gz;\
                        ftp.sra.ebi.ac.uk/vol1/r_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/r_2.fastq.gz\n";
        let found = parse_report(unpaired, "SRR111111").unwrap_or_else(|e| unreachable!("{e}"));
        assert_eq!(found.urls.len(), 2, "the unnumbered file is set aside: {:?}", found.urls);
        let none = "run_accession\tfastq_ftp\nSRR111111\t\n";
        assert!(matches!(parse_report(none, "SRR111111"), Err(RemoteRefusal::NoFastq(_))));
        assert!(matches!(parse_report("", "SRR111111"), Err(RemoteRefusal::NotFound(_))));
    }
}
