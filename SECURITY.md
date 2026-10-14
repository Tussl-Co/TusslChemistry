# Security

tchem reads files you give it, and with an accession it fetches reads from the public archive ENA.
It treats every input as untrusted. A malformed file, name or stream gets a typed refusal, never a
crash. Every read has a size limit (CONTRIBUTING.md, "Engineering rules"). The readers and
detection are fuzzed (`fuzz/`).

## Reporting a vulnerability

Please report it privately through GitHub's private vulnerability reporting: the repository's
**Security** tab, then **Report a vulnerability**. Do not open a public issue for it. Include the
tchem version (`tchem --version`), what you ran, and an input that shows it if you can share one.

We aim to acknowledge a report within 5 working days and to agree a fix and disclosure date with
you.

## Supported versions

Fixes go into the latest release. tchem is at 0.x: there are no long-term support branches yet.
