# CHANGELOG

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project adheres to [Semantic Versioning](https://semver.org/).

## [1.1.0] - 2025-09-25

### Added

- Add TSS (Transcription Start Site) QC plot generation to assess signal around TSS regions.
- Parallelize footprint calling (fixed to 10 cores) to reduce runtime on multi-core systems and improve throughput on large datasets.

### Fixed

- Pin `igvtools` to version `2.16.2` (downgraded from `2.17.3`) to resolve TDF-generation issues observed with `2.17.3`.
