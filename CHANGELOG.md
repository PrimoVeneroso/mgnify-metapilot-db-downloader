# Changelog

## 1.0.0 — 2026-11-07

- Bash and Python implementations with the same options (`-o -c -v -j -n -V`).
- Species list read from `Species_rep` in `genomes-all_metadata.tsv` (includes IDs with `.1` suffix).
- Parallel download, resumable, gzip integrity check before installing files.
- Validation of MetaPilot import rules 1–6 and exact ID comparison with the catalogue.
- Tested on human-gut v3.0 (6,461 species).
