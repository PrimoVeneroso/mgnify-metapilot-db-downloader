# mgnify-metapilot-db

Download an [MGnify](https://www.ebi.ac.uk/metagenomics/browse/genomes) genome catalogue
(default: **human-gut v3.0**, 6,461 species) and turn it into a folder that
[MetaPilot](https://github.com/cksakura/MetaPilot-Releases) accepts in **Add MAG database**.

Also validates the result against MetaPilot's import rules, so you know *before*
importing whether something is missing or corrupted.

## Why

MetaPilot expects plain `.faa` and `_eggNOG.tsv` files in two flat folders plus
`genomes-all_metadata.tsv`. Since v3.0, EBI ships these files gzip-compressed inside
thousands of nested folders, and species IDs are not sequential (representatives were
replaced, see `replaced_species_reps.tsv` on EBI). Crawling the FTP index folder by
folder takes ~a day; this tool reads the species list from the metadata file, builds
the exact URLs and downloads in parallel (~15–30 min on a good connection).

## Output

```text
MGnify_MAG_DB/human-gut/v3.0/
├── genomes-all_metadata.tsv         # needs columns Genome and Lineage
├── original_db/<genomeID>.faa       # 6,461 protein FASTA files
└── eggNOG/<genomeID>_eggNOG.tsv     # 6,461 eggNOG-mapper annotations
```

Disk space: ~5.5 GB downloaded, ~15 GB after decompression.
`Marker4.fasta` is **not** created here: MetaPilot generates it when you add the folder.

## Usage

Two equivalent implementations; use whichever fits your system.

| Script | Platforms | Requirements |
| --- | --- | --- |
| `mgnify_metapilot_db.py` | Windows, macOS, Linux | Python ≥ 3.8 (standard library only) |
| `download_mgnify_mag_db.sh` | Linux, macOS, WSL, HPC | bash, curl, gzip, awk, GNU find |

```bash
# Quick test with 2 species
python mgnify_metapilot_db.py -o ./MGnify_TEST -n 2

# Full download + validation
python mgnify_metapilot_db.py -o ./MGnify_MAG_DB -j 8

# Validate an existing folder only (no network)
python mgnify_metapilot_db.py -o ./MGnify_MAG_DB -V
```

The Bash script takes the same options (`./download_mgnify_mag_db.sh -o DIR -j 8`, `-n 2`, `-V`).

| Option | Meaning | Default |
| --- | --- | --- |
| `-o` | output root folder | `MGnify_MAG_DB` |
| `-c` | catalogue (e.g. `human-gut`, `mouse-gut`) | `human-gut` |
| `-v` | catalogue version | `v3.0` |
| `-j` | parallel downloads (keep ≤ 16 to be polite to EBI) | `8` |
| `-n` | only the first N species (test run) | all |
| `-V` | validate only, do not download | off |

Interrupted? Just run the same command again: complete files are skipped.

## Validation checks

| Rule | Check |
| --- | --- |
| 1 | same number of files in `original_db/` and `eggNOG/` |
| 2 | no empty files, unexpected extensions or leftover `.part`/`.gz` |
| 3 | every genome has both files; IDs match exactly the `Species_rep` list in the metadata |
| 4 | FASTA headers are `>genomeID_proteinID protein name`, sequences present and valid |
| 5 | eggNOG files are tab-separated, have data rows, consistent column count, queries start with `genomeID_` |
| 6 | `genomes-all_metadata.tsv` has `Genome` and `Lineage` filled on every row |

Passing these checks does not guarantee MetaPilot will accept the folder; the final
test is importing it with **Add MAG database**.

## Tested

- human-gut v3.0: 6,461 / 6,461 species downloaded, all checks passed (Linux, Nov 2026).
- Other catalogues/versions use the same EBI layout but have not been tested.
  v2.0.2 ships uncompressed files and is **not** supported by the current code.

## Data license and citation

This tool only downloads public data; it does not redistribute it. MGnify data are
provided by EMBL-EBI under their [terms of use](https://www.ebi.ac.uk/about/terms-of-use).
If you use the catalogue, cite MGnify and the catalogue release (see `README_v3.0.txt`
on the EBI FTP).

## License

[MIT](LICENSE) © 2026 Primo Veneroso, for the code in this repository.
