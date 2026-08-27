# data/raw — never edited, checksummed on download

Raw source files are **not** committed (see `.gitignore`). Drop each source into
its own subdirectory; loaders in `src/warhead/io/` read from here and raise an
actionable error naming the expected file if it is absent.

## Expected files (G1 -> G2b slice)

### `depmap/`  (depmap.org/portal/download, release 26Q1)
- `Model.csv` — model metadata; must carry a doubling-time column
  (`doubling_time_hours` / `DoublingTime` / …). If doubling time ships separately,
  add `DoublingTime.csv` with `ModelID` + a doubling-time column.
- `OmicsExpressionProteinCodingGenesTPMLogp1.csv` — expression (G2a: ABCB1/ABCG2).
- `CRISPRGeneEffect.csv` — Chronos dependency (G2c).

### `prism/`  (DepMap portal, PRISM Repurposing **secondary** screen)
- `prism_secondary_long.parquet` (or `.csv`) — tidy long dose-response with
  columns `[compound_id, ModelID, dose_M, viability]`. Only the secondary screen
  (8-pt dose-response) is usable for G1; the primary screen is single-dose.

## Other sources (deferred gates)
See `WARHEAD.md` §2 and the docstring of each loader in `src/warhead/io/` for the
source URL, access route, and target tidy schema (CTRP, GDSC, NCI-60, Tahoe,
LINCS, JUMP, ADCdb, ChEMBL, COCONUT, NPAtlas, GTEx, HPA, FAERS, TCGA).

## Integrity
`warhead.io.download_file(url, dest, sha256=...)` verifies a checksum on download
and refuses a rotated/corrupt file. Record checksums when you pull a release.
