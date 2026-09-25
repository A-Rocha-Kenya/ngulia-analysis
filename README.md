# Ngulia analyses

Scientific analyses of the curated Ngulia ringing dataset, organized by question. See [the analysis guide](docs/analysis.md) for scope, run order, and interpretation limits.

The canonical dataset and its Zenodo and GBIF publication workflows are maintained in [ngulia-dataset](https://github.com/A-Rocha-Kenya/ngulia-dataset). These analyses read the four curated CSV files and the shared capture-group rule from that project. Set `NGULIA_DATASET_DIR` to the local dataset checkout before running scripts; otherwise the scripts expect a sibling directory named `ngulia-dataset`.

For reproducible results, record the dataset's Zenodo version DOI and Git release tag with each analysis once the dataset is published. Until then, record the exact dataset Git commit and curated CSV checksums.

Scripts run from this repository root and write results under `outputs/analysis/`. Generated outputs are kept locally and are not tracked in Git.
