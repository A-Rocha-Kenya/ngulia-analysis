# Ngulia analyses

Scientific analyses of the curated Ngulia ringing dataset, organized by question. See [the analysis guide](docs/analysis.md) for scope, run order, and interpretation limits.

The canonical dataset and its Zenodo and GBIF publication workflows are maintained in [ngulia-dataset](https://github.com/A-Rocha-Kenya/ngulia-dataset). These analyses read the four curated CSV files from that project. Set `NGULIA_DATASET_DIR` to the local dataset checkout before running scripts; otherwise the scripts expect a sibling directory named `ngulia-dataset`.

For reproducible results, record the dataset's Zenodo version DOI and Git release tag with each analysis. The dataset concept DOI is [10.5281/zenodo.21395879](https://doi.org/10.5281/zenodo.21395879); it points to the latest version rather than identifying a fixed input version.

Scripts run from this repository root and write results under `outputs/analysis/`. Generated outputs are kept locally and are not tracked in Git.
