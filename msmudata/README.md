# msmudata — Hierarchical Feature Index for MuData

Related [GitHub Issue](https://github.com/scverse/mudata/issues/111) | [alphabase scverse tutorial](https://alphabase.readthedocs.io/en/latest/tutorials/tutorial_scverse_compatibility.html)

## Problem

In LC/MS proteomics, quantification is **hierarchical**: mass spectrometers measure *precursors* (~100k features), which search engines aggregate into *proteins* (~3–10k features). The relationship is **N:M** — many precursors map to one protein, and shared peptides can map a single precursor to multiple proteins.

[MuData](https://mudata.readthedocs.io/) can store these modalities side-by-side, but provides **no way to query the mapping between them**. You cannot currently ask:
- "Which precursors belong to protein P04637?"
- "Plot the intensity distribution of all precursors that map to protein MAPK"

## Solution

This package introduces two classes that add a **queryable hierarchical feature index** on top of MuData:

| Class | Role |
|---|---|
| `FeatureMapping` | Sparse, bidirectional N:M graph between two feature sets (e.g. precursors ↔ proteins). Backed by a DataFrame of edges + a scipy sparse adjacency matrix. |
| `HierarchicalMuData` | Wraps a standard `MuData` and attaches one or more `FeatureMapping`s. Provides cross-modality feature lookups and data slicing. |

## Quickstart

```python
from src import HierarchicalMuData, make_synthetic_proteomics_data

# Generate synthetic data (6 samples, 200 proteins, ~1800 precursors, N:M mapping)
hmdata = make_synthetic_proteomics_data()
print(hmdata)
# HierarchicalMuData with 6 obs
#   modalities: proteins: 6 x 200, precursors: 6 x 1827
#   mappings: precursors→proteins (1918 edges)

# Forward query: which proteins does a precursor map to?
hmdata.get_related_features("precursors", "proteins", "prec_000000")
# ['PROT_0000']

# Reverse query: which precursors map to a protein?
hmdata.get_related_features("proteins", "precursors", "PROT_0042")
# ['prec_000394', 'prec_000395', ..., 'prec_000405']

# Data slicing: get precursor-level AnnData for a specific protein
adata_slice = hmdata.slice_data("precursors", by="proteins", features="PROT_0042")
# AnnData view: 6 samples × 12 precursors

# Multi-protein query (e.g. a pathway)
hmdata.slice_data("precursors", by="proteins", features=["PROT_0001", "PROT_0010", "PROT_0042"])

# Serialize mappings into MuData.uns (survives disk I/O)
hmdata.store_mappings_in_uns()

# Reconstruct from a saved MuData
hmdata_restored = HierarchicalMuData.from_mudata(hmdata.mdata)
```

## Project structure

```
msmudata/
├── README.md
├── src/
│   ├── __init__.py
│   ├── feature_mapping.py      # FeatureMapping class
│   ├── hierarchical_mudata.py  # HierarchicalMuData class
│   └── synthetic.py            # Synthetic data generator
├── notebooks/
│   └── demo_hierarchical_mudata.ipynb   # Interactive demo
└── rfc/
    └── .gitkeep
```

## Setup

```bash
# Create conda environment
conda create -n msmudata python=3.12 -y
conda activate msmudata
pip install anndata mudata numpy pandas scipy matplotlib scanpy ipykernel

# Register Jupyter kernel
python -m ipykernel install --user --name msmudata --display-name "msmudata"
```

Then open `notebooks/demo_hierarchical_mudata.ipynb` and select the **msmudata** kernel.

## API overview

### `FeatureMapping`

| Method | Description |
|---|---|
| `targets_of(source)` | All targets linked to a source feature |
| `sources_of(target)` | All sources linked to a target feature |
| `targets_of_many(sources)` | Union of targets for multiple sources |
| `sources_of_many(targets)` | Union of sources for multiple targets |
| `subgraph(sources, targets)` | New mapping restricted to a subset |
| `from_long_dataframe(df, ...)` | Construct from a PSM-report-style DataFrame |
| `.adjacency` | Sparse boolean matrix (n_sources × n_targets) |

### `HierarchicalMuData`

| Method | Description |
|---|---|
| `get_related_features(from_mod, to_mod, features)` | Cross-modality feature lookup (auto-resolves direction) |
| `slice_data(data_mod, by=mod, features=...)` | Slice an AnnData to features related to feature(s) in another modality |
| `add_mapping(source_mod, target_mod, mapping)` | Register a new mapping (validates against var_names) |
| `get_mapping(source_mod, target_mod)` | Retrieve a registered mapping |
| `mapping_summary()` | Summary table of all mappings |
| `store_mappings_in_uns()` | Persist mappings into `MuData.uns` |
| `from_mudata(mdata)` | Reconstruct from a MuData with mappings in `.uns` |

## Roadmap

- [ ] Implement an RFC (see `rfc/RFC.md`)
- [x] Implement the prototype data structure (scverse/anndata/mudata compatible)
- [x] Downstream analysis demo (plot precursor intensities for a protein)
- [ ] Data ingestion reader from a search engine output (e.g. alphadia, DIA-NN)
- [ ] Support deeper hierarchies (fragments → precursors → peptides → proteins → genes)
- [ ] Disk I/O roundtrip tests (`.h5mu` write/read)
- [ ] Discussion with mudata maintainers on upstream integration

## References

- [mudata docs](https://mudata.readthedocs.io/stable/notebooks/nuances.html)
- [QFeatures](https://rformassspectrometry.github.io/QFeatures/articles/QFeatures.html) — conceptually similar R package
- [alphaquant](https://github.com/MannLabs/alphaquant.git) — potential application of the mapping approach
- [alphabase scverse tutorial](https://alphabase.readthedocs.io/en/latest/tutorials/tutorial_scverse_compatibility.html) — PSM reader → AnnData/MuData
