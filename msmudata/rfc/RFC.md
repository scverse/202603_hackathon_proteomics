# RFC: Hierarchical Feature Mappings for MuData

**Status**: Draft  
**Authors**: msmudata hackathon team  
**Related**: [scverse/mudata#111](https://github.com/scverse/mudata/issues/111), [alphabase scverse tutorial](https://alphabase.readthedocs.io/en/latest/tutorials/tutorial_scverse_compatibility.html)  
**Prototype**: [`src/`](../src/)

---

## 1. Motivation

### 1.1 The problem

In LC/MS proteomics, quantification is inherently hierarchical:

```
fragments → precursors → peptides → proteins → genes
```

Each level has its own feature space and intensity matrix. The relationships between levels are **N:M** — a precursor may map to multiple proteins (shared peptides / protein groups), and a protein is quantified from many precursors.

Today, [MuData](https://mudata.readthedocs.io/) can store each level as a separate modality:

```python
mdata = md.MuData({"proteins": protein_adata, "precursors": precursor_adata})
```

But MuData provides **no mechanism to express or query the mapping between their feature indices**. This means common analytical tasks like "plot all precursor intensities for protein X" require users to maintain and join external lookup tables — error-prone, not portable, and not standardized.

### 1.2 Prior art

| Project | Approach | Limitation |
|---|---|---|
| [QFeatures](https://rformassspectrometry.github.io/QFeatures/articles/QFeatures.html) (R/Bioconductor) | `AssayLinks` connecting `SummarizedExperiment` layers via a DAG | R-only, not scverse-compatible |
| [alphaquant](https://github.com/MannLabs/alphaquant) | Internal mapping tables (DataFrames) | Not formalized, not portable |
| [alphabase scverse tutorial](https://alphabase.readthedocs.io/en/latest/tutorials/tutorial_scverse_compatibility.html) | Reads PSM reports into AnnData/MuData | No cross-modality mapping |
| mudata#111 proposals by @mffrank | DAG mapping proof-of-concept | Not integrated |
| mudata#111 proposals by @ilan-gold | Pandas extension index arrays | Design-stage |

### 1.3 Scope

This RFC proposes a **minimal, composable extension** to MuData for representing directed, queryable N:M feature mappings between modalities. The design is motivated by MS proteomics but generalizes to any multi-modal setting with hierarchical or relational feature structure (e.g. ATAC peaks → genes, spatial transcriptomics spots → cells).

---

## 2. Design

### 2.1 Core concepts

**`FeatureMapping`** — a directed bipartite graph between two ordered feature sets.

- **source_index**: `pd.Index` — feature identifiers in the source modality (e.g. precursor IDs)
- **target_index**: `pd.Index` — feature identifiers in the target modality (e.g. protein IDs)
- **edges**: `pd.DataFrame` with at minimum columns `["source", "target"]`. Additional columns (e.g. `"weight"`, `"score"`) are preserved as edge attributes.
- **adjacency**: lazily-computed `scipy.sparse.csr_matrix` of shape `(n_sources, n_targets)` for efficient queries.

**`HierarchicalMuData`** — a wrapper around `MuData` that holds one or more `FeatureMapping`s, keyed by `(source_mod, target_mod)` tuples.

```
┌──────────────────────────────────────────────────┐
│  HierarchicalMuData                              │
│                                                  │
│  ┌─────────────┐         ┌─────────────┐         │
│  │  AnnData    │         │  AnnData    │         │
│  │ "precursors"│         │ "proteins"  │         │
│  │  6 × 1827   │         │  6 × 200    │         │
│  └──────┬──────┘         └──────┬──────┘         │
│         │    FeatureMapping     │                 │
│         │  (1918 N:M edges)     │                 │
│         └───────────────────────┘                 │
│                                                  │
│  mdata.uns["mapping_precursors__proteins"] = {   │
│    source_index: [...],                          │
│    target_index: [...],                          │
│    edges: {source: [...], target: [...]},        │
│  }                                               │
└──────────────────────────────────────────────────┘
```

### 2.2 Relationship to MuData

`HierarchicalMuData` does **not** subclass `MuData`. It wraps it:

```python
class HierarchicalMuData:
    _mdata: md.MuData
    _mappings: dict[tuple[str, str], FeatureMapping]
```

**Rationale**: Subclassing `MuData` would couple tightly to its internal API and complicate serialization. A wrapper keeps the mapping layer cleanly separated and allows the prototype to evolve independently. If this design proves useful, the mapping infrastructure could eventually be proposed as a first-class MuData feature.

### 2.3 Direction and auto-resolution

Mappings are stored with an explicit direction: `(source_mod, target_mod)`. However, all query methods **auto-resolve the reverse direction**. If only `("precursors", "proteins")` is registered, both of these work:

```python
hmdata.get_related_features("precursors", "proteins", "prec_000042")  # forward
hmdata.get_related_features("proteins", "precursors", "PROT_0010")    # reverse (auto-resolved)
```

This is possible because the sparse adjacency matrix supports efficient lookups in both row and column directions.

### 2.4 Validation

When a `FeatureMapping` is registered on a `HierarchicalMuData`, the following invariants are enforced:

1. Both `source_mod` and `target_mod` must exist in `mdata.mod`.
2. Every feature in `mapping.source_index` must exist in `mdata.mod[source_mod].var_names`.
3. Every feature in `mapping.target_index` must exist in `mdata.mod[target_mod].var_names`.

The mapping indices are allowed to be a **subset** of the modality's `var_names` (not every feature needs to participate in a mapping).

---

## 3. API

### 3.1 `FeatureMapping`

#### Construction

```python
# From explicit indices + edges
mapping = FeatureMapping(
    source_index=["prec_1", "prec_2", "prec_3"],
    target_index=["prot_A", "prot_B"],
    edges=pd.DataFrame({"source": ["prec_1", "prec_2", "prec_3"],
                         "target": ["prot_A", "prot_A", "prot_B"]}),
)

# From a long-format DataFrame (e.g. PSM report columns)
mapping = FeatureMapping.from_long_dataframe(
    psm_report,
    source_col="precursor_id",
    target_col="protein_id",
)
```

#### Queries

| Method | Returns | Description |
|---|---|---|
| `targets_of(source)` | `list[str]` | All targets linked to one source |
| `sources_of(target)` | `list[str]` | All sources linked to one target |
| `targets_of_many(sources)` | `list[str]` | Union of targets for multiple sources |
| `sources_of_many(targets)` | `list[str]` | Union of sources for multiple targets |
| `subgraph(sources?, targets?)` | `FeatureMapping` | Restricted sub-mapping |

#### Properties

| Property | Type | Description |
|---|---|---|
| `.source_index` | `pd.Index` | Source feature identifiers |
| `.target_index` | `pd.Index` | Target feature identifiers |
| `.edges` | `pd.DataFrame` | Edge list with optional attributes |
| `.adjacency` | `csr_matrix` | Sparse boolean matrix (lazy) |
| `.n_sources`, `.n_targets`, `.n_edges` | `int` | Counts |

### 3.2 `HierarchicalMuData`

#### Construction

```python
hmdata = HierarchicalMuData(
    mdata=mdata,
    mappings={("precursors", "proteins"): mapping},
)

# Or add mappings incrementally
hmdata = HierarchicalMuData(mdata)
hmdata.add_mapping("precursors", "proteins", mapping)
```

#### Cross-modality queries

```python
# Feature lookup (works in both directions)
hmdata.get_related_features("precursors", "proteins", "prec_000042")
hmdata.get_related_features("proteins", "precursors", "PROT_0010")
hmdata.get_related_features("proteins", "precursors", ["PROT_0001", "PROT_0002"])

# Data slicing — returns an AnnData view
adata_slice = hmdata.slice_data("precursors", by="proteins", features="PROT_0010")
# → AnnData with only the precursor columns that map to PROT_0010
```

#### Serialization

```python
# Persist to MuData.uns (survives .h5mu write/read)
hmdata.store_mappings_in_uns()

# Reconstruct from a MuData loaded from disk
hmdata = HierarchicalMuData.from_mudata(mdata)
```

---

## 4. Storage format

### 4.1 In-memory

Mappings are held as `FeatureMapping` instances in a Python `dict` keyed by `(source_mod, target_mod)`.

### 4.2 On-disk (via `MuData.uns`)

Each mapping is stored as a dictionary in `mdata.uns` under the key `mapping_{source_mod}__{target_mod}`:

```python
mdata.uns["mapping_precursors__proteins"] = {
    "source_mod": "precursors",          # str
    "target_mod": "proteins",            # str
    "source_index": ["prec_0", ...],     # list[str]
    "target_index": ["PROT_0", ...],     # list[str]
    "edges": {                           # dict of lists (DataFrame-like)
        "source": ["prec_0", ...],
        "target": ["PROT_0", ...],
        # optional additional columns:
        # "weight": [0.95, ...],
    },
}
```

This format is chosen because:

1. `MuData.uns` accepts nested dicts and lists natively.
2. It serializes to `.h5mu` / HDF5 without custom infrastructure.
3. It is human-inspectable.
4. The `mapping_` prefix provides a namespace to avoid collisions.

### 4.3 Future considerations

- **Sparse matrix storage**: For very large mappings (>1M edges), storing the sparse adjacency directly in HDF5 may be more space-efficient than the edge list. This would use the same `indptr`/`indices`/`data` format that AnnData already uses for sparse `X`.
- **Dedicated `.h5mu` group**: If mappings become a first-class MuData concept, they could get their own HDF5 group (e.g. `/mappings/precursors__proteins/`) rather than living in `uns`.

---

## 5. Design decisions and alternatives

### 5.1 Wrapper vs. subclass

**Chosen**: Wrapper (`HierarchicalMuData` holds a `MuData`).  
**Alternative**: Subclass `MuData`.

A wrapper was chosen because:
- MuData's internal API is not designed for extension (view mechanics, update logic).
- A wrapper makes the mapping layer opt-in and composable.
- It avoids breakage when MuData's internals change.

If this pattern proves useful, the long-term path is to propose mappings as a core MuData feature (akin to how `obsm`/`varm` were added to AnnData).

### 5.2 Sparse matrix vs. DataFrame-only

**Chosen**: Both — DataFrame for construction/inspection, sparse matrix for queries.  
**Alternative**: DataFrame-only, or sparse-only.

The dual representation gives the best of both worlds:
- DataFrames are easy to construct from PSM reports and can carry edge attributes.
- Sparse matrices enable O(nnz) lookups via row/column slicing.
- The sparse matrix is lazily computed and cached, so there is no cost if only DataFrame-based operations are used.

### 5.3 Edge attributes

The `edges` DataFrame can carry additional columns beyond `source` and `target` (e.g. identification score, spectral angle). These are preserved through construction and subgraph operations. The sparse adjacency matrix is boolean-only (edge presence). Weighted adjacency could be added if use cases emerge.

### 5.4 Mapping direction

Mappings have a canonical direction (source → target) but queries work in both directions. An alternative would be to store two separate mappings for each relationship, but this wastes memory and creates synchronization issues. The single-mapping + auto-resolve approach is simpler.

---

## 6. Example: end-to-end workflow

```python
import anndata as ad
import mudata as md
import pandas as pd
from src import FeatureMapping, HierarchicalMuData

# --- Step 1: Read PSM report from a search engine ---
psm_report = pd.read_parquet("alphadia_output.parquet")

# --- Step 2: Build AnnData for each level ---
protein_adata = build_protein_adata(psm_report)   # 6 × 9880
precursor_adata = build_precursor_adata(psm_report)  # 6 × 117615

# --- Step 3: Create MuData ---
mdata = md.MuData({"proteins": protein_adata, "precursors": precursor_adata})

# --- Step 4: Extract the mapping from the PSM report ---
mapping = FeatureMapping.from_long_dataframe(
    psm_report,
    source_col="precursor_id",
    target_col="protein_group",
)

# --- Step 5: Wrap in HierarchicalMuData ---
hmdata = HierarchicalMuData(mdata, mappings={("precursors", "proteins"): mapping})

# --- Step 6: Query across the hierarchy ---
# "Which precursors contribute to protein P04637?"
precs = hmdata.get_related_features("proteins", "precursors", "P04637")

# "Give me the intensity matrix for those precursors"
adata_slice = hmdata.slice_data("precursors", by="proteins", features="P04637")
# → AnnData view: 6 samples × N precursors

# --- Step 7: Save (mappings survive disk I/O) ---
hmdata.store_mappings_in_uns()
mdata.write("experiment.h5mu")

# --- Step 8: Load and reconstruct ---
mdata_loaded = md.read("experiment.h5mu")
hmdata_loaded = HierarchicalMuData.from_mudata(mdata_loaded)
```

---

## 7. Roadmap

### Implemented (this prototype)

- [x] `FeatureMapping` with bidirectional N:M queries, subgraph extraction, and `from_long_dataframe` constructor
- [x] `HierarchicalMuData` with cross-modality feature lookups, data slicing, and `MuData.uns` serialization
- [x] Synthetic data generator for testing
- [x] Demo notebook with intensity distribution plots

### Next steps

- [ ] **Data ingestion**: Reader that constructs `HierarchicalMuData` from an alphadia / DIA-NN / QuantMS PSM report in a single call
- [ ] **Disk I/O roundtrip**: Verify `.h5mu` write → read → `from_mudata()` preserves mappings exactly
- [ ] **Deeper hierarchies**: Chain multiple `FeatureMapping`s (e.g. precursors → peptides → proteins) and support transitive queries (`hmdata.get_related_features("precursors", "proteins", ...)` when only precursor→peptide and peptide→protein mappings are registered)
- [ ] **Integration with scanpy/scverse tools**: Ensure `slice_data` views compose well with scanpy preprocessing, plotting, and statistical testing
- [ ] **Upstream proposal**: Open a design discussion with MuData maintainers on whether mappings should become a first-class `.mappings` attribute (similar to `.obsm`/`.varm`)

---

## 8. Open questions

1. **Should unmapped features be allowed?** Currently, mapping indices can be a subset of `var_names`. Should we also support features in `var_names` that have no mapping at all (orphans), and if so, how should queries handle them?

2. **Transitive queries**: If we have `precursors → peptides` and `peptides → proteins`, should `get_related_features("precursors", "proteins", ...)` automatically compose the two mappings? This is useful but introduces complexity (sparse matrix multiplication, potential combinatorial explosion).

3. **Weighted mappings**: Some use cases may benefit from weighted edges (e.g. posterior probability that a precursor belongs to a protein). Should the adjacency matrix support weights, or should this remain a DataFrame-only edge attribute?

4. **MuData core vs. extension**: Is this best implemented as a standalone package that wraps MuData, or should the mapping infrastructure be proposed for inclusion in MuData itself? The former is more pragmatic; the latter would enable tighter integration (e.g. mappings that auto-update when modalities are subsetted).

5. **Naming**: `HierarchicalMuData` implies a strict hierarchy (DAG), but the current implementation supports arbitrary N:M bipartite graphs. Should the naming reflect this generality (e.g. `MappedMuData`, `LinkedMuData`)?
