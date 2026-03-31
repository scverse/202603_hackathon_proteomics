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

| Project                                                                                                                | Approach                                                        | Limitation                     |
| ---------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------- | ------------------------------ |
| [QFeatures](https://rformassspectrometry.github.io/QFeatures/articles/QFeatures.html) (R/Bioconductor)                 | `AssayLinks` connecting `SummarizedExperiment` layers via a DAG | R-only, not scverse-compatible |
| [alphaquant](https://github.com/MannLabs/alphaquant)                                                                   | Internal mapping tables (DataFrames)                            | Not formalized, not portable   |
| [alphabase scverse tutorial](https://alphabase.readthedocs.io/en/latest/tutorials/tutorial_scverse_compatibility.html) | Reads PSM reports into AnnData/MuData                           | No cross-modality mapping      |
| mudata#111 proposals by @mffrank                                                                                       | DAG mapping proof-of-concept                                    | Not integrated                 |
| mudata#111 proposals by @ilan-gold                                                                                     | Pandas extension index arrays                                   | Design-stage                   |

### 1.3 Scope

This RFC proposes a **minimal, composable extension** to MuData for representing directed, queryable N:M feature mappings between modalities. The design is motivated by MS proteomics but generalizes to any multi-modal setting with hierarchical or relational feature structure (e.g. ATAC peaks → genes, spatial transcriptomics spots → cells).

---

## Core concepts

### Store a feature mapping in mudata.varp

Use the `mudata.varp` attribute to formalize the inherent tree structure of MS data as DAG, as proposed by [Ilia Kats](https://github.com/scverse/mudata/issues/111#issuecomment-3772341577)

Given a mdata object with a total number of features $p=p_0 + p_1 + ... + p_n$ from feature **levels** `0`, `1`, ..., `N` The `varp` attribute has the shape `p x p`

A feature mapping in the `varp` attribute represents the _binary adjacency matrix_ of shape `p x p` of a a tree-like directed acylic graph that stores the mapping between features of the different levels. The indices of the adjacency matrix are aligned with `mdata.var_names`. A `1` in the the adjacency matrix at position `i, j` represents a connection between a feature `mdata.var_names[i]` and `mdata.var_names[j]`.

**Advantages** This approach naturally extends to any number of feature levels and `n:m` mappings. The support of sparse matrices in `mdata.varp` makes the storage of this feature mapping very efficient.

**Alternatives**:
Pandas index extension - The provided approach could be further extended with a custom `pandas.Index` extension that faciliates the querying.

**New attribute** Add a `FeatureMapping` attribute, similar to the `QFeatures` implementation in R: This approach does not leverage the existing infrastructure of the `mudata` object (e.g. with respect to serialization) and is harder to extend to an arbitrary number of feature levels.

### Querying

### Validation

In future extensions, the obtained graph could be validated. E.g. that mappings between feature-levels follow further constraints

- `N:1`, `1:N` mapping
- unidirectonality across feature levels
- existence of unmapped features

### Storage

The object can be readily serialized with the existing infrastructure

## Ecosystem

### Data ingestion

**Data ingestion**: Reader that constructs mudata from an alphadia / DIA-NN / QuantMS PSM report in a single call

_Work in progress_

### Plotting

For a parent node, plot the values of all child nodes.

_Work in progress_

## Future developments

- [ ] **Transitive queries**: If we have `precursors → peptides` and `peptides → proteins`, should `get_related_features("precursors", "proteins", ...)` automatically compose the two mappings? This is useful but introduces complexity (sparse matrix multiplication, potential combinatorial explosion).

- [ ] **Weighted mappings**: Some use cases may benefit from weighted edges (e.g. posterior probability that a precursor belongs to a protein). Should the adjacency matrix support weights, or should this remain a DataFrame-only edge attribute?

- [ ] **MuData core vs. extension**: Is this best implemented as a standalone package that wraps MuData, or should the mapping infrastructure be proposed for inclusion in MuData itself? The former is more pragmatic; the latter would enable tighter integration (e.g. mappings that auto-update when modalities are subsetted).

- [ ] **Naming**
