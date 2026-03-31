import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import shortest_path


def _to_csr(mat):
    if sp.issparse(mat):
        return mat.tocsr()
    return sp.csr_matrix(mat)


def _get_modality_feature_indices(mdata, modality_name):
    """
    Return global feature indices in mdata.var_names for one modality.
    """
    global_index = {name: i for i, name in enumerate(mdata.var_names)}
    feats = list(mdata.mod[modality_name].var_names)
    return np.array([global_index[f] for f in feats], dtype=int)


def _extract_submapping(mdata, source_mod, target_mod, mapping_key="feature_mapping"):
    """
    Get direct source->target mapping block from mdata.varp[mapping_key].
    Returns sparse matrix of shape (n_source_features, n_target_features).
    """
    fm = _to_csr(mdata.varp[mapping_key])
    src_idx = _get_modality_feature_indices(mdata, source_mod)
    tgt_idx = _get_modality_feature_indices(mdata, target_mod)
    sub = fm[src_idx][:, tgt_idx]
    return sub


def _compose_path_mapping(mdata, path, mapping_key="feature_mapping", binarize=True):
    """
    Compose mapping along a modality path, e.g.
    ['precursors', 'proteins', 'genes'].

    Returns matrix of shape:
      (#features in first modality, #features in last modality)

    If binarize=True, any nonzero composed value becomes 1.
    """
    if len(path) < 2:
        raise ValueError("path must contain at least 2 modalities")

    cur = _extract_submapping(mdata, path[0], path[1], mapping_key=mapping_key)

    for i in range(1, len(path) - 1):
        nxt = _extract_submapping(mdata, path[i], path[i + 1], mapping_key=mapping_key)
        cur = cur @ nxt

    if binarize:
        cur = cur.copy()
        cur.data[:] = 1.0
        cur.eliminate_zeros()

    return cur.tocsr()


def _find_modality_path_from_mapping(mdata, source_mod, target_mod, mapping_key="feature_mapping"):
    """
    Infer a modality-level path using existence of nonzero mapping blocks.
    Example:
      precursors -> proteins -> genes

    Returns a list like ['precursors', 'proteins', 'genes'].

    This uses the modality graph induced by nonzero block mappings.
    """
    mods = list(mdata.mod.keys())
    n = len(mods)

    # Build adjacency between modalities based on nonzero direct mappings
    adj = np.full((n, n), np.inf, dtype=float)
    np.fill_diagonal(adj, 0.0)

    for i, src in enumerate(mods):
        for j, tgt in enumerate(mods):
            if i == j:
                continue
            block = _extract_submapping(mdata, src, tgt, mapping_key=mapping_key)
            if block.nnz > 0:
                adj[i, j] = 1.0

    dist, predecessors = shortest_path(
        adj, directed=True, return_predecessors=True
    )

    src_i = mods.index(source_mod)
    tgt_i = mods.index(target_mod)

    if np.isinf(dist[src_i, tgt_i]):
        raise ValueError(f"No modality path found from '{source_mod}' to '{target_mod}'")

    # Reconstruct path
    path_idx = []
    cur = tgt_i
    while cur != src_i:
        path_idx.append(cur)
        cur = predecessors[src_i, cur]
        if cur < 0:
            raise ValueError(f"Failed to reconstruct path from '{source_mod}' to '{target_mod}'")
    path_idx.append(src_i)
    path_idx.reverse()

    return [mods[i] for i in path_idx]


def _aggregate_matrix(X, mapping, agg="sum", split_shared=False):
    """
    Aggregate X using source->target mapping.

    Parameters
    ----------
    X : array-like or sparse, shape (n_obs, n_source_features)
    mapping : sparse matrix, shape (n_source_features, n_target_features)
        source feature -> target feature
    agg : {'sum', 'mean', 'median', 'max', 'min'}
    split_shared : bool
        If True, divide each source feature contribution equally across all mapped targets.

    Returns
    -------
    out : np.ndarray, shape (n_obs, n_target_features)
    """
    if sp.issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X)

    mapping = _to_csr(mapping).copy()

    if split_shared:
        row_sums = np.asarray(mapping.sum(axis=1)).ravel()
        nz = row_sums > 0
        if nz.any():
            scale = np.ones_like(row_sums, dtype=float)
            scale[nz] = 1.0 / row_sums[nz]
            mapping = sp.diags(scale) @ mapping

    n_obs, n_src = X.shape
    if mapping.shape[0] != n_src:
        raise ValueError(
            f"Shape mismatch: X has {n_src} source features but mapping has {mapping.shape[0]}"
        )

    n_tgt = mapping.shape[1]

    if agg == "sum":
        return X @ mapping.toarray()

    out = np.full((n_obs, n_tgt), np.nan, dtype=float)

    mapping_csc = mapping.tocsc()

    for j in range(n_tgt):
        src_ids = mapping_csc.indices[mapping_csc.indptr[j]:mapping_csc.indptr[j + 1]]
        if len(src_ids) == 0:
            continue

        vals = X[:, src_ids]

        if split_shared and agg != "sum":
            # For non-sum aggregations, weighted split is not very natural.
            # We apply per-source scaling before aggregation if requested.
            weights = np.asarray(mapping[src_ids, j].todense()).ravel()
            vals = vals * weights[None, :]

        if agg == "mean":
            out[:, j] = np.mean(vals, axis=1)
        elif agg == "median":
            out[:, j] = np.median(vals, axis=1)
        elif agg == "max":
            out[:, j] = np.max(vals, axis=1)
        elif agg == "min":
            out[:, j] = np.min(vals, axis=1)
        else:
            raise ValueError(f"Unsupported agg='{agg}'")

    return out


def aggregate_between_modalities(
    mdata,
    source_mod,
    target_mod,
    agg="sum",
    mapping_key="feature_mapping",
    path=None,
    infer_path=False,
    binarize_composed=True,
    split_shared=False,
    source_layer=None,
    target_layer=None,
    write_to_X=False,
    return_mapping=False,
):
    """
    Aggregate values from one modality to another using mdata.varp[mapping_key].

    Parameters
    ----------
    mdata : MuData
    source_mod : str
    target_mod : str
    agg : str
        'sum', 'mean', 'median', 'max', 'min'
    mapping_key : str
        Key in mdata.varp
    path : list[str] | None
        Explicit modality path, e.g. ['precursors', 'proteins', 'genes'].
        If None and source_mod -> target_mod is direct, direct mapping is used.
    infer_path : bool
        If True and path is None, infer a modality path automatically.
    binarize_composed : bool
        If composing multiple mappings, binarize nonzero entries to 1.
    split_shared : bool
        Split one source feature's contribution across multiple targets.
        Most natural for 'sum'.
    source_layer : str | None
        If set, read from mdata.mod[source_mod].layers[source_layer] instead of X.
    target_layer : str | None
        If set, write result into mdata.mod[target_mod].layers[target_layer].
    write_to_X : bool
        If True, also write result into mdata.mod[target_mod].X.
    return_mapping : bool
        If True, also return the mapping used.

    Returns
    -------
    out : np.ndarray
    mapping : sparse matrix (optional)
    """
    if source_layer is None:
        X = mdata.mod[source_mod].X
    else:
        X = mdata.mod[source_mod].layers[source_layer]

    # Determine mapping
    if path is not None:
        if path[0] != source_mod or path[-1] != target_mod:
            raise ValueError("Explicit path must start with source_mod and end with target_mod")
        mapping = _compose_path_mapping(
            mdata, path, mapping_key=mapping_key, binarize=binarize_composed
        )
    else:
        direct = _extract_submapping(mdata, source_mod, target_mod, mapping_key=mapping_key)
        if direct.nnz > 0:
            mapping = direct
        elif infer_path:
            found_path = _find_modality_path_from_mapping(
                mdata, source_mod, target_mod, mapping_key=mapping_key
            )
            mapping = _compose_path_mapping(
                mdata, found_path, mapping_key=mapping_key, binarize=binarize_composed
            )
        else:
            raise ValueError(
                f"No direct mapping from '{source_mod}' to '{target_mod}'. "
                f"Provide path=... or set infer_path=True."
            )

    out = _aggregate_matrix(X, mapping, agg=agg, split_shared=split_shared)

    target_adata = mdata.mod[target_mod]
    if out.shape != target_adata.shape:
        raise ValueError(
            f"Aggregated output shape {out.shape} does not match target modality shape {target_adata.shape}"
        )

    if target_layer is not None:
        target_adata.layers[target_layer] = out

    if write_to_X:
        target_adata.X = out

    if return_mapping:
        return out, mapping
    return out