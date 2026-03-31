import numpy as np
import mudata as md

def extract_feature_bounds_from_mudata(
    mdata: md.MuData,
) -> dict[str, tuple[int, int]]:
    bounds = {}
    current_index = 0
    for key, mod in mdata.mod.items():
        n_vars = mod.n_vars
        bounds[key] = (current_index, current_index + n_vars)
        current_index += n_vars
    return bounds

def feature_index_to_adjacency_index(
    query_feature_index: int,
    feature_level: str,
    feature_bounds: dict[str, tuple[int, int]],
) -> int:
    """Return actual index of query in adjacency matrix."""
    # Find the start of the requested feature level
    if feature_level not in feature_bounds:
        raise ValueError(f"Feature level '{feature_level}' not found in feature bounds")
    
    start, end = feature_bounds[feature_level]
    
    level_size = end - start
    if query_feature_index < 0 or query_feature_index >= level_size:
        raise ValueError(f"Query index {query_feature_index} out of range for level '{feature_level}' (size: {level_size})")
    
    return start + query_feature_index

def adjacency_index_to_feature_index(
    adjacency_index: int,
    feature_bounds: dict[str, tuple[int, int]],
) -> tuple[str, int]:
    """Return feature level and index of query in adjacency matrix."""
    for level, (start, end) in feature_bounds.items():
        if start <= adjacency_index < end:
            return level, adjacency_index - start
    raise ValueError(f"Adjacency index {adjacency_index} out of bounds for feature bounds.")

def slice_associated_features(
    query_feature_index: int,
    feature_level: str,
    feature_bounds: dict[str, tuple[int, int]],
    adjacency_matrix: np.ndarray,
) -> dict[str, list[int]]:
    """Convert a query feature index into a map of associated features across all levels.
    
    Works with both dense numpy arrays and sparse scipy matrices.
    
    Returns:
        Dictionary mapping feature level names to lists of feature indices that are
        connected to the query feature in the adjacency matrix.
    """
    from scipy import sparse
    
    # Convert query feature index to adjacency matrix index
    query_adjacency_index = feature_index_to_adjacency_index(
        query_feature_index=query_feature_index,
        feature_level=feature_level,
        feature_bounds=feature_bounds,
    )

    # Slice adjacency matrix to get associated features (where value is 1)
    # Handle both dense and sparse matrices
    if sparse.issparse(adjacency_matrix):
        # For sparse matrices, use getrow() to maintain 2D shape
        row = adjacency_matrix.getrow(query_adjacency_index)
        # Convert to dense 1D array for np.where
        row_dense = row.toarray().flatten()
        associated_adjacency_indices = np.where(row_dense == 1)[0]
    else:
        # For dense matrices, slice normally
        associated_adjacency_indices = np.where(adjacency_matrix[query_adjacency_index, :] == 1)[0]
    
    # Convert adjacency indices back to feature indices grouped by level
    associated_features = {level: [] for level in feature_bounds.keys()}
    
    for adj_idx in associated_adjacency_indices:
        level, feature_idx = adjacency_index_to_feature_index(
            adjacency_index=adj_idx,
            feature_bounds=feature_bounds,
        )
        associated_features[level].append(feature_idx)
    
    return associated_features