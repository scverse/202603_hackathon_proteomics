import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import alphapepttools as at
from scipy.sparse import csr_matrix
from itertools import combinations
from scipy import sparse

def supported_feature_level_parameters() -> dict[str, str]:
    """
    Return a dictionary of supported feature level parameters for the PSM to MuData conversion process.

    """
    return {
        "diann": {
            "precursor": {
                "intensity_column": "Precursor.Quantity",
                "feature_id_column": "Precursor.Id",
                "level": "psm"
            },
            "protein": {
                "intensity_column": "PG.MaxLFQ",
                "feature_id_column": "Protein.Group",
                "level": "protein"
            },
            "gene": {
                "intensity_column": "Genes.MaxLFQ",
                "feature_id_column": "Genes",
                "level": "gene"
            }
        }
    }
    

def get_unique_mappings(psm_path: str, feature_level_names: list[str]) -> pd.DataFrame:
    """
    Get unique mappings from PSM table for specified feature levels.
    
    Parameters:
    - psm_table: DataFrame containing PSM data with columns for each feature level.
    - feature_level_names: List of column names corresponding to feature levels (e.g., ["Precursor", "Protein"]).
    
    Returns:
    - DataFrame with unique mappings between the specified feature levels.
    """
    # Select only the relevant columns for mapping
    psm_df = pd.read_csv(psm_path, sep="\t")
    mapping_df = psm_df.loc[:, feature_level_names].drop_duplicates()
    
    return mapping_df

def sparse_matrix_mapping(unique_mapping_df: pd.DataFrame) -> csr_matrix:
    """
    Create a square sparse adjacency matrix for varp from a feature-level mapping.

    Parameters
    ----------
    mapping_df : pd.DataFrame
        DataFrame where the index contains source features (e.g., precursors)
        and each column contains target features (e.g., protein groups).
        Produced by get_unique_mappings().

    Returns
    -------
    csr_matrix
        Square adjacency matrix of shape (n_total, n_total) where
        n_total = n_source + n_target features.
    """
    all_values = pd.unique(unique_mapping_df.values.ravel())                                                                                   
    value_to_idx = {v: i for i, v in enumerate(all_values)}                                                                     
    n = len(all_values)                                                                                                         
                                                                                                                                
    rows, cols = [], []                                                                                                         
    for _, row in unique_mapping_df.iterrows():                                                                                                
        for v1, v2 in combinations(row.values, 2):                                                                                     
            i, j = value_to_idx[v1], value_to_idx[v2]                                                                           
            rows.extend([i, j])                                                                                                 
            cols.extend([j, i])                                                                                                 
                                                                                                                                
    adj = sparse.coo_matrix(                                                                                                    
        (np.ones(len(rows), dtype=np.float64), (rows, cols)),
        shape=(n, n),                                                                                                           
    ).tocsr()                                                                                                                     
    adj.data = np.ones_like(adj.data)  # collapse summed duplicates to 1 
    adj = pd.DataFrame.sparse.from_spmatrix(adj, index=all_values, columns=all_values)                                                                                                                      
    return adj

def create_mudata_diann(psm_path: str, feature_level_names: list[str]) -> md.MuData:
    """
    create  a MuData object from the PSM table, including unique mappings between feature levels.
   
   Parameters:
    - psm_path: Path to the Diann PSM table file (e.g., TSV or CSV).
    - feature_level_names: List of column names in the PSM table that represent the feature levels to map (e.g., ["Precursor.Id", "Protein.Group"]).
    Returns:
    - DataFrame with unique mappings between the specified feature levels.
    """
    # Get supported levels
    level_parameters = supported_feature_level_parameters()['diann']
    supported_levels = level_parameters.keys()

    # TODO: shift this over to AlphaBase standardization for reading
    # For now, get actual data column for requested levels from supported levels dic
    actual_levels = []
    for level in feature_level_names:
        if level not in supported_levels:
            raise ValueError(f"Unsupported feature level '{level}'. Supported levels for Diann are: {supported_levels}")
        actual_levels.append(level_parameters[level]['feature_id_column'])

    #Create mapping between feature levels
    mapping_df = get_unique_mappings(psm_path, actual_levels)
    matrix_varp = sparse_matrix_mapping(mapping_df)

    # iterate requested levels and load data
    layer_dict = {}
    for level in feature_level_names:
        params = level_parameters[level]
        adata = at.io.read_psm_table(psm_path, level=params['level'], search_engine="diann", 
                                    intensity_column=params['intensity_column'], feature_id_column=params['feature_id_column'],
                                    sample_id_column="Run"
                                    )
        layer_dict[f"{level}_level"] = adata

    mudata = md.MuData(
        layer_dict
    )

    mudata.varp["feature_mapping"] = matrix_varp.reindex(index=mudata.var_names, columns=mudata.var_names)

    return mudata