import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import alphapepttools as at
from scipy.sparse import csr_matrix
from itertools import combinations
from scipy import sparse

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
    #Create mapping between feature levels
    mapping_df = get_unique_mappings(psm_path, feature_level_names)
    matrix_varp = sparse_matrix_mapping(mapping_df)

    # load precursor and protein data from diann with alphapepttools functions
    prec_adata = at.io.read_psm_table(psm_path, level="psm", search_engine="diann", 
                                    intensity_column="Precursor.Quantity", feature_id_column="Precursor.Id",
                                    sample_id_column="Run", var_columns=feature_level_names
                                    )

    prot_adata = at.io.read_psm_table(psm_path, level="protein", search_engine="diann",
                                    intensity_column="PG.MaxLFQ", feature_id_column="Protein.Group",
                                    sample_id_column="Run"
                                    )
    
    #create mudata object
    mudata = md.MuData(
    # These are the raw data levels
        {
            "protein_level": prot_adata,
            #"peptide_level": ad.AnnData(...),
            "precursor_level": prec_adata,
        },
    )

    mudata.varp["feature_mapping"] = matrix_varp.reindex(index=mudata.var_names, columns=mudata.var_names)

    return mudata