import pandas as pd 

import pandas as pd                                                                                                             
import numpy as np                                                                                                            
from scipy import sparse                                                                                                        
from itertools import combinations

def get_unique_mappings(psm_table: pd.DataFrame, feature_level_names: list[str]) -> pd.DataFrame:
    """Get all unique mappings"""
    return psm_table.loc[:, feature_level_names].drop_duplicates(subset=feature_level_names)
    

def adjacency_matrix_from_mapping(df: pd.DataFrame) -> pd.DataFrame:                                                            
    """Convert a mapping DataFrame to a sparse undirected adjacency matrix.                                                     
                                                                                                                                
    Each row creates edges between all pairs of values across columns,                                                               
    representing hierarchical relationships.                                                   
    """                                                                                                                         
    all_values = pd.unique(df.values.ravel())                                                                                   
    value_to_idx = {v: i for i, v in enumerate(all_values)}                                                                     
    n = len(all_values)                                                                                                         
                                                                                                                                
    rows, cols = [], []                                                                                                         
    for _, row in df.iterrows():                                                                                                
        for v1, v2 in combinations(row.values, 2):                                                                                     
            i, j = value_to_idx[v1], value_to_idx[v2]                                                                           
            rows.extend([i, j])                                                                                                 
            cols.extend([j, i])                                                                                                 
                                                                                                                                
    adj = sparse.coo_matrix(                                                                                                    
        (np.ones(len(rows), dtype=np.float64), (rows, cols)),
        shape=(n, n),                                                                                                           
    ).tocsr()                                                                                                                     
    adj.data = np.ones_like(adj.data)  # collapse summed duplicates to 1                                                                                                                        
    return pd.DataFrame.sparse.from_spmatrix(adj, index=all_values, columns=all_values)