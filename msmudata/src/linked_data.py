# MuData wrapper class to enable linked proteomics data construction, filtering and export

from __future__ import annotations

from typing import Any

import anndata as ad
import mudata as md
import pandas as pd

from .adjancency_matrix_construct import get_unique_mappings, adjacency_matrix_from_mapping
from .adjacency_matrix_subset import extract_feature_bounds_from_mudata, feature_index_to_adjacency_index, adjacency_index_to_feature_index, slice_associated_features

class LinkedData:
    """Linked data object for MS-based proteomics. 
    
    To support the multi-layered feature structure of MS-proteomics data,
    this class provides a wrapper around MuData with explicit feature mappings between
    feature levels (precursors, proteins, genes). 

    """

    def __init__(
        self,
        mdata: md.MuData,
    ) -> None:
        pass
