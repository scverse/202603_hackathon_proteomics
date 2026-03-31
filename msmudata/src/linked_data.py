# MuData wrapper class to enable linked proteomics data construction, filtering and export

from __future__ import annotations

from typing import Any

import anndata as ad
import mudata as md
import pandas as pd
import logging

from .adjancency_matrix_construct import get_unique_mappings, adjacency_matrix_from_mapping
from .adjacency_matrix_subset import extract_feature_bounds_from_mudata, feature_index_to_adjacency_index, adjacency_index_to_feature_index, slice_associated_features

logger = logging.getLogger(__name__)   
class LinkedData:
    """Linked data object for MS-based proteomics. 
    
    To support the multi-layered feature structure of MS-proteomics data,
    this class provides a wrapper around MuData with explicit feature mappings between
    feature levels (precursors, proteins, genes). 

    """

    def __init__(
        self,
    ) -> None:
        pass

    def validate_mapping(
        self,
    ) -> None:
        """Validate that the mapping features are present in the MuData modalities."""
        logger.info("NOT IMPLEMENTED - validate that mapping and mudata match up")

    def adjacency_matrix_from_report(
        self,
        psm_report: pd.DataFrame,
        feature_level_names: list[str],
    ) -> None:
        """Construct the adjacency matrix from a PSM report"""

        # Construct the mapping
        self._mapping_df = get_unique_mappings(
            psm_table = psm_report,
            feature_level_names = feature_level_names,
        )

        # Construct the adjacency matrix
        self._adjacency_matrix = adjacency_matrix_from_mapping(self._mapping_df)

        # validate adjacency matrix against existing mudata
        #TODO: handle this better
        if hasattr(self, "_mdata"):
            self.validate_mapping()

    def add_mudata(
        self,
        mdata: md.MuData,
    ) -> None:
        """Add a MuData object to the LinkedData instance, and extract feature bounds."""
        self._mdata = mdata
        self._feature_bounds = extract_feature_bounds_from_mudata(mdata) 

        # validate adjacency matrix against new mudata
        #TODO: handle this better
        if hasattr(self, "_adjacency_matrix"):
            self.validate_mapping()



