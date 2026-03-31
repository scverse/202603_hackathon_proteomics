# MuData wrapper class to enable linked proteomics data construction, filtering and export

from __future__ import annotations

from typing import Optional, Union
import warnings

import anndata as ad
import mudata as md
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix
import logging

from .adjancency_matrix_construct import get_unique_mappings, adjacency_matrix_from_mapping
from .adjacency_matrix_subset import extract_feature_bounds_from_mudata, slice_associated_features

logger = logging.getLogger(__name__)


class LinkedData:
    """Linked data object for MS-based proteomics.

    To support the multi-layered feature structure of MS-proteomics data,
    this class provides a wrapper around MuData with explicit feature mappings between
    feature levels (precursors, proteins, genes).

    The LinkedData class can be initialized in three ways:
    1. With a MuData that already has adjacency matrix in varp
    2. With a MuData and a separate PSM report to construct the adjacency matrix
    3. Empty, then add MuData and/or adjacency matrix later

    Attributes
    ----------
    mdata : md.MuData
        The underlying MuData object
    varp : dict
        Mirror of the MuData varp, providing access to the adjacency matrix
    feature_bounds : dict
        Boundaries of each modality in the adjacency matrix

    Examples
    --------
    >>> # Initialize with MuData that has adjacency matrix
    >>> ld = LinkedData(mdata)

    >>> # Initialize with MuData and construct adjacency from report
    >>> ld = LinkedData(mdata, psm_report=report, feature_levels=['genes', 'proteins', 'precursors'])

    >>> # Initialize empty and add components
    >>> ld = LinkedData()
    >>> ld.add_mudata(mdata)
    >>> ld.construct_adjacency_from_report(report, ['genes', 'proteins', 'precursors'])
    """

    def __init__(
        self,
        mdata: Optional[md.MuData] = None,
        psm_report: Optional[pd.DataFrame] = None,
        feature_levels: Optional[list[str]] = None,
        adjacency_key: str = "feature_mapping",
    ) -> None:
        """Initialize LinkedData object.

        Parameters
        ----------
        mdata : Optional[md.MuData]
            MuData object to wrap. Can contain adjacency matrix in varp.
        psm_report : Optional[pd.DataFrame]
            PSM report to construct adjacency matrix from
        feature_levels : Optional[list[str]]
            Feature level names for adjacency matrix construction
        adjacency_key : str, default="feature_mapping"
            Key to use for adjacency matrix in varp
        """
        self._mdata = None
        self._adjacency_matrix = None
        self._feature_bounds = None
        self._adjacency_key = adjacency_key
        self._filter_dict = None

        # Initialize with MuData if provided
        if mdata is not None:
            self.add_mudata(mdata)

        # Construct adjacency from report if provided
        if psm_report is not None:
            if feature_levels is None:
                raise ValueError("feature_levels must be provided with psm_report")
            self.construct_adjacency_from_report(psm_report, feature_levels)

    @property
    def mdata(self) -> Optional[md.MuData]:
        """Get the underlying MuData object."""
        return self._mdata

    @property
    def varp(self) -> dict:
        """Access varp dictionary, mirroring MuData varp.

        Returns the varp from MuData if it exists, otherwise returns
        a dictionary with the adjacency matrix if constructed separately.
        """
        if self._mdata is not None and hasattr(self._mdata, 'varp'):
            return self._mdata.varp
        elif self._adjacency_matrix is not None:
            return {self._adjacency_key: self._adjacency_matrix}
        else:
            return {}

    @property
    def feature_bounds(self) -> Optional[dict]:
        """Get feature boundaries for the adjacency matrix."""
        if self._feature_bounds is None and self._mdata is not None:
            self._feature_bounds = extract_feature_bounds_from_mudata(self._mdata)
        return self._feature_bounds

    @property
    def adjacency_matrix(self) -> Optional[Union[np.ndarray, sparse.spmatrix]]:
        """Get the adjacency matrix from varp."""
        return self.varp.get(self._adjacency_key)

    def add_mudata(
        self,
        mdata: md.MuData,
        overwrite: bool = False,
    ) -> None:
        """Add or update the MuData object.

        Parameters
        ----------
        mdata : md.MuData
            MuData object to add
        overwrite : bool, default=False
            Whether to overwrite existing MuData
        """
        if self._mdata is not None and not overwrite:
            raise ValueError("MuData already exists. Set overwrite=True to replace.")

        self._mdata = mdata
        self._feature_bounds = extract_feature_bounds_from_mudata(mdata)

        # If MuData has adjacency matrix in varp, use it
        if self._adjacency_key in mdata.varp:
            logger.info(f"Using existing adjacency matrix from MuData varp['{self._adjacency_key}']")
            self._adjacency_matrix = mdata.varp[self._adjacency_key]

        # Validate if we have both MuData and adjacency matrix
        if self._adjacency_matrix is not None:
            self._validate_adjacency_against_mudata()

    def construct_adjacency_from_report(
        self,
        psm_report: pd.DataFrame,
        feature_level_names: list[str],
        sparse_matrix: bool = True,
    ) -> None:
        """Construct adjacency matrix from a PSM report.

        Parameters
        ----------
        psm_report : pd.DataFrame
            PSM report with feature mappings
        feature_level_names : list[str]
            Column names in psm_report for each feature level
        sparse_matrix : bool, default=True
            Whether to store as sparse matrix
        """
        # Construct the mapping
        mapping_df = get_unique_mappings(
            psm_table=psm_report,
            feature_level_names=feature_level_names,
        )

        # Construct the adjacency matrix
        adjacency = adjacency_matrix_from_mapping(mapping_df)

        # Convert to sparse if requested
        if sparse_matrix and not sparse.issparse(adjacency):
            adjacency = csr_matrix(adjacency)

        self._adjacency_matrix = adjacency

        # If we have MuData, update its varp
        if self._mdata is not None:
            # Reindex to match MuData var names
            if isinstance(adjacency, pd.DataFrame):
                adjacency = adjacency.reindex(
                    index=self._mdata.var_names,
                    columns=self._mdata.var_names
                )
            self._mdata.varp[self._adjacency_key] = adjacency
            self._validate_adjacency_against_mudata()

        logger.info(f"Constructed adjacency matrix with shape {adjacency.shape}")

    def _validate_adjacency_against_mudata(self) -> None:
        """Validate that adjacency matrix dimensions match MuData."""
        if self._mdata is None or self._adjacency_matrix is None:
            return

        n_vars = self._mdata.n_vars
        adj_shape = self._adjacency_matrix.shape

        if adj_shape[0] != n_vars or adj_shape[1] != n_vars:
            warnings.warn(
                f"Adjacency matrix shape {adj_shape} doesn't match "
                f"MuData n_vars ({n_vars}). Consider reconstructing."
            )

    def _query_feature_to_index(
        self,
        query_feature: Union[int, str],
        feature_level: str,
    ) -> int:
        """Convert query feature to the corresponding feature index from the AnnData.

        Parameters
        ----------
        query_feature : Union[int, str]
            Either the feature name (str) or index (int)
        feature_level : str
            The modality name

        Returns
        -------
        int
            The numerical index of the feature
        """
        if isinstance(query_feature, int):
            # Already an index, just return it
            return query_feature
        elif isinstance(query_feature, str):
            # Convert feature name to index
            feature_adata = self._mdata.mod[feature_level]
            return feature_adata.var_names.get_loc(query_feature)
        else:
            raise TypeError(f"query_feature must be int or str, got {type(query_feature)}")

    def get_associated_features(
        self,
        query_feature: Union[int, str],
        feature_level: str,
    ) -> dict[str, list[int]]:
        """Get features associated with a query feature across all levels and store as filter.

        Parameters
        ----------
        query_feature : Union[int, str]
            Either the feature name (e.g., 'gene0', 'protein1') or its index (0, 1, 2...)
        feature_level : str
            Name of the modality containing the query feature (e.g., 'genes', 'proteins')

        Returns
        -------
        dict[str, list[int]]
            Dictionary mapping modality names to lists of associated feature indices

        Examples
        --------
        >>> # Query by feature name
        >>> ld.get_associated_features('gene0', 'genes')
        {'genes': [0], 'proteins': [0, 1], 'precursors': [0, 1, 2]}

        >>> # Query by feature index
        >>> ld.get_associated_features(0, 'genes')
        {'genes': [0], 'proteins': [0, 1], 'precursors': [0, 1, 2]}
        """
        if self.adjacency_matrix is None:
            raise ValueError("No adjacency matrix available")

        if self.feature_bounds is None:
            raise ValueError("No feature bounds available")

        # Convert feature name to index if necessary
        query_feature_index = self._query_feature_to_index(query_feature, feature_level)

        # Get associated features and store as filter
        self._filter_dict = slice_associated_features(
            query_feature_index=query_feature_index,
            feature_level=feature_level,
            feature_bounds=self.feature_bounds,
            adjacency_matrix=self.adjacency_matrix
        )

        return self._filter_dict

    def to_anndata(
        self,
        feature_level: str,
    ) -> ad.AnnData:
        """Extract filtered AnnData for a specific modality.

        Uses the filter dictionary set by get_associated_features to return
        a filtered copy of the AnnData for the specified modality.

        Parameters
        ----------
        feature_level : str
            Name of the modality to extract (e.g., 'genes', 'proteins', 'precursors')

        Returns
        -------
        ad.AnnData
            Filtered copy of the AnnData containing only the features in the filter

        Examples
        --------
        >>> # First set the filter
        >>> ld.get_associated_features('gene0', 'genes')
        >>> # Then extract filtered data
        >>> filtered_proteins = ld.to_anndata('proteins')
        >>> filtered_precursors = ld.to_anndata('precursors')
        """
        if self._filter_dict is None:
            raise ValueError("No filter set. Call get_associated_features first.")

        if self._mdata is None:
            raise ValueError("No MuData available")

        if feature_level not in self._mdata.mod:
            raise ValueError(f"Modality '{feature_level}' not found. Available: {list(self._mdata.mod.keys())}")

        # Get the indices for this feature level from the filter
        if feature_level not in self._filter_dict:
            raise ValueError(f"Feature level '{feature_level}' not in filter dict")

        feature_indices = self._filter_dict[feature_level]

        if not feature_indices:
            logger.warning(f"No features to filter for {feature_level}")
            # Return empty AnnData with same obs
            adata = self._mdata.mod[feature_level]
            return ad.AnnData(
                X=adata.X[:, :0],
                obs=adata.obs.copy(),
                var=pd.DataFrame(index=pd.Index([], name=adata.var.index.name))
            )

        # Get the AnnData and filter by features
        adata = self._mdata.mod[feature_level]
        filtered_adata = adata[:, feature_indices].copy()

        return filtered_adata

    def __repr__(self) -> str:
        """String representation of LinkedData object."""
        parts = ["LinkedData object"]

        if self._mdata is not None:
            parts.append(f"  MuData: {self._mdata.n_obs} obs × {self._mdata.n_vars} vars")
            parts.append(f"  Modalities: {', '.join(self._mdata.mod.keys())}")
        else:
            parts.append("  MuData: None")

        if self.adjacency_matrix is not None:
            parts.append(f"  Adjacency matrix: {self.adjacency_matrix.shape}")
            if sparse.issparse(self.adjacency_matrix):
                parts.append(f"    Sparse with {self.adjacency_matrix.nnz} connections")
        else:
            parts.append("  Adjacency matrix: None")

        return "\n".join(parts)



