"""FeatureMapping: stores and queries N:M relationships between two feature sets.

This is the core data structure for representing hierarchical feature relationships
in MS proteomics, e.g. precursor <-> protein mappings.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd
from scipy import sparse


class FeatureMapping:
    """Bidirectional N:M mapping between two indexed feature sets.

    Parameters
    ----------
    source_index
        Unique identifiers for source features (e.g. precursor IDs).
    target_index
        Unique identifiers for target features (e.g. protein IDs).
    edges
        DataFrame with columns ``["source", "target"]`` containing pairs of
        feature identifiers that are related. Optionally, additional columns
        (e.g. ``"weight"``) are preserved as edge attributes.

    Examples
    --------
    >>> mapping = FeatureMapping(
    ...     source_index=["prec_1", "prec_2", "prec_3"],
    ...     target_index=["prot_A", "prot_B"],
    ...     edges=pd.DataFrame({
    ...         "source": ["prec_1", "prec_1", "prec_2", "prec_3"],
    ...         "target": ["prot_A", "prot_B", "prot_A", "prot_B"],
    ...     }),
    ... )
    >>> mapping.targets_of("prec_1")
    ['prot_A', 'prot_B']
    >>> mapping.sources_of("prot_B")
    ['prec_1', 'prec_3']
    """

    def __init__(
        self,
        source_index: pd.Index | Sequence[str],
        target_index: pd.Index | Sequence[str],
        edges: pd.DataFrame,
    ) -> None:
        self._source_index = pd.Index(source_index)
        self._target_index = pd.Index(target_index)

        if not {"source", "target"}.issubset(edges.columns):
            raise ValueError("edges must contain 'source' and 'target' columns")

        unknown_sources = set(edges["source"]) - set(self._source_index)
        if unknown_sources:
            raise ValueError(
                f"{len(unknown_sources)} source(s) in edges not found in source_index, "
                f"e.g. {list(unknown_sources)[:5]}"
            )
        unknown_targets = set(edges["target"]) - set(self._target_index)
        if unknown_targets:
            raise ValueError(
                f"{len(unknown_targets)} target(s) in edges not found in target_index, "
                f"e.g. {list(unknown_targets)[:5]}"
            )

        self._edges = edges.reset_index(drop=True)
        self._adjacency: sparse.csr_matrix | None = None

    # ---- properties ----

    @property
    def source_index(self) -> pd.Index:
        return self._source_index

    @property
    def target_index(self) -> pd.Index:
        return self._target_index

    @property
    def edges(self) -> pd.DataFrame:
        return self._edges

    @property
    def n_sources(self) -> int:
        return len(self._source_index)

    @property
    def n_targets(self) -> int:
        return len(self._target_index)

    @property
    def n_edges(self) -> int:
        return len(self._edges)

    # ---- sparse adjacency (lazy) ----

    @property
    def adjacency(self) -> sparse.csr_matrix:
        """Sparse boolean adjacency matrix of shape (n_sources, n_targets)."""
        if self._adjacency is None:
            source_codes = self._source_index.get_indexer(self._edges["source"])
            target_codes = self._target_index.get_indexer(self._edges["target"])
            data = np.ones(len(source_codes), dtype=bool)
            self._adjacency = sparse.csr_matrix(
                (data, (source_codes, target_codes)),
                shape=(self.n_sources, self.n_targets),
            )
        return self._adjacency

    def _invalidate_cache(self) -> None:
        self._adjacency = None

    # ---- query API ----

    def targets_of(self, source: str) -> list[str]:
        """Return all target features linked to a given source feature."""
        idx = self._source_index.get_loc(source)
        cols = self.adjacency[idx].nonzero()[1]
        return self._target_index[cols].tolist()

    def sources_of(self, target: str) -> list[str]:
        """Return all source features linked to a given target feature."""
        idx = self._target_index.get_loc(target)
        rows = self.adjacency[:, idx].nonzero()[0]
        return self._source_index[rows].tolist()

    def targets_of_many(self, sources: Sequence[str]) -> list[str]:
        """Return the union of target features for multiple source features."""
        idxs = self._source_index.get_indexer(sources)
        if (idxs == -1).any():
            missing = [s for s, i in zip(sources, idxs) if i == -1]
            raise KeyError(f"Sources not found: {missing}")
        cols = self.adjacency[idxs].nonzero()[1]
        return self._target_index[np.unique(cols)].tolist()

    def sources_of_many(self, targets: Sequence[str]) -> list[str]:
        """Return the union of source features for multiple target features."""
        idxs = self._target_index.get_indexer(targets)
        if (idxs == -1).any():
            missing = [t for t, i in zip(targets, idxs) if i == -1]
            raise KeyError(f"Targets not found: {missing}")
        rows = self.adjacency[:, idxs].nonzero()[0]
        return self._source_index[np.unique(rows)].tolist()

    def subgraph(
        self,
        sources: Sequence[str] | None = None,
        targets: Sequence[str] | None = None,
    ) -> FeatureMapping:
        """Return a new FeatureMapping restricted to the given sources/targets."""
        mask = pd.Series(True, index=self._edges.index)
        new_sources = self._source_index
        new_targets = self._target_index

        if sources is not None:
            new_sources = pd.Index(sources)
            mask &= self._edges["source"].isin(new_sources)
        if targets is not None:
            new_targets = pd.Index(targets)
            mask &= self._edges["target"].isin(new_targets)

        return FeatureMapping(
            source_index=new_sources,
            target_index=new_targets,
            edges=self._edges.loc[mask],
        )

    # ---- construction helpers ----

    @classmethod
    def from_long_dataframe(
        cls,
        df: pd.DataFrame,
        source_col: str,
        target_col: str,
        *,
        extra_cols: Sequence[str] | None = None,
    ) -> FeatureMapping:
        """Build a FeatureMapping from a long-format DataFrame.

        Useful for constructing mappings directly from PSM report columns.
        """
        edges_cols = ["source", "target"]
        edges = pd.DataFrame(
            {"source": df[source_col].values, "target": df[target_col].values}
        )
        if extra_cols:
            for col in extra_cols:
                edges[col] = df[col].values
                edges_cols.append(col)

        edges = edges.drop_duplicates(subset=["source", "target"]).reset_index(drop=True)

        return cls(
            source_index=pd.Index(edges["source"].unique()),
            target_index=pd.Index(edges["target"].unique()),
            edges=edges,
        )

    # ---- display ----

    def __repr__(self) -> str:
        return (
            f"FeatureMapping(n_sources={self.n_sources}, n_targets={self.n_targets}, "
            f"n_edges={self.n_edges})"
        )
