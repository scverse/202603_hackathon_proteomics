"""HierarchicalMuData: a MuData wrapper with queryable cross-modality feature mappings.

Extends MuData with explicit, directed feature-level relationships between
modalities — designed for the hierarchical quantification structure in
MS-based proteomics (precursors → peptides → proteins).
"""

from __future__ import annotations

from typing import Any

import anndata as ad
import mudata as md
import pandas as pd

from .feature_mapping import FeatureMapping


class HierarchicalMuData:
    """MuData with explicit hierarchical feature mappings between modalities.

    Parameters
    ----------
    mdata
        A regular MuData object containing the modalities.
    mappings
        Dictionary keyed by ``(source_mod, target_mod)`` tuples, with
        :class:`FeatureMapping` values encoding the feature relationships.

    Examples
    --------
    >>> hmdata = HierarchicalMuData(
    ...     mdata=mdata,
    ...     mappings={("precursors", "proteins"): precursor_protein_mapping},
    ... )
    >>> hmdata.get_targets("precursors", "proteins", feature="PEPTIDEK_2")
    ['P04637', 'Q9Y6K9']
    >>> precursor_adata = hmdata.slice_data("precursors", by="proteins", feature="P04637")
    """

    def __init__(
        self,
        mdata: md.MuData,
        mappings: dict[tuple[str, str], FeatureMapping] | None = None,
    ) -> None:
        self._mdata = mdata
        self._mappings: dict[tuple[str, str], FeatureMapping] = {}

        if mappings:
            for key, mapping in mappings.items():
                self.add_mapping(key[0], key[1], mapping)

    # ---- MuData delegation ----

    @property
    def mdata(self) -> md.MuData:
        """The underlying MuData object."""
        return self._mdata

    @property
    def mod(self) -> dict[str, ad.AnnData]:
        return self._mdata.mod

    @property
    def obs(self) -> pd.DataFrame:
        return self._mdata.obs

    @property
    def var(self) -> pd.DataFrame:
        return self._mdata.var

    @property
    def n_obs(self) -> int:
        return self._mdata.n_obs

    def __getitem__(self, key: str) -> ad.AnnData:
        return self._mdata[key]

    # ---- mapping management ----

    @property
    def mappings(self) -> dict[tuple[str, str], FeatureMapping]:
        """All registered feature mappings."""
        return dict(self._mappings)

    def add_mapping(
        self,
        source_mod: str,
        target_mod: str,
        mapping: FeatureMapping,
    ) -> None:
        """Register a feature mapping between two modalities.

        Validates that the mapping indices are consistent with the
        corresponding modality's ``var_names``.
        """
        if source_mod not in self._mdata.mod:
            raise ValueError(f"Source modality '{source_mod}' not in MuData")
        if target_mod not in self._mdata.mod:
            raise ValueError(f"Target modality '{target_mod}' not in MuData")

        source_var = set(self._mdata.mod[source_mod].var_names)
        mapped_sources = set(mapping.source_index)
        if not mapped_sources.issubset(source_var):
            diff = mapped_sources - source_var
            raise ValueError(
                f"{len(diff)} mapping source(s) not in '{source_mod}'.var_names, "
                f"e.g. {list(diff)[:5]}"
            )

        target_var = set(self._mdata.mod[target_mod].var_names)
        mapped_targets = set(mapping.target_index)
        if not mapped_targets.issubset(target_var):
            diff = mapped_targets - target_var
            raise ValueError(
                f"{len(diff)} mapping target(s) not in '{target_mod}'.var_names, "
                f"e.g. {list(diff)[:5]}"
            )

        self._mappings[(source_mod, target_mod)] = mapping

    def get_mapping(self, source_mod: str, target_mod: str) -> FeatureMapping:
        """Retrieve a mapping, automatically reversing direction if needed."""
        if (source_mod, target_mod) in self._mappings:
            return self._mappings[(source_mod, target_mod)]

        if (target_mod, source_mod) in self._mappings:
            return self._mappings[(target_mod, source_mod)]

        raise KeyError(
            f"No mapping registered between '{source_mod}' and '{target_mod}'"
        )

    # ---- cross-modality queries ----

    def _resolve_direction(
        self, source_mod: str, target_mod: str
    ) -> tuple[FeatureMapping, bool]:
        """Return (mapping, is_reversed) for a query direction."""
        if (source_mod, target_mod) in self._mappings:
            return self._mappings[(source_mod, target_mod)], False
        if (target_mod, source_mod) in self._mappings:
            return self._mappings[(target_mod, source_mod)], True
        raise KeyError(
            f"No mapping registered between '{source_mod}' and '{target_mod}'"
        )

    def get_targets(
        self, source_mod: str, target_mod: str, feature: str
    ) -> list[str]:
        """Get all features in *target_mod* that map to *feature* in *source_mod*."""
        mapping, reversed_ = self._resolve_direction(source_mod, target_mod)
        if reversed_:
            return mapping.sources_of(feature)
        return mapping.targets_of(feature)

    def get_sources(
        self, target_mod: str, source_mod: str, feature: str
    ) -> list[str]:
        """Get all features in *source_mod* that map to *feature* in *target_mod*."""
        return self.get_targets(target_mod, source_mod, feature)

    def get_related_features(
        self,
        query_mod: str,
        return_mod: str,
        features: str | list[str],
    ) -> list[str]:
        """Flexible query: given feature(s) in *query_mod*, return related features in *return_mod*."""
        if isinstance(features, str):
            return self.get_targets(query_mod, return_mod, features)

        mapping, reversed_ = self._resolve_direction(query_mod, return_mod)
        if reversed_:
            return mapping.sources_of_many(features)
        return mapping.targets_of_many(features)

    def slice_data(
        self,
        data_mod: str,
        *,
        by: str,
        features: str | list[str],
    ) -> ad.AnnData:
        """Slice an AnnData modality to features related to feature(s) in another modality.

        Parameters
        ----------
        data_mod
            The modality whose data to return.
        by
            The modality that *features* belong to.
        features
            One or more feature identifiers in the *by* modality.

        Returns
        -------
        A view of the AnnData from *data_mod*, filtered to only the features
        that map to the given *features* in the *by* modality.

        Examples
        --------
        >>> # Get precursor-level data for protein "P04637"
        >>> adata_slice = hmdata.slice_data("precursors", by="proteins", features="P04637")
        """
        related = self.get_related_features(by, data_mod, features)
        adata = self._mdata.mod[data_mod]
        mask = adata.var_names.isin(related)
        return adata[:, mask]

    def mapping_summary(self) -> pd.DataFrame:
        """Summary table of all registered mappings."""
        rows = []
        for (src, tgt), m in self._mappings.items():
            rows.append(
                {
                    "source_modality": src,
                    "target_modality": tgt,
                    "n_source_features": m.n_sources,
                    "n_target_features": m.n_targets,
                    "n_edges": m.n_edges,
                    "mean_targets_per_source": m.adjacency.sum(axis=1).mean(),
                    "mean_sources_per_target": m.adjacency.sum(axis=0).mean(),
                }
            )
        return pd.DataFrame(rows)

    # ---- serialization via MuData.uns ----

    def store_mappings_in_uns(self) -> None:
        """Persist all mappings into ``mdata.uns`` so they survive disk I/O."""
        for (src, tgt), mapping in self._mappings.items():
            key = f"mapping_{src}__{tgt}"
            self._mdata.uns[key] = {
                "source_mod": src,
                "target_mod": tgt,
                "edges": mapping.edges.to_dict(orient="list"),
                "source_index": mapping.source_index.tolist(),
                "target_index": mapping.target_index.tolist(),
            }

    @classmethod
    def from_mudata(cls, mdata: md.MuData) -> HierarchicalMuData:
        """Reconstruct a HierarchicalMuData from a MuData that has mappings stored in ``.uns``."""
        mappings = {}
        for key, val in mdata.uns.items():
            if not key.startswith("mapping_"):
                continue
            edges = pd.DataFrame(val["edges"])
            mapping = FeatureMapping(
                source_index=pd.Index(val["source_index"]),
                target_index=pd.Index(val["target_index"]),
                edges=edges,
            )
            mappings[(val["source_mod"], val["target_mod"])] = mapping
        return cls(mdata=mdata, mappings=mappings)

    # ---- display ----

    def __repr__(self) -> str:
        mod_str = ", ".join(
            f"{name}: {a.n_obs} x {a.n_vars}" for name, a in self._mdata.mod.items()
        )
        mapping_str = ", ".join(
            f"{s}→{t} ({m.n_edges} edges)" for (s, t), m in self._mappings.items()
        )
        return (
            f"HierarchicalMuData with {self._mdata.n_obs} obs\n"
            f"  modalities: {mod_str}\n"
            f"  mappings: {mapping_str or 'none'}"
        )
