"""Synthetic MS proteomics data generator for testing and demos.

Creates realistic-ish hierarchical data: precursors → proteins with N:M
relationships, mirroring the structure of real LC/MS proteomics datasets.
"""

from __future__ import annotations

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd

from .feature_mapping import FeatureMapping
from .hierarchical_mudata import HierarchicalMuData


def make_synthetic_proteomics_data(
    n_samples: int = 6,
    n_proteins: int = 200,
    n_precursors_per_protein: tuple[int, int] = (3, 15),
    shared_precursor_fraction: float = 0.05,
    missing_value_rate: float = 0.1,
    seed: int = 42,
) -> HierarchicalMuData:
    """Generate a synthetic HierarchicalMuData mimicking LC/MS proteomics.

    Parameters
    ----------
    n_samples
        Number of samples (runs).
    n_proteins
        Number of distinct proteins.
    n_precursors_per_protein
        (min, max) precursors per protein, drawn uniformly.
    shared_precursor_fraction
        Fraction of precursors that map to >1 protein (creates N:M edges).
    missing_value_rate
        Proportion of intensity values set to NaN.
    seed
        Random seed.

    Returns
    -------
    A :class:`HierarchicalMuData` with ``"proteins"`` and ``"precursors"``
    modalities and a registered mapping between them.
    """
    rng = np.random.default_rng(seed)

    protein_ids = [f"PROT_{i:04d}" for i in range(n_proteins)]
    sample_ids = [f"sample_{i:02d}" for i in range(n_samples)]

    precursor_list: list[str] = []
    edges: list[dict[str, str]] = []
    precursor_counter = 0

    for prot_id in protein_ids:
        n_prec = rng.integers(
            n_precursors_per_protein[0], n_precursors_per_protein[1] + 1
        )
        for _ in range(n_prec):
            prec_id = f"prec_{precursor_counter:06d}"
            precursor_list.append(prec_id)
            edges.append({"source": prec_id, "target": prot_id})
            precursor_counter += 1

    n_shared = int(len(precursor_list) * shared_precursor_fraction)
    shared_indices = rng.choice(len(precursor_list), size=n_shared, replace=False)
    for idx in shared_indices:
        prec_id = precursor_list[idx]
        other_prot = rng.choice(protein_ids)
        edges.append({"source": prec_id, "target": other_prot})

    edges_df = pd.DataFrame(edges).drop_duplicates().reset_index(drop=True)

    n_precursors = len(precursor_list)

    protein_intensities = rng.lognormal(mean=20, sigma=2, size=(n_samples, n_proteins))
    protein_mask = rng.random((n_samples, n_proteins)) < missing_value_rate
    protein_intensities[protein_mask] = np.nan

    precursor_intensities = rng.lognormal(
        mean=18, sigma=2.5, size=(n_samples, n_precursors)
    )
    precursor_mask = rng.random((n_samples, n_precursors)) < missing_value_rate
    precursor_intensities[precursor_mask] = np.nan

    protein_adata = ad.AnnData(
        X=protein_intensities.astype(np.float32),
        obs=pd.DataFrame(index=sample_ids),
        var=pd.DataFrame(index=protein_ids),
    )

    precursor_adata = ad.AnnData(
        X=precursor_intensities.astype(np.float32),
        obs=pd.DataFrame(index=sample_ids),
        var=pd.DataFrame(index=precursor_list),
    )

    mdata = md.MuData({"proteins": protein_adata, "precursors": precursor_adata})

    mapping = FeatureMapping(
        source_index=pd.Index(precursor_list),
        target_index=pd.Index(protein_ids),
        edges=edges_df,
    )

    return HierarchicalMuData(
        mdata=mdata,
        mappings={("precursors", "proteins"): mapping},
    )
