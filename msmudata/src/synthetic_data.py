"""Synthetic data generation for hierarchical proteomics MuData."""

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse import csr_matrix
import anndata as ad
import mudata as md

from .adjancency_matrix_construct import get_unique_mappings, adjacency_matrix_from_mapping


def generate_synthetic_mudata(
    n_samples: int = 5,
    n_genes: int = 2,
    n_proteins: int = 3,
    n_precursors: int = 6,
    sparse_matrix: bool = True,
    random_seed: int = 42
) -> md.MuData:
    """Generate example hierarchical proteomics MuData with feature mappings.

    Creates a MuData object with three modalities (genes, proteins, precursors)
    and a feature mapping adjacency matrix in varp.

    Parameters
    ----------
    n_samples : int, default=5
        Number of samples (observations)
    n_genes : int, default=2
        Number of genes
    n_proteins : int, default=3
        Number of proteins
    n_precursors : int, default=6
        Number of precursors
    sparse_matrix : bool, default=True
        Whether to use sparse matrix for adjacency matrix
    random_seed : int, default=42
        Random seed for reproducible data generation

    Returns
    -------
    md.MuData
        MuData object with three modalities and feature mappings

    Examples
    --------
    >>> mdata = generate_synthetic_mudata()
    >>> print(mdata)
    MuData object with n_obs × n_vars = 5 × 11
      varp:	'feature_mapping'
      3 modalities
        precursors:	5 x 6
        proteins:	5 x 3
        genes:	5 x 2
    """
    np.random.seed(random_seed)

    # Create hierarchical mapping for precursors
    var_mapping = pd.DataFrame(
        {
            "genes": [*(["gene0", "gene1", "gene1"]*2)][:n_precursors],
            "proteins": [*(["protein0", "protein1", "protein2"]*2)][:n_precursors],
            "precursors": [str(idx) for idx in range(1, n_precursors + 1)]
        },
        index=[str(idx) for idx in range(1, n_precursors + 1)]
    )

    # Create precursors AnnData with some random data
    precursors = ad.AnnData(
        X=np.random.randn(n_samples, n_precursors),
        obs=pd.DataFrame(index=[f"sample{idx}" for idx in range(n_samples)]),
        var=var_mapping
    )

    # Create proteins AnnData
    protein_names = [f"protein{i}" for i in range(n_proteins)]
    gene_mapping = ["gene0", "gene1", "gene1"][:n_proteins]

    proteins = ad.AnnData(
        X=np.random.randn(n_samples, n_proteins),
        obs=pd.DataFrame(index=[f"sample{idx}" for idx in range(n_samples)]),
        var=pd.DataFrame(
            data={
                "proteins": protein_names,
                "genes": gene_mapping
            },
            index=protein_names
        )
    )

    # Create genes AnnData
    gene_names = [f"gene{i}" for i in range(n_genes)]

    genes = ad.AnnData(
        X=np.random.randn(n_samples, n_genes),
        obs=pd.DataFrame(index=[f"sample{idx}" for idx in range(n_samples)]),
        var=pd.DataFrame(index=gene_names)
    )

    # Create feature mapping for adjacency matrix
    levels = ["genes", "proteins", "precursors"]
    psm_subset = pd.concat({"A": var_mapping, "B": var_mapping}).reset_index(level=0, names="sample")

    df = get_unique_mappings(psm_subset, feature_level_names=levels)
    adjacency = adjacency_matrix_from_mapping(df).sort_index(axis=0).sort_index(axis=1)

    # Create MuData
    mdata = md.MuData(
        data={"precursors": precursors, "proteins": proteins, "genes": genes},
    )

    # Add adjacency matrix to varp
    adjacency_reindexed = adjacency.reindex(index=mdata.var_names, columns=mdata.var_names)

    if sparse_matrix:
        mdata.varp["feature_mapping"] = csr_matrix(adjacency_reindexed.values)
    else:
        mdata.varp["feature_mapping"] = adjacency_reindexed.values

    return mdata


def generate_simple_test_mudata() -> md.MuData:
    """Generate a minimal test MuData for unit testing.

    Returns
    -------
    md.MuData
        Minimal MuData with known structure for testing
    """
    return generate_synthetic_mudata(
        n_samples=3,
        n_genes=2,
        n_proteins=2,
        n_precursors=4,
        sparse_matrix=True,
        random_seed=0
    )