"""Visualization functions for hierarchical MuData structures."""

import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import mudata as md


def show_adjacency_matrix(mdata: md.MuData) -> None:
    """Visualize the adjacency matrix from MuData's varp.

    Parameters
    ----------
    mdata : md.MuData
        MuData object containing the adjacency matrix in varp

    Examples
    --------
    >>> from src.synthetic_data import generate_synthetic_mudata
    >>> mdata = generate_synthetic_mudata()
    >>> show_adjacency_matrix(mdata)
    """
    # Get the adjacency matrix from varp
    adj_matrix = mdata.varp['feature_mapping']

    # Convert to dense array if sparse
    if sparse.issparse(adj_matrix):
        adj_array = adj_matrix.toarray()
    else:
        adj_array = np.array(adj_matrix)

    # Get feature names
    features = list(mdata.var_names)

    # Create the figure
    plt.figure(figsize=(7, 7))

    # Create heatmap
    plt.imshow(adj_array, cmap='viridis', aspect='auto')

    # Set ticks and labels
    plt.xticks(ticks=np.arange(len(features)), labels=features, rotation=90)
    plt.yticks(ticks=np.arange(len(features)), labels=features)

    # Set title and labels
    plt.title('Feature Association Adjacency Matrix')
    plt.xlabel('Features')
    plt.ylabel('Features')

    # Adjust layout and show
    plt.tight_layout()
    plt.show()