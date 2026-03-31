"""Visualization functions for hierarchical MuData structures."""

import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import mudata as md
import networkx as nx 
from scipy.sparse import triu
import seaborn as sns 

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


def plot_peptide_intensities(
    mudata: md.MuData,
    obs_group: str,
    level: str = "precursors",
    varp_key: str = "feature_mapping",
    figsize: tuple = (7, 8),
):
    """
    Generate a box plot of intensities between multiple conditions for all entities in an anndata object.

    Parameters:
        mudata: MuData
            Filtered MuData object, containing intensity data at different levels.
        obs_group: str
            The observation group to plot (must be a column name in mudata[level].obs).
        level:
            Name of the level in `mudata` that should be visualized. Default: "precursors".
    """
    fig, (ax_graph, ax_box) = plt.subplots(2, 1, figsize=figsize)
    plt.subplots_adjust(hspace=0)

    ### visualize graph on top
    # create a directed graph from the triangular adjacency matrix
    adj_directed = triu(mudata.varp[varp_key], k=-1).tocsr()
    G = nx.from_scipy_sparse_array(adj_directed, create_using=nx.DiGraph)
    # label the nodes according to the index of mudata.var
    mapping = {i: mudata.var.index[i] for i in range(mudata.varp[varp_key].shape[0])}
    G = nx.relabel_nodes(G, mapping, copy=False)

    pos = nx.nx_agraph.graphviz_layout(G, prog="dot", args="-Grankdir=BT")
    nx.draw_networkx(G, pos, with_labels=True, node_size=1000, node_color="grey", arrows=False, ax=ax_graph)
    ax_graph.axis("off")

    ### boxblot below
    # use dummy data instead of what's in mudata["precursors"] because the example data only contains 0 intensity values
    # TODO: Fix
    precursor_order = np.argsort(np.array(list(pos.values())[:mudata[level].shape[1]])[:, 0])
    df = mudata[level][:, precursor_order].to_df() # reorder anndata so that the boxes are in the same order as the leaves
    df["group"] = mudata[level].obs[obs_group]
    df_melt = df.melt(id_vars="group", var_name="peptide", value_name="intensity")

    sns.boxplot(x="peptide", y="intensity", hue="group", data=df_melt, ax=ax_box)
    ax_box.tick_params("x", rotation=90)

    return fig