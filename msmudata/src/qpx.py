"""Load QPX Parquet datasets into a MuData object.

Each modality keeps its own quantitative obs x var matrix in ``AnnData.X``.
Cross-modality relationships are stored in ``mdata.varp["feature_mapping"]``
as a symmetric sparse adjacency matrix over all var IDs.

Modalities produced
-------------------
precursors
    obs = runs, var = unique peptidoform|charge, X = LFQ intensity
proteins
    obs = protein-level runs from pg.parquet, var = unique anchor_protein,
    X = LFQ intensity

Cross-modality mapping
----------------------
mdata.varp["feature_mapping"]
    Symmetric boolean CSR matrix, shape ``(n_vars, n_vars)``, indexed by
    ``mdata.var_names`` (precursor IDs followed by anchor protein IDs).
    A non-zero entry at ``[i, j]`` indicates that feature ``i`` and feature
    ``j`` are directly connected in the hierarchy.  Every precursor has
    exactly one protein edge; the matrix is symmetric so the inverse
    direction is implicit.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import polars as pl
import scipy.sparse as sp


_PRECURSOR_ID_EXPR = (
    pl.col("peptidoform") + "|" + pl.col("charge").cast(pl.String)
).alias("precursor_id")


def from_qpx(
    data_dir: Path | str,
    *,
    intensity_label: str = "LFQ",
) -> md.MuData:
    """Build a MuData object from a QPX Parquet directory.

    Cross-modality relationships are stored in ``mdata.varp["feature_mapping"]``
    as a symmetric sparse adjacency matrix indexed by ``mdata.var_names``.
    """
    data_dir = Path(data_dir)

    run_index = _load_run_index(data_dir)

    precursor_adata = _build_precursor_adata(data_dir, run_index, intensity_label)
    protein_adata = _build_protein_adata(data_dir, intensity_label)

    mdata = md.MuData(
        {
            "precursors": precursor_adata,
            "proteins": protein_adata,
        }
    )

    mdata.varp["feature_mapping"] = _build_feature_mapping_varp(mdata, precursor_adata)

    return mdata


def _find_parquet(data_dir: Path, table: str) -> Path:
    matches = sorted(data_dir.glob(f"*.{table}.parquet"))
    if not matches:
        raise FileNotFoundError(f"No *.{table}.parquet found in {data_dir}")
    return matches[0]


def _load_run_index(data_dir: Path) -> pd.Index:
    run_df = pl.read_parquet(_find_parquet(data_dir, "run"))
    return pd.Index(
        sorted(run_df["run_file_name"].drop_nulls().unique().to_list()),
        name="run_file_name",
    )


def _load_run_obs(data_dir: Path, run_index: pd.Index) -> pd.DataFrame:
    """Load scalar run metadata aligned to ``run_index``."""
    run_df = pl.read_parquet(_find_parquet(data_dir, "run"))

    nested_types = (pl.List, pl.Struct, pl.Array)
    simple_cols = [
        col_name
        for col_name in run_df.columns
        if col_name != "run_file_name"
        and not isinstance(run_df.schema[col_name], nested_types)
    ]

    return (
        run_df.select(["run_file_name", *simple_cols])
        .to_pandas()
        .set_index("run_file_name")
        .reindex(run_index)
    )


def _extract_intensity_label(lf: pl.LazyFrame, label: str) -> pl.LazyFrame:
    return (
        lf.explode("intensities")
        .filter(pl.col("intensities").struct.field("label") == label)
        .with_columns(
            pl.col("intensities").struct.field("intensity").alias("intensity")
        )
        .drop("intensities")
    )


def _long_to_csr(
    lf: pl.LazyFrame,
    *,
    run_col: str,
    var_col: str,
    value_col: str,
    run_index: pd.Index,
    var_index: pd.Index,
    dtype: type = np.float32,
) -> sp.csr_matrix:
    run_idx_df = pl.DataFrame(
        {run_col: list(run_index), "__row__": list(range(len(run_index)))}
    )
    var_idx_df = pl.DataFrame(
        {var_col: list(var_index), "__col__": list(range(len(var_index)))}
    )

    df = (
        lf.collect()
        .join(run_idx_df, on=run_col, how="inner")
        .join(var_idx_df, on=var_col, how="inner")
        .select(["__row__", "__col__", value_col])
        .group_by(["__row__", "__col__"])
        .agg(pl.sum(value_col))
    )

    row = df["__row__"].to_numpy()
    col = df["__col__"].to_numpy()
    values = df[value_col].cast(pl.Float64).to_numpy().astype(dtype)

    return sp.coo_matrix(
        (values, (row, col)),
        shape=(len(run_index), len(var_index)),
        dtype=dtype,
    ).tocsr()


def _build_precursor_adata(
    data_dir: Path, run_index: pd.Index, intensity_label: str
) -> ad.AnnData:
    feat_path = _find_parquet(data_dir, "feature")

    var_df = (
        pl.scan_parquet(feat_path)
        .with_columns(_PRECURSOR_ID_EXPR)
        .group_by("precursor_id")
        .agg(
            [
                pl.first("sequence"),
                pl.first("peptidoform"),
                pl.first("charge"),
                pl.first("missed_cleavages"),
                pl.first("calculated_mz"),
                pl.first("anchor_protein"),
                pl.first("unique"),
            ]
        )
        .collect()
        .sort("precursor_id")
        .to_pandas()
        .set_index("precursor_id")
    )
    var_df.index.name = "precursor_id"
    var_index = pd.Index(var_df.index, name="precursor_id")

    intensity_lf = (
        pl.scan_parquet(feat_path)
        .with_columns(_PRECURSOR_ID_EXPR)
        .select(["run_file_name", "precursor_id", "intensities"])
        .pipe(_extract_intensity_label, intensity_label)
        .select(["run_file_name", "precursor_id", "intensity"])
    )
    X = _long_to_csr(
        intensity_lf,
        run_col="run_file_name",
        var_col="precursor_id",
        value_col="intensity",
        run_index=run_index,
        var_index=var_index,
    )

    return ad.AnnData(X=X, obs=_load_run_obs(data_dir, run_index), var=var_df)


def _build_protein_adata(data_dir: Path, intensity_label: str) -> ad.AnnData:
    pg_path = _find_parquet(data_dir, "pg")
    pg_run_index = pd.Index(
        sorted(
            pl.read_parquet(pg_path)["run_file_name"]
            .drop_nulls()
            .unique()
            .to_list()
        ),
        name="run_file_name",
    )

    var_df = (
        pl.scan_parquet(pg_path)
        .filter(pl.col("anchor_protein").is_not_null())
        .group_by("anchor_protein")
        .agg(
            [
                pl.first("pg_names").list.join("; ").alias("pg_names"),
                pl.first("gg_accessions").list.join("; ").alias("gg_accessions"),
                pl.first("gg_names").list.join("; ").alias("gg_names"),
                pl.first("molecular_weight"),
                pl.first("sequence_coverage"),
                pl.first("global_qvalue"),
            ]
        )
        .collect()
        .sort("anchor_protein")
        .to_pandas()
        .set_index("anchor_protein")
    )
    var_df.index.name = "anchor_protein"
    var_index = pd.Index(var_df.index, name="anchor_protein")

    intensity_lf = (
        pl.scan_parquet(pg_path)
        .filter(pl.col("anchor_protein").is_not_null())
        .select(["run_file_name", "anchor_protein", "intensities"])
        .pipe(_extract_intensity_label, intensity_label)
        .select(["run_file_name", "anchor_protein", "intensity"])
    )
    X = _long_to_csr(
        intensity_lf,
        run_col="run_file_name",
        var_col="anchor_protein",
        value_col="intensity",
        run_index=pg_run_index,
        var_index=var_index,
    )

    return ad.AnnData(X=X, obs=pd.DataFrame(index=pg_run_index), var=var_df)


def _build_feature_mapping_varp(
    mdata: md.MuData,
    precursor_adata: ad.AnnData,
) -> pd.DataFrame:
    """Build a symmetric feature adjacency indexed by ``mdata.var_names``.

    Connects each precursor var to its anchor protein var with a 1
    (and the symmetric reverse edge).
    """
    mapping_df = (
        precursor_adata.var[["anchor_protein"]]
        .dropna()
        .reset_index()[["precursor_id", "anchor_protein"]]
    )

    all_vars = mdata.var_names
    var_to_idx = pd.Series(range(len(all_vars)), index=all_vars)
    n = len(all_vars)

    src_idx = mapping_df["precursor_id"].map(var_to_idx)
    tgt_idx = mapping_df["anchor_protein"].map(var_to_idx)
    valid = src_idx.notna() & tgt_idx.notna()
    i = src_idx[valid].astype(int).to_numpy()
    j = tgt_idx[valid].astype(int).to_numpy()

    half = sp.coo_matrix(
        (np.ones(len(i), dtype=np.float64), (i, j)), shape=(n, n)
    ).tocsr()
    adj = (half + half.T).tocsr()

    return pd.DataFrame.sparse.from_spmatrix(adj, index=all_vars, columns=all_vars)
