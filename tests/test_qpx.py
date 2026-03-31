from pathlib import Path
import sys

import polars as pl
import pytest


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "msmudata"))


def _write_qpx_fixture(tmp_path: Path) -> Path:
    data_dir = tmp_path / "qpx"
    data_dir.mkdir()

    pl.DataFrame(
        {
            "run_file_name": ["run_1", "run_2"],
            "run_accession": ["run_1", "run_2"],
        }
    ).write_parquet(data_dir / "toy.run.parquet")

    pl.DataFrame(
        {
            "sequence": ["PEPA", "PEPA", "PEPB"],
            "peptidoform": ["PEPA", "PEPA", "PEPB"],
            "charge": [2, 2, 2],
            "missed_cleavages": [0, 0, 0],
            "calculated_mz": [500.0, 500.0, 600.0],
            "anchor_protein": ["PROT_A", "PROT_A", "PROT_B"],
            "unique": [True, True, True],
            "run_file_name": ["run_1", "run_2", "run_1"],
            "intensities": [
                [{"label": "LFQ", "intensity": 10.0}],
                [{"label": "LFQ", "intensity": 12.0}],
                [{"label": "LFQ", "intensity": 20.0}],
            ],
        }
    ).write_parquet(data_dir / "toy.feature.parquet")

    pl.DataFrame(
        {
            "run_file_name": ["run_1", "run_2", "run_1"],
            "anchor_protein": ["PROT_A", "PROT_A", "PROT_B"],
            "pg_names": [["Protein A"], ["Protein A"], ["Protein B"]],
            "gg_accessions": [["GENE_A"], ["GENE_A"], ["GENE_B"]],
            "gg_names": [["Gene A"], ["Gene A"], ["Gene B"]],
            "molecular_weight": [1000.0, 1000.0, 2000.0],
            "sequence_coverage": [0.5, 0.5, 0.7],
            "global_qvalue": [0.001, 0.001, 0.002],
            "intensities": [
                [{"label": "LFQ", "intensity": 100.0}],
                [{"label": "LFQ", "intensity": 110.0}],
                [{"label": "LFQ", "intensity": 200.0}],
            ],
        }
    ).write_parquet(data_dir / "toy.pg.parquet")

    return data_dir


@pytest.fixture()
def qpx_dir(tmp_path: Path) -> Path:
    return _write_qpx_fixture(tmp_path)


def test_qpx_builds_precursor_protein_mapping(qpx_dir: Path) -> None:
    from src.qpx import from_qpx

    mdata = from_qpx(qpx_dir)

    assert set(mdata.mod.keys()) == {"precursors", "proteins"}
    assert "feature_mapping" in mdata.varp

    fm = mdata.varp["feature_mapping"]  # CSR matrix, shape (n_vars, n_vars)
    n_vars = len(mdata.var_names)
    assert fm.shape == (n_vars, n_vars)

    # Build a lookup from var name → position in mdata.var_names
    var_pos = {v: i for i, v in enumerate(mdata.var_names)}

    def edge(src: str, tgt: str) -> float:
        return fm[var_pos[src], var_pos[tgt]]

    # Symmetric adjacency: precursor ↔ protein edges
    assert edge("PEPA|2", "PROT_A") == 1.0
    assert edge("PROT_A", "PEPA|2") == 1.0
    assert edge("PEPB|2", "PROT_B") == 1.0
    assert edge("PROT_B", "PEPB|2") == 1.0

    # No spurious edges
    assert edge("PEPA|2", "PROT_B") == 0.0
    assert edge("PEPA|2", "PEPB|2") == 0.0

    # 2 edges × 2 directions = 4 non-zeros
    assert fm.nnz == 4


def test_qpx_matrices_are_quantitative(qpx_dir: Path) -> None:
    from src.qpx import from_qpx

    mdata = from_qpx(qpx_dir)

    precursors = mdata["precursors"]
    assert precursors.X.shape == (2, 2)
    assert precursors.X.toarray().sum() == pytest.approx(42.0)

    proteins = mdata["proteins"]
    assert proteins.X.shape == (2, 2)
    assert proteins.X.toarray().sum() == pytest.approx(410.0)
