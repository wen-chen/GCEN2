"""Every `gcen2 plot` command, run on the outputs of the golden pipeline
(tests/golden); a new plot command fails this test until it is added here."""

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("matplotlib")
pytest.importorskip("anndata")

from gcen2.core import _HAVE_CORE  # noqa: E402

pytestmark = pytest.mark.skipif(not _HAVE_CORE, reason="compiled extension not available")


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    from tests.golden.golden import INPUTS, run_steps

    out = tmp_path_factory.mktemp("plots")
    run_steps(out)
    import anndata as ad

    cells = ad.read_h5ad(INPUTS / "sc.h5ad", backed="r").obs_names
    xy = np.random.default_rng(0).normal(size=(len(cells), 2))
    pd.DataFrame(xy, index=cells, columns=["x", "y"]).to_csv(out / "umap.tsv", sep="\t", lineterminator="\n")
    return out, INPUTS


def _cases(o, i):
    o, i = str(o), str(i)
    cta = f"{o}/ctm/A"
    return {
        "soft-power": ["-i", f"{o}/wgcna_soft_power.tsv"],
        "dendrogram": ["-i", f"{o}/wgcna_dendrogram.json"],
        "eigengene-network": ["-i", f"{cta}/eigengenes.tsv"],
        "module-trait": ["-i", f"{o}/wgcna_trait_cor.tsv"],
        "kme": ["-i", f"{cta}/kme.tsv"],
        "mm-gs": ["-i", f"{o}/wgcna_kme.tsv", "--gs", f"{o}/wgcna_gs.tsv", "-m", "ME1",
                  "--trait", "weight"],
        "hub-network": ["-t", f"{cta}/tom_hubs.tsv", "-k", f"{cta}/kme.tsv"],
        "hub-network-all": ["-t", f"{cta}/tom_hubs.tsv", "-k", f"{cta}/kme.tsv"],
        "module-umap": ["-t", f"{cta}/tom_hubs.tsv", "-k", f"{cta}/kme.tsv"],
        "tom-heatmap": ["-e", f"{cta}/metacells.tsv", "-m", f"{cta}/modules.tsv",
                        "--power", "12"],
        "sample-tree": ["-e", f"{o}/filt.tsv"],
        "module-sizes": ["-m", f"{cta}/modules.tsv"],
        "me-embedding": ["-i", f"{cta}/cell_eigengenes.tsv", "-c", f"{o}/umap.tsv"],
        "me-dotplot": ["-i", f"{cta}/cell_eigengenes.tsv", "-g", f"{i}/sc.h5ad",
                       "--groupby", "batch"],
        "me-violin": ["-i", f"{cta}/cell_eigengenes.tsv", "-g", f"{i}/sc.h5ad",
                      "--groupby", "batch"],
        "me-radar": ["-i", f"{o}/spm/spot_eigengenes.tsv", "-g", f"{i}/sp.h5ad",
                     "--groupby", "region"],  # a radar needs 3 or more groups
        "dme-volcano": ["-i", f"{cta}/dmes.tsv"],
        "dme-lollipop": ["-i", f"{cta}/dmes.tsv"],
        "celltype-overview": ["-i", f"{o}/ctm/celltypes.tsv"],
        "metacells": ["-i", f"{o}/ctm", "-c", f"{o}/umap.tsv"],
        "set-eigengenes": ["-i", f"{o}/ctm_cons/A/set_eigengenes.tsv"],
        "preservation": ["-i", f"{o}/preservation/preservation_Z.tsv"],
        "moran": ["-i", f"{o}/sp_moran_w.tsv"],
        "expression-qc": ["-e", f"{i}/bulk.tsv"],
        "mean-variance": ["-e", f"{o}/norm_tmm.tsv"],
        "sample-correlation": ["-e", f"{o}/filt.tsv"],
        "degree-distribution": ["-i", f"{o}/pear.network"],
        "edge-weights": ["-i", f"{o}/spear.network"],
        "network-view": ["-i", f"{o}/pear.network"],
        "module-size-hist": ["-m", f"{o}/jaccard.txt"],
        "celltype-networks": ["-i", f"{o}/celltype"],
        "null-distribution": ["-i", f"{o}/pear.network", "-s", f"{o}/shuf.network"],
        "merge-overlap": ["-i", f"{o}/pear.network", "-i", f"{o}/spear.network", "-m",
                          f"{o}/merged.network"],
        "extract-overview": ["-i", f"{o}/pear.network", "-g", f"{i}/seeds.txt"],
        "enrich-bar": ["-i", f"{o}/enrich.tsv"],
        "enrich-dot": ["-i", f"{o}/annot_m"],
        "enrich-heatmap": ["-i", f"{o}/annot_m"],
        "rwr-rank": ["-i", f"{o}/rwr.tsv"],
        "rwr-network": ["-n", f"{o}/pear.network", "-i", f"{o}/rwr.tsv", "-g",
                        f"{i}/seeds.txt"],
        "diff-overview": ["-i", f"{o}/diff.tsv"],
        "diff-network": ["-i", f"{o}/diff.tsv"],
        "trajectory": ["-i", f"{o}/trajectory"],
        "fusion": ["-i", f"{o}/fused.tsv"],
    }


def _plot_commands():
    import typer

    from gcen2.cli import plot_app

    return sorted(typer.main.get_command(plot_app).commands)


def test_every_plot_command_is_tested(tmp_path):
    assert sorted(_cases(tmp_path, tmp_path)) == _plot_commands()


@pytest.mark.parametrize("name", _plot_commands())
def test_plot_command(run, name, tmp_path):
    from typer.testing import CliRunner

    from gcen2.cli import app

    out, inputs = run
    args = _cases(out, inputs)[name]
    r = CliRunner().invoke(app, ["plot", name, *args, "-o", str(tmp_path / "f"), "--format",
                                 "png", "--format", "svg"])
    assert r.exit_code == 0, r.output
    # hub-network: one figure per module (f_ME<k>)
    for ext in ("png", "svg"):
        files = sorted(tmp_path.glob(f"f*.{ext}"))
        assert files and all(f.stat().st_size > 1000 for f in files), ext
