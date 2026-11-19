"""Publication figures for every GCEN2 analysis.

All plots share one style (:mod:`gcen2.plot._style`): journal widths
(``size="single"`` 89 mm, ``"onehalf"`` 120 mm, ``"double"`` 183 mm),
5-7 pt sans-serif text kept editable in PDF / SVG, and one colour per module
in every figure (WGCNA's standardColors by default; ``set_defaults(palette=
"okabe-ito")``). Each function takes GCEN2 result files or the objects
themselves, draws on a new figure (or the given ``ax``), saves it when
``output`` is given (a stem: PDF, PNG and SVG; or one file by its
extension) and returns the figure.

The figures are also available from the command line (``gcen2 plot ...``)
and in the HTML / PDF report.
"""

from gcen2.plot._save import save_figure
from gcen2.plot._style import (
    figsize,
    get_defaults,
    module_color,
    module_label,
    set_defaults,
    style,
)
from gcen2.plot.diff import diff_network, diff_overview
from gcen2.plot.enrichment import enrich_bar, enrich_dot, enrich_heatmap
from gcen2.plot.function import rwr_rank, rwr_subnetwork
from gcen2.plot.multiomics import fusion_overview
from gcen2.plot.network import (
    celltype_networks,
    degree_distribution,
    edge_weights,
    extract_overview,
    merge_overlap,
    module_size_hist,
    network_view,
    null_distribution,
)
from gcen2.plot.qc import expression_distribution, mean_variance, sample_correlation
from gcen2.plot.single_cell import (
    celltype_overview,
    dme_lollipop,
    dme_volcano,
    me_dotplot,
    me_embedding,
    me_radar,
    me_violin,
    metacells,
    moran,
    preservation,
    preservation_lollipop,
    set_eigengenes,
    spatial_me,
)
from gcen2.plot.trajectory import trajectory_overview
from gcen2.plot.wgcna import (
    dendrogram,
    eigengene_network,
    hub_network,
    hub_network_all,
    kme,
    mm_gs,
    module_sizes,
    module_trait,
    module_umap,
    sample_tree,
    soft_power,
    tom_heatmap,
)

__all__ = ["celltype_networks", "degree_distribution", "diff_network", "diff_overview",
           "edge_weights", "enrich_bar", "enrich_dot", "enrich_heatmap",
           "expression_distribution", "extract_overview", "fusion_overview", "mean_variance",
           "merge_overlap", "metacells", "module_size_hist", "moran",
           "network_view", "null_distribution", "rwr_rank", "rwr_subnetwork",
           "sample_correlation", "trajectory_overview",
           "celltype_overview", "dme_lollipop", "dme_volcano", "me_dotplot", "me_embedding",
           "me_radar", "me_violin", "preservation", "preservation_lollipop", "set_eigengenes",
           "spatial_me",
           "dendrogram", "eigengene_network", "figsize", "get_defaults", "hub_network",
           "hub_network_all", "kme", "mm_gs", "module_color", "module_label", "module_sizes",
           "module_trait", "module_umap", "sample_tree", "save_figure", "set_defaults",
           "soft_power", "style", "tom_heatmap"]
