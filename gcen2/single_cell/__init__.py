"""Single-cell specific functionality for GCEN2."""

from gcen2.single_cell.core import (
    adata_matrix,
    build_by_celltype,
    build_dropout_aware,
    effective_samples,
    is_anndata,
    metacell_effective_samples,
    pseudobulk,
    write_celltype_networks,
)
from gcen2.single_cell.dme import find_all_dmes, find_dmes, module_trait_cor
from gcen2.single_cell.metacell import metacells
from gcen2.single_cell.modules import (
    cell_eigengenes,
    expressed_fraction,
    fraction_genes,
    iter_modules_by_celltype,
    modules_by_celltype,
    write_celltype_modules,
)
from gcen2.single_cell.projection import project_modules

__all__ = ["adata_matrix", "is_anndata", "metacells", "pseudobulk", "build_by_celltype",
           "build_dropout_aware", "effective_samples", "metacell_effective_samples",
           "write_celltype_networks",
           "cell_eigengenes", "expressed_fraction", "fraction_genes", "iter_modules_by_celltype",
           "modules_by_celltype", "write_celltype_modules", "find_dmes", "find_all_dmes",
           "module_trait_cor", "project_modules"]
