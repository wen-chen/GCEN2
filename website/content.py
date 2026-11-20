"""Content of the GCEN2 website that is not kept elsewhere in the repository.

Facts that are kept elsewhere (version, benchmark table, number of commands,
citation) are read by ``build.py`` from the repository itself.
"""

# True once GCEN2 is on PyPI and Bioconda. False shows a pre-release notice
# and gives the installation from source first.
RELEASED = True

SITE_URL = "https://www.biochen.org/gcen2/"
GITHUB = "https://github.com/wen-chen/GCEN2"
GCEN_V1_SITE = "https://www.biochen.org/gcen"
GCEN_V1_GITHUB = "https://github.com/wen-chen/GCEN"
BIOCHEN = "https://www.biochen.org/"
CONTACT_EMAIL = "chenwen@biochen.org"

# Analysis kinds: the lanes of the workflow figure, with its colours.
KINDS = {
    "networks": ("Networks", "#0072B2"),
    "wgcna": ("WGCNA modules", "#56B4E9"),
    "single-cell": ("Single cell", "#D55E00"),
    "spatial": ("Spatial", "#009E73"),
    "function": ("Function", "#E69F00"),
    "multiomics": ("Trajectories & multi-omics", "#B4527F"),
}

CAPABILITIES = [
    {
        "kind": "networks",
        "title": "Co-expression networks",
        "text": "Everything GCEN v1 does, with the same commands, files and results: "
                "normalization, filtering, Pearson or Spearman networks with p-values and "
                "FDR, modules, merge, extract, shuffle, statistics. Differential "
                "co-expression between two conditions with DGCA's test.",
        "commands": ["normalize", "filter", "build", "module", "diff"],
        "tutorial": "docs/tutorials/bulk.html",
    },
    {
        "kind": "wgcna",
        "title": "WGCNA modules",
        "text": "Soft power, topological overlap, dynamic tree cut, eigengenes, kME, hub "
                "genes and module–trait correlations, consensus modules over several data "
                "sets and module preservation, giving the same modules as WGCNA.",
        "commands": ["module -e", "consensus-modules", "module-preservation"],
        "tutorial": "docs/tutorials/bulk.html",
    },
    {
        "kind": "single-cell",
        "title": "Single cells",
        "text": "The hdWGCNA workflow from an AnnData file: metacells per cell type, "
                "modules, cell eigengenes, Harmony, differential module eigengenes and "
                "module projection; networks per cell type with p-values from the "
                "effective number of samples. Atlas scale, from backed files.",
        "commands": ["celltype-modules", "celltype", "project-modules"],
        "tutorial": "docs/tutorials/single_cell.html",
    },
    {
        "kind": "spatial",
        "title": "Spatial transcriptomes",
        "text": "Visium, Visium HD and Slide-seq: metaspots on the grid or by nearest "
                "neighbours, modules and spot eigengenes on the tissue, module–region "
                "correlations, spatially smoothed networks and Moran's I.",
        "commands": ["spatial-modules", "spatial-build"],
        "tutorial": "docs/tutorials/spatial.html",
    },
    {
        "kind": "multiomics",
        "title": "Trajectories and multi-omics",
        "text": "Networks along pseudotime with the edges gained, lost and transient; "
                "fusion of co-expression with prior networks such as STRING and TRRUST; "
                "random walk with restart on multiplex networks.",
        "commands": ["trajectory", "fuse", "rwr-multiplex"],
        "tutorial": "docs/tutorials/trajectory_multiomics.html",
    },
    {
        "kind": "function",
        "title": "Gene function, figures and reports",
        "text": "GO and KEGG annotation and enrichment, gene ranking from seed genes, "
                "{n_plots} kinds of publication figures (PDF, PNG, SVG) and HTML or PDF "
                "reports of whole analyses.",
        "commands": ["annotate", "enrich", "rwr", "plot", "report"],
        "tutorial": "docs/tutorials/figures.html",
    },
]

# Headline numbers on the home page: (row of the README benchmark table whose
# first cell starts with this text, the figure shown, what it means).
HEADLINES = [
    ("GCEN v1", "6–89×", "faster than GCEN v1, with the same edges and modules"),
    ("WGCNA", "4.0×", "faster than WGCNA on GTEx whole blood, 2.7× less memory"),
    ("hdWGCNA 0.4.12 | Human Lung", "8 min", "for all 49 cell types of a 584,944-cell atlas"),
    ("DGCA", "185×", "faster than DGCA for Spearman differential co-expression"),
]

# The gallery. Each figure is made by the setup commands of its group and
# then its own command, all run by gallery/make_gallery.py in one directory
# (EX, SC and SP point to data/example/). The page shows `command`, which is
# exactly what made the image (its PNG is the `-o` path plus ".png").
GALLERY_SETUP = {
    "bulk": [
        "gcen2 normalize -i $EX/gene_expr.tsv -o norm.tsv -m tmm",
        "gcen2 filter -i norm.tsv -o filt.tsv -p 0.75",
        "gcen2 build -i filt.tsv -o net.network -m spearman -c 0.8 -p 0.001 -f -t 8",
        "gcen2 module -i net.network -o module.txt -s 0.5 -t 8",
        "gcen2 rwr -n net.network -g $EX/rwr_seed_genes.list -o rwr.tsv",
        "gcen2 enrich -e $EX/enrichment_gene.list -b $EX/background_gene.list "
        "-g $EX/go-basic.obo.gz -a $EX/gene_go.assoc.gz -o enrich.tsv",
        "gcen2 filter -i norm.tsv -o var.tsv -p 0.75 --percent-sd 0.2",
        "gcen2 module -e var.tsv -o wgcna/modules.tsv -l log2 --signed -t 8 "
        "--me wgcna/eigengenes.tsv --kme wgcna/kme.tsv --hubs wgcna/hubs.tsv "
        "--soft-power wgcna/soft_power.tsv --dendrogram wgcna/dendrogram.json "
        "--tom-hubs wgcna/tom_hubs.tsv",
    ],
    "single-cell": [
        "gcen2 celltype-modules -i $SC --groupby cell_type -o modules/ -t 8 "
        "--write-metacells --write-tom --dme-by cluster --trait total_counts",
    ],
    "spatial": [
        "gcen2 spatial-modules -i $SP --groupby region -o spatial/ -t 8 --write-tom",
        "gcen2 spatial-build -i $SP -m weighted -o smoothed.network -c 0.6 -t 8 "
        "--moran moran_smoothed.tsv",
    ],
}

GALLERY = [
    # id, group, kind, title, description, command
    ("degree", "bulk", "networks", "Degree distribution",
     "Node degrees of the co-expression network on log–log axes.",
     "gcen2 plot degree-distribution -i net.network -o fig/degree"),
    ("network", "bulk", "networks", "Network view",
     "The network with genes coloured by module.",
     "gcen2 plot network-view -i net.network -m module.txt -o fig/network"),
    ("soft_power", "bulk", "wgcna", "Soft-threshold power",
     "Scale-free fit and mean connectivity for each power; the chosen power marked.",
     "gcen2 plot soft-power -i wgcna/soft_power.tsv -o fig/soft_power"),
    ("dendrogram", "bulk", "wgcna", "Gene dendrogram",
     "Hierarchical clustering of 1 − TOM with the module colours below.",
     "gcen2 plot dendrogram -i wgcna/dendrogram.json -o fig/dendrogram"),
    ("eigengenes", "bulk", "wgcna", "Eigengene network",
     "Correlations between module eigengenes, clustered.",
     "gcen2 plot eigengene-network -i wgcna/eigengenes.tsv -o fig/eigengenes"),
    ("rwr", "bulk", "function", "Random walk with restart",
     "Genes ranked by their proximity to the seed genes in the network.",
     "gcen2 plot rwr-rank -i rwr.tsv -g $EX/rwr_seed_genes.list -o fig/rwr"),
    ("enrich", "bulk", "function", "GO enrichment",
     "The most enriched GO terms of a gene list.",
     "gcen2 plot enrich-bar -i enrich.tsv -o fig/enrich"),
    ("celltypes", "single-cell", "single-cell", "Modules per cell type",
     "Metacells, genes, power and modules of each cell type.",
     "gcen2 plot celltype-overview -i modules/celltypes.tsv -o fig/celltypes"),
    ("b_eigengenes", "single-cell", "single-cell", "Cell eigengenes",
     "Module eigengenes of B cells on the UMAP.",
     "gcen2 plot me-embedding -i modules/B/cell_eigengenes.tsv -c $SC -o fig/b_eigengenes"),
    ("mono_dotplot", "single-cell", "single-cell", "Eigengenes by cluster",
     "Mean eigengene and the fraction of cells above zero in each cluster.",
     "gcen2 plot me-dotplot -i modules/CD14_Mono/cell_eigengenes.tsv -g $SC "
     "--groupby cluster -o fig/mono_dotplot"),
    ("mono_dmes", "single-cell", "single-cell", "Differential module eigengenes",
     "Modules whose eigengene differs between clusters.",
     "gcen2 plot dme-volcano -i modules/CD14_Mono/dmes.tsv -o fig/mono_dmes"),
    ("b_hub_me1", "single-cell", "single-cell", "Hub genes",
     "The hub genes of one B-cell module and their strongest TOM links.",
     "gcen2 plot hub-network -t modules/B/tom_hubs.tsv -k modules/B/kme.tsv -m ME1 "
     "-o fig/b_hub_me1"),
    ("spatial_me", "spatial", "spatial", "Spot eigengenes",
     "Module eigengenes on the tissue (mouse brain, Visium).",
     "gcen2 plot me-embedding -i spatial/spot_eigengenes.tsv -c $SP --basis spatial "
     "--spatial -o fig/spatial_me"),
    ("regions", "spatial", "spatial", "Modules and regions",
     "Correlation of each module with each brain region.",
     "gcen2 plot module-trait -i spatial/group_cor.tsv --annotate stars -o fig/regions"),
    ("spatial_kme", "spatial", "wgcna", "Module membership",
     "The genes of each spatial module ranked by kME, hub genes named.",
     "gcen2 plot kme -i spatial/kme.tsv -o fig/spatial_kme"),
    ("moran", "spatial", "spatial", "Spatial autocorrelation",
     "Moran's I of the genes, against permuted spots.",
     "gcen2 plot moran -i moran_smoothed.tsv -o fig/moran"),
]

FAQ = [
    ("What is GCEN2?",
     "<p>GCEN2 is a toolkit for gene co-expression network analysis: a command-line "
     "program, <code>gcen2</code>, and a Python package, with C++ kernels for speed. It "
     "covers bulk transcriptomes, single cells, spatial transcriptomes, trajectories and "
     "multi-omics networks, and makes publication figures and reports.</p>"),
    ("How does GCEN2 differ from GCEN?",
     "<p>GCEN2 is the successor of <a href=\"{v1_site}\">GCEN</a>. It keeps every GCEN "
     "command and file format, with the same results, and adds WGCNA and hdWGCNA module "
     "analysis, single-cell and spatial workflows, differential co-expression, trajectory "
     "and multi-omics networks, figures and reports. GCEN was a set of C++ programs; GCEN2 "
     "is one Python package with the same C++ speed where it matters.</p>"),
    ("Do my GCEN pipelines still work?",
     "<p>Yes, with the <code>gcen2</code> prefix: <code>network_build</code> becomes "
     "<code>gcen2 build</code>, <code>module_identify</code> becomes <code>gcen2 "
     "module</code>, and so on, with the same options and output files. The "
     "<a href=\"docs/migration.html\">migration guide</a> maps every command, including "
     "the utilities, and lists the few deliberate differences.</p>"),
    ("Command-line scares you off?",
     "<p>Each analysis is one command with a few options, and the "
     "<a href=\"docs/tutorials/bulk.html\">tutorials</a> run on example data that come "
     "with GCEN2. Everything is also available from Python for those who prefer "
     "notebooks.</p>"),
    ("What input does GCEN2 take?",
     "<p>Expression matrices as tab-separated text (genes in rows, samples in columns, "
     "as in GCEN), and single-cell or spatial data as AnnData (<code>.h5ad</code>) files, "
     "including Visium, Visium HD and Slide-seq. Networks, modules and annotations use "
     "GCEN's text formats. CSV files can be read in Python with "
     "<code>gcen2.read_csv</code> and written as TSV with <code>gcen2.write_tsv</code>.</p>"),
    ("Does it run on Windows?",
     "<p>Yes. GCEN2 is built and tested on Linux (x86-64), Windows (x86-64) and macOS "
     "(Apple Silicon). Use <code>pip</code> on Windows; Bioconda does not build Windows "
     "packages.</p>"),
    ("Are the results reproducible?",
     "<p>Yes. The results do not depend on the number of threads, and the test suite "
     "checks that they are the same on Linux, Windows and macOS. Random steps such as "
     "metacells take a seed.</p>"),
    ("How large can the data be?",
     "<p>The single-cell workflow ran on the 584,944 cells of the Human Lung Cell Atlas, "
     "all 49 cell types with at least 100 cells, in about eight minutes on one 16-core "
     "machine. Large <code>.h5ad</code> files can be read in backed mode, without loading "
     "the whole matrix.</p>"),
    ("Where do I get GO and KEGG annotations?",
     "<p>GO annotations for most model organisms are available from the "
     "<a href=\"https://geneontology.org/\">Gene Ontology Consortium</a>; KEGG "
     "annotations of protein-coding genes can be obtained with "
     "<a href=\"https://www.genome.jp/kegg/kaas/\">KAAS</a>, and InterProScan gives both "
     "for proteins. The GCEN <a href=\"{v1_site}/download\">download page</a> provides "
     "annotation files for human, mouse and zebrafish, which GCEN2 reads unchanged.</p>"),
    ("How do I make figures?",
     "<p><code>gcen2 plot</code> has {n_plots} figures for the files the analyses write, "
     "in journal column widths, as PDF, PNG and SVG with editable text; see the "
     "<a href=\"gallery.html\">gallery</a>. <code>gcen2 report</code> collects an analysis "
     "into one HTML or PDF report.</p>"),
    ("How do I cite GCEN2?",
     "<p>See <a href=\"index.html#cite\">Cite</a>. Please also cite WGCNA and hdWGCNA when "
     "you use their methods.</p>"),
    ("What if I can't solve a problem?",
     "<p>Open an issue on <a href=\"{github}/issues\">GitHub</a> or write to "
     "<a href=\"mailto:{email}\">{email}</a>. Questions and suggestions are welcome.</p>"),
]
