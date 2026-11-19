Small synthetic GO annotation fixture with the output of a compiled GCEN v1
0.6.3 `annotate`:

    annotate -g go.obo -a gene_go.assoc -n net.network  -o v1_network --pval 1 -t 1
    annotate -g go.obo -a gene_go.assoc -m modules.txt  -o v1_module  --pval 1 -t 1

(`--pval 1` keeps every term with p < 1; v1 0.6.3 does not accept the short
`-p` flag.) `tests/test_annotation.py::TestV1AnnotateFixture` checks that GCEN2
reproduces the queries, counts, names and namespaces exactly, and that the
two-sided Fisher p-values computed from GCEN2's counts match v1's.
