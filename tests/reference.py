"""Reference implementations the tests compare GCEN2 with."""

from scipy.stats import hypergeom


def hypergeometric_pvalue(pop_n: int, pop_k: int, study_n: int, study_k: int) -> float:
    """P(X >= study_k) for X ~ Hypergeometric(pop_n, pop_k, study_n), by scipy."""
    return float(hypergeom.sf(study_k - 1, pop_n, pop_k, study_n))
