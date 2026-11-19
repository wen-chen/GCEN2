"""GCEN2 exception hierarchy."""


class GCENError(Exception):
    """Base exception for all GCEN2 errors."""


class GCENInputError(GCENError):
    """Raised for bad input files or invalid parameters."""


class GCENComputeError(GCENError):
    """Raised for numerical failures in C++ kernels."""


class GCENIOError(GCENError):
    """Raised for file read/write failures."""
