API Reference
=============

gcen2.preprocessing
-------------------

Data normalization and filtering.

.. automodule:: gcen2.preprocessing
   :members:
   :undoc-members: False
   :show-inheritance:

gcen2.network
-------------

Network construction and utilities.

.. automodule:: gcen2.network
   :members:
   :undoc-members: False
   :show-inheritance:

gcen2.module
------------

Gene module identification.

.. automodule:: gcen2.module
   :members:
   :undoc-members: False
   :show-inheritance:

gcen2.tom_modules
-----------------

WGCNA-style modules from the topological overlap (``gcen2 module -e``):
soft-thresholding power, TOM, tree cut, eigengenes, kME, trait correlations.

.. automodule:: gcen2.tom_modules
   :members:
   :undoc-members: False
   :show-inheritance:

gcen2.annotation
----------------

GO and KEGG functional annotation.

.. automodule:: gcen2.annotation
   :members:
   :undoc-members: False
   :show-inheritance:

gcen2.rwr
---------

Random Walk with Restart for gene function prediction.

.. automodule:: gcen2.rwr
   :members:
   :undoc-members: False
   :show-inheritance:

gcen2.single_cell
-----------------

Single-cell analysis utilities.

.. automodule:: gcen2.single_cell.core
   :members:
   :undoc-members: False
   :show-inheritance:

.. automodule:: gcen2.single_cell.metacell
   :members: metacells

Cell-type modules (the hdWGCNA workflow, ``gcen2 celltype-modules``),
differential module eigengenes and module-trait correlations, and module
projection:

.. automodule:: gcen2.single_cell.modules
   :members:

.. automodule:: gcen2.single_cell.dme
   :members:

.. automodule:: gcen2.single_cell.projection
   :members:

gcen2.spatial
-------------

Spatial transcriptomics: spatial modules (the hdWGCNA spatial workflow,
``gcen2 spatial-modules``), networks with spatial information
(``gcen2 spatial-build``) and Moran's I.

.. automodule:: gcen2.spatial.modules
   :members:

.. automodule:: gcen2.spatial.core
   :members:

.. automodule:: gcen2.spatial.neighbors
   :members:

gcen2.plot
----------

Publication figures for every analysis (one style, journal widths, editable vector text).

.. automodule:: gcen2.plot
   :members:

gcen2.consensus
---------------

Consensus modules of several expression sets (WGCNA blockwiseConsensusModules).

.. automodule:: gcen2.consensus
   :members: consensus_tom_modules, HDWGCNA_CONSENSUS

gcen2.preservation
------------------

Module preservation (WGCNA modulePreservation).

.. automodule:: gcen2.preservation
   :members: module_preservation

gcen2.harmony
-------------

Harmony batch correction of an embedding (used for harmonized module eigengenes).

.. automodule:: gcen2.harmony
   :members: harmony_integrate

gcen2.trajectory
----------------

Co-expression networks in sliding windows along pseudotime.

.. automodule:: gcen2.trajectory
   :members:
   :undoc-members: False
   :show-inheritance:

gcen2.multilayer
----------------

Multi-omics network integration: weighted fusion and multiplex RWR.

.. automodule:: gcen2.multilayer
   :members:
   :undoc-members: False
   :show-inheritance:

gcen2.differential
------------------

Differential co-expression analysis.

.. automodule:: gcen2.differential
   :members:
   :undoc-members: False
   :show-inheritance:

gcen2.plotting
--------------

Visualization functions.

.. automodule:: gcen2.plotting
   :members:
   :undoc-members: False
   :show-inheritance:

gcen2.report
------------

Automated research report generation.

.. automodule:: gcen2.report
   :members:
   :undoc-members: False
   :show-inheritance:

gcen2.models
------------

Core data models.

.. automodule:: gcen2.models
   :members:
   :undoc-members: False
   :show-inheritance:

gcen2.io
--------

File I/O utilities.

.. automodule:: gcen2.io
   :members:
   :undoc-members: False
   :show-inheritance:

gcen2.exceptions
----------------

Exception hierarchy.

.. automodule:: gcen2.exceptions
   :members:
   :show-inheritance:

gcen2.config
------------

Package configuration.

.. automodule:: gcen2.config
   :members:
