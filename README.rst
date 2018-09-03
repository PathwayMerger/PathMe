ComPath Reloaded
================

The primary goal of this package is to convert KEGG, Reactome, and WikiPathways to Biological Expression Language (BEL). ComPath Reloaded is the continuation of the ComPath web application aimed at exploring, analyzing, and curating pathway knowledge in a more simplistic gene-centric view. This different approach involves converting all the pathways to BEL as a pivotal integration schema and evaluating consensus and gaps in pathway knowledge. Additionally, ComPath Reloaded is complemented with PathMe, a web application that enables the exploration of all the pathways from these resources using the mappings curated from ComPath.

Installation
------------
1. ``compath-reloaded`` can be installed with the following commmands:

.. code-block:: sh

    python3 -m pip install git+https://github.com/ComPath/ComPath-Reloaded.git@master

2. or in editable mode with:

.. code-block:: sh

    git clone https://github.com/ComPath/ComPath-Reloaded.git

.. code-block:: sh

    cd compath-reloaded

.. code-block:: sh

    python3 -m pip install -e .
    
How to use
----------

1. Download content

ComPath Reloaded first requires to download the raw files from the original pathway databases.

.. code-block:: python

    python3 -m compath_reloaded download
    
2. Generate BEL Graphs

.. code-block:: python

    python3 -m compath_reloaded populate

Alternatively, you can do any of these two steps for a particular database by the following command:

.. code-block:: python

    python3 -m compath_reloaded database_name action

Example:

.. code-block:: python

    python3 -m compath_reloaded kegg download











