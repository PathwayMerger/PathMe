PathMe |build|
==============

The primary goal of this package is to convert KEGG, Reactome, and WikiPathways to Biological Expression Language (BEL). PathMe is the continuation of the ComPath web application aimed at exploring, analyzing, and curating pathway knowledge in a more simplistic gene-centric view. This different approach involves converting all the pathways to BEL as a pivotal integration schema and evaluating consensus and gaps in pathway knowledge. Additionally, ComPath Reloaded is complemented with PathMe, a web application that enables the exploration of all the pathways from these resources using the mappings curated from ComPath.

Installation
------------
1. ``pathme`` can be installed with the following commmands:

.. code-block:: sh

    python3 -m pip install git+https://github.com/ComPath/PathMe.git@master

2. or in editable mode with:

.. code-block:: sh

    git clone https://github.com/ComPath/PathMe.git

.. code-block:: sh

    cd pathme

.. code-block:: sh

    python3 -m pip install -e .
    
How to use
----------

1. Download content

PathMe first requires to download the raw files from the original pathway databases.

.. code-block:: python

    python3 -m pathme download
    
2. Generate BEL Graphs

.. code-block:: python

    python3 -m pathme populate

Alternatively, you can do any of these two steps for a particular database by the following command:

.. code-block:: python

    python3 -m pathme database_name action

Example:

.. code-block:: python

    python3 -m pathme kegg download



.. |build| image:: https://travis-ci.org/ComPath/PathMe.svg?branch=master
    :target: https://travis-ci.org/ComPath/PathMe
    :alt: Build Status










