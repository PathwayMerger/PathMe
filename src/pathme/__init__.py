# -*- coding: utf-8 -*-

"""
The goal of this package is to facilitate the evaluation of pathway knowledge across three of the major pathway databases by harmozing and consolidating different formats.
PathMe does that by converting KEGG, Reactome, and WikiPathways to Biological Expression Language (BEL).
Once the three databases are harmonized into BEL, we can evaluate the consensus and gaps in pathway knowledge.
For that, PathMe is complemented with a web application (`PathMe Viewer <https://github.com/ComPath/PathMe-Viewer>`_)
that enables the exploration of all the pathways from these three resources. PathMe is the follow-up of the ComPath web
application which is aimed at exploring, analyzing, and curating pathway knowledge in a more simplistic, gene-centric
view.

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
"""

import logging

log = logging.getLogger(__name__)

__version__ = '0.0.3'

__title__ = 'pathme'
__description__ = "Harmonizing pathway databases using Biological Expression Language (BEL)"
__url__ = 'https://github.com/ComPath/PathMe'

__author__ = 'Daniel Domingo-Fernández, Sarah Mubeen, and Josep Marín-Llaó'
__email__ = 'daniel.domingo.fernandez@scai.fraunhofer.de'

__license__ = 'MIT License'
__copyright__ = 'Copyright (c) 2018 Daniel Domingo-Fernández, Sarah Mubeen, and Josep Marín-Llaó'
