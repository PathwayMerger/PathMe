# -*- coding: utf-8 -*-

"""PathMe.

PathMe facilitates a systematic comparison across three of the major pathway databases.

The primary goal of this package is to convert KEGG, Reactome, and WikiPathways to Biological Expression Language (BEL).
PathMe is the followw-up of the ComPath web application which is aimed at exploring, analyzing, and curating pathway
knowledge in a more simplistic, gene-centric view. This approach involves converting all pathways to BEL as a pivotal
integration schema and evaluating consensus and gaps in pathway knowledge. Additionally, PathMe is complemented with
a web application (`PathMe Viewer <https://github.com/ComPath/PathMe-Viewer>`_) that enables the exploration of all the
pathways from these three resources.

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

__version__ = '0.0.1'

__title__ = 'pathme'
__description__ = "A systematic comparison across pathway databases"
__url__ = 'https://gitlab.scai.fraunhofer.de/sarah.mubeen/PathMe'

__author__ = 'Sarah Mubeen, Daniel Domingo-Fernández & Josep Marín-Llaó'
__email__ = 'daniel.domingo.fernandez@scai.fraunhofer.de'

__license__ = 'MIT License'
__copyright__ = 'Copyright (c) 2018 Sarah Mubeen, Daniel Domingo-Fernández & Josep Marín-Llaó'
