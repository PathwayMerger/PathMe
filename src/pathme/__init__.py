# -*- coding: utf-8 -*-

"""Harmonizing pathway knowledge across three of the major pathway databases by consolidating different formats.

PathMe does that by converting KEGG, Reactome, and WikiPathways to Biological Expression Language (BEL).
Once the three databases are harmonized into BEL, we can evaluate the consensus and gaps in pathway knowledge.
For that, PathMe is complemented with a web application (`PathMe Viewer <https://github.com/ComPath/PathMe-Viewer>`_)
that enables the exploration of all the pathways from these three resources. PathMe is the follow-up of the ComPath web
application which is aimed at exploring, analyzing, and curating pathway knowledge in a more simplistic, gene-centric
view.

It can be used by doing the following:

1. Download content with ``pathme download``
2. Generate BEL Graphs with ``pathme populate``. Alternatively, you can do the
   first two steps for a particular database with ``pathme <database_name> <action>``
"""
