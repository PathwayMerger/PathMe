PathMe |build|
==============

The primary goal of this package is to convert KEGG, Reactome, and WikiPathways (see References below) to Biological Expression Language (BEL). PathMe is the continuation of the ComPath web application aimed at exploring, analyzing, and curating pathway knowledge in a more simplistic gene-centric view. This different approach involves converting all the pathways to BEL as a pivotal integration schema and evaluating consensus and gaps in pathway knowledge. Additionally, ComPath Reloaded is complemented with PathMe, a web application that enables the exploration of all the pathways from these resources using the mappings curated from ComPath.

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

1. **Download content**

PathMe first requires to download the raw files from the original pathway databases.

.. code-block:: python

    python3 -m pathme download
    
2. **Generate BEL Graphs**

.. code-block:: python

    python3 -m pathme populate

Alternatively, you can do any of these two steps for a particular database by the following command:

.. code-block:: python

    python3 -m pathme database_name action

Example:

.. code-block:: python

    python3 -m pathme kegg download
    
References
----------

KEGG
~~~~
- Kanehisa, Furumichi, M., Tanabe, M., Sato, Y., and Morishima, K.; KEGG: new perspectives on genomes,
  pathways, diseases and drugs. Nucleic Acids Res. 45, D353-D361 (2017).
- Kanehisa, M., Sato, Y., Kawashima, M., Furumichi, M., and Tanabe, M.; KEGG as a reference resource
  for gene and protein annotation. Nucleic Acids Res. 44, D457-D462 (2016).
- Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 27-30 (2000).

Reactome
~~~~~~~~
- Fabregat, Antonio et al. “The Reactome Pathway Knowledgebase.” Nucleic Acids Research 44.Database issue (2016): D481–D487. PMC. Web. 6 Oct. 2017.
- Croft, David et al. “The Reactome Pathway Knowledgebase.” Nucleic Acids Research 42.Database issue (2014): D472–D477. PMC. Web. 6 Oct. 2017.

WikiPathways
~~~~~~~~~~~~
- Slenter, D.N., et al WikiPathways: a multifaceted pathway database bridging metabolomics to other omics research
  Nucleic Acids Research, (2017) doi.org/10.1093/nar/gkx1064
- Kutmon, M., et al. WikiPathways: capturing the full diversity of pathway knowledge Nucl. Acids Res., 44, D488-D494
  (2016) doi:10.1093/nar/gkv1024
- Kelder, T., et al. WikiPathways: building research communities on biological pathways. Nucleic Acids Res. 2012
  Jan;40(Database issue):D1301-7

.. |build| image:: https://travis-ci.org/ComPath/PathMe.svg?branch=master
    :target: https://travis-ci.org/ComPath/PathMe
    :alt: Build Status
