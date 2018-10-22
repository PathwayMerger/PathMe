PathMe |build| |coverage| |docs| |zenodo|
=========================================
PathMe is a Python package aimed to convert KEGG, Reactome, and WikiPathways (see References below) to Biological
Expression Language (BEL). It is the continuation of the ComPath web application aimed at exploring, analyzing,
and curating pathway knowledge in a gene-centric view. This different approach involves converting
all the pathways in these resources into BEL as a pivotal integration schema to harmonize entities and relationships in order across these multiple resources; thus, enabling a more comprehensive evaluation of pathway cross-talks, consensus, and boundaries. Additionally, PathMe is complemented with `PathMe-Viewer <https://github.com/ComPath/PathMe-Viewer>`_ , a web application that enables querying, browsing, and navigating  pathway knowledge assisted by a user-friendly visualization.

Installation |pypi_version| |python_versions| |pypi_license|
------------------------------------------------------------
1. ``pathme`` can be installed with the following commands:

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

Each database has three main commands: "download", "bel", and "summarize":

1. **Download content**

PathMe first requires to download the raw files from the original pathway databases. This can be accomplished by
running the command ('database' can be either KEGG, Reactome, or WikiPathways):

.. code-block:: python

    python3 -m pathme 'database' download

2. **Generate BEL Graphs**

Once the raw files are downloaded, you can run the following to command to generate BELGraphs that will be exported
as Python pickles files for further analysis. Furthermore, the conversion to BEL can be tuned differently for each
database by using specific commands. For example, KEGG parameters are shown when running
`python3 -m pathme kegg bel --help`.

.. code-block:: python

    python3 -m pathme 'database' bel

2. **Summarize**

Summarizes the result of the conversion to BEL.

.. code-block:: python

    python3 -m pathme 'database' summarize

Advance parameters
------------------

KEGG Functionalities
~~~~~~~~~~~~~~~~~~~~

The KEGG module of PathMe is able to handle KGML differently depending on the goal. By default, KEGG groups
together the complex of nodes (e.g., gene families) into one node as it is depicted in the KEGG cartoons and
represented in the KGML files. However, this behavior can be modified by adding the parameter `--flatten=True`
in the exporting command. Example:

.. code-block:: python

    python3 -m pathme kegg bel --flatten=True


References
----------

KEGG
~~~~
PathMe makes use of KEGG KGML files that are downloaded via the KEGG API for academic purposes (see `KEGG Terms and
conditions <https://www.kegg.jp/kegg/rest/>`_.).

- Kanehisa, Furumichi, M., Tanabe, M., Sato, Y., and Morishima, K.; KEGG: new perspectives on genomes,
  pathways, diseases and drugs. Nucleic Acids Res. 45, D353-D361 (2017).
- Kanehisa, M., Sato, Y., Kawashima, M., Furumichi, M., and Tanabe, M.; KEGG as a reference resource
  for gene and protein annotation. Nucleic Acids Res. 44, D457-D462 (2016).
- Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 27-30 (2000).

Reactome
~~~~~~~~
- Fabregat, Antonio et al. “The Reactome Pathway Knowledgebase.” Nucleic Acids Research 44.Database issue (2016):
  D481–D487. PMC. Web. 6 Oct. 2017.
- Croft, David et al. “The Reactome Pathway Knowledgebase.” Nucleic Acids Research 42.Database issue (2014):
  D472–D477. PMC. Web. 6 Oct. 2017.

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

.. |coverage| image:: https://codecov.io/gh/ComPath/PathMe/coverage.svg?branch=master
    :target: https://codecov.io/gh/ComPath/PathMe?branch=master
    :alt: Coverage Status

.. |docs| image:: http://readthedocs.org/projects/pathme/badge/?version=latest
    :target: https://pathme.readthedocs.io/en/latest/
    :alt: Documentation Status

.. |climate| image:: https://codeclimate.com/github/compath/pathme/badges/gpa.svg
    :target: https://codeclimate.com/github/compath/pathme
    :alt: Code Climate

.. |python_versions| image:: https://img.shields.io/pypi/pyversions/pathme.svg
    :alt: Stable Supported Python Versions

.. |pypi_version| image:: https://img.shields.io/pypi/v/pathme.svg
    :alt: Current version on PyPI

.. |pypi_license| image:: https://img.shields.io/pypi/l/pathme.svg
    :alt: Apache-2.0

.. |zenodo| image:: https://zenodo.org/badge/146161418.svg
    :target: https://zenodo.org/badge/latestdoi/146161418
