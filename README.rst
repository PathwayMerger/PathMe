PathMe |build| |coverage| |docs| |zenodo|
=========================================
PathMe is a Python package aimed to convert KEGG [2]_ [3]_ [4]_, Reactome [5]_ [6]_, and WikiPathways [7]_ [8]_ [9]_ to
Biological Expression Language (BEL).

This project is the continuation of the ComPath web application aimed at exploring, analyzing,
and curating pathway knowledge in a gene-centric view. This different approach involves converting
all the pathways in these resources into BEL as a pivotal integration schema to harmonize entities and relationships in
order across these multiple resources; thus, enabling a more comprehensive evaluation of pathway cross-talks, consensus,
and boundaries. Additionally, PathMe is complemented with the
`PathMe-Viewer <https://github.com/ComPath/PathMe-Viewer>`_, a web application that enables querying, browsing, and
navigating  pathway knowledge assisted by a user-friendly visualization.

Database Versions
-----------------
PathMe currently uses the following versions of the databases:

- **KEGG**: Up-to-date (KEGG does not have tag its releases)
- **Reactome**: 67 Release
- **WikiPathways**: September 2019 Release

Citation
--------
If you use PathMe in your work, please consider citing:

.. [1] Domingo-Fernández, D., *et al.* (2019). `PathMe: Merging and exploring mechanistic pathway knowledge
    <https://doi.org/10.1186/s12859-019-2863-9>`_. *BMC Bioinformatics*, 20:243.

Installation |pypi_version| |python_versions| |pypi_license|
------------------------------------------------------------
``pathme`` can be directly installed from PyPi with pip:

.. code-block:: sh

    $ python3 -m pip install pathme

To use the latest version install directly from GitHub:

.. code-block:: sh

    $ python3 -m pip install git+https://github.com/PathwayMerger/PathMe.git

2. or in editable mode with:

.. code-block:: sh

    $ git clone https://github.com/PathwayMerger/PathMe.git
    $ cd pathme
    $ python3 -m pip install -e .

How to Use
----------
Before using PathMe, make sure you have installed and populated the `Bio2BEL HGNC <https://github.com/bio2bel/hgnc>`_
and `Bio2BEL ChEBI <https://github.com/bio2bel/chebi>`_ databases (Simple run:"python3 -m bio2bel_hgnc populate" and
"python3 -m bio2bel_chebi populate") in your favourite terminal.

Each database has three main commands: ``download``, ``bel``, and ``summarize``:

1. **Download content**

PathMe first requires to download the raw files from the original pathway databases. This can be accomplished by
running the command ('database' can be either KEGG, Reactome, or WikiPathways). E.g., python3 -m pathme kegg download

.. code-block:: sh

    $ python3 -m pathme <database> download

2. **Generate BEL Graphs**

Once the raw files are downloaded, you can run the following to command to generate BELGraphs that will be exported
as Python pickles files for further analysis. Furthermore, the conversion to BEL can be tuned differently for each
database by using specific commands. For example, KEGG parameters are shown when running "python3 -m pathme kegg bel
--help". Finally, please bear in mind that converting the Reactome files take up to 8 hours due to the large amount of
its RDF file.

.. code-block:: sh

    $ python3 -m pathme <database> bel

2. **Summarize**

Summarizes the result of the conversion to BEL.

.. code-block:: sh

    $ python3 -m pathme <database> summarize

Advanced Parameters
-------------------
KEGG Functionalities
~~~~~~~~~~~~~~~~~~~~
The KEGG module of PathMe is able to handle KGML differently depending on the goal. By default, KEGG groups
together the complex of nodes (e.g., gene families) into one node as it is depicted in the KEGG cartoons and
represented in the KGML files. However, this behavior can be modified by adding the parameter `--flatten=True`
in the exporting command. Example:

.. code-block:: bash

    $ python3 -m pathme kegg bel --flatten

References
----------
KEGG
~~~~
PathMe makes use of KEGG KGML files that are downloaded via the KEGG API for academic purposes (please make sure you comply their `Terms and Conditions <https://www.kegg.jp/kegg/rest/>`_).

.. [2] Kanehisa, *et al.* (2017) KEGG: new perspectives on genomes, pathways, diseases and drugs. Nucleic Acids Res. 45,
       D353-D361.
.. [3] Kanehisa, M., *et al.* (2016). KEGG as a reference resource
       for gene and protein annotation. Nucleic Acids Res. 44, D457-D462.
.. [4] Kanehisa, M. and Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 27-30.

Reactome
~~~~~~~~
.. [5] Fabregat, A *et al.* (2016). The Reactome Pathway Knowledgebase. Nucleic Acids Research 44. Database issue:
       D481–D487.
.. [6] Croft, D *et al.* (2014). The Reactome Pathway Knowledgebase. *Nucleic Acids Research* 42.Database issue:
       D472–D477.

WikiPathways
~~~~~~~~~~~~
.. [7] Slenter, D.N.,  *et al.* (2017). WikiPathways: a multifaceted pathway database bridging metabolomics to other omics
       research. *Nucleic Acids Research*, doi.org/10.1093/nar/gkx1064
.. [8] Kutmon, M., *et al.* (2016). WikiPathways: capturing the full diversity of pathway knowledge Nucl. Acids Res., 44,
       D488-D494.
.. [9] Kelder, T., *et al.* (201). WikiPathways: building research communities on biological pathways. Nucleic Acids Res.
       Jan;40(Database issue):D1301-7

.. |build| image:: https://travis-ci.com/PathwayMerger/PathMe.svg?branch=master
    :target: https://travis-ci.com/PathwayMerger/PathMe
    :alt: Build Status

.. |coverage| image:: https://codecov.io/gh/PathwayMerger/PathMe/coverage.svg?branch=master
    :target: https://codecov.io/gh/PathwayMerger/PathMe?branch=master
    :alt: Coverage Status

.. |docs| image:: http://readthedocs.org/projects/pathme/badge/?version=latest
    :target: https://pathme.readthedocs.io/en/latest/
    :alt: Documentation Status

.. |climate| image:: https://codeclimate.com/github/pathwaymerger/pathme/badges/gpa.svg
    :target: https://codeclimate.com/github/pathwaymerger/pathme
    :alt: Code Climate

.. |python_versions| image:: https://img.shields.io/pypi/pyversions/pathme.svg
    :alt: Stable Supported Python Versions

.. |pypi_version| image:: https://img.shields.io/pypi/v/pathme.svg
    :alt: Current version on PyPI

.. |pypi_license| image:: https://img.shields.io/pypi/l/pathme.svg
    :alt: Apache-2.0

.. |zenodo| image:: https://zenodo.org/badge/146161418.svg
    :target: https://zenodo.org/badge/latestdoi/146161418
