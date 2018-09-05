# -*- coding: utf-8 -*-

"""This module contains all the constants used in ComPath-reloaded repo."""

import os

from bio2bel.utils import get_connection


def get_data_dir(module_name):
    """Ensures the appropriate ComPath data directory exists for the given module, then returns the file path

    :param str module_name: The name of the module. Ex: 'compath_reloaded'
    :return: The module's data directory
    :rtype: str
    """
    module_name = module_name.lower()
    data_dir = os.path.join(COMPATH_DIR, module_name)
    os.makedirs(data_dir, exist_ok=True)
    return data_dir


MODULE_NAME = 'compath_reloaded'
COMPATH_DIR = os.environ.get('COMPATH_DIRECTORY', os.path.join(os.path.expanduser('~'), '.compath'))
DATA_DIR = get_data_dir(MODULE_NAME)
DEFAULT_CACHE_CONNECTION = get_connection(MODULE_NAME)

KEGG = 'kegg'
REACTOME = 'reactome'
WIKIPATHWAYS = 'wikipathways'

HGNC = 'HGNC'
ENSEMBL = 'ENSEMBL'
UNIPROT = 'UNIPROT'
CHEBI = 'ChEBI'

KEGG_MODIFICATIONS = {
        'phosphorylation':'Ph',
        'glycosylation':'Glyco',
        'ubiquitination':'Ub',
        'meythylation':'Me'
        }
KEGG_CITATION = '10592173'

KEGG_WIKIPATHWAYS_MAPPINGS = 'https://github.com/ComPath/curation/raw/master/mappings/kegg_wikipathways.xlsx'
KEGG_REACTOME_MAPPINGS = 'https://github.com/ComPath/curation/raw/master/mappings/kegg_reactome.xlsx'
WIKIPATHWAYS_REACTOME_MAPPINGS = 'https://github.com/ComPath/curation/raw/master/mappings/wikipathways_reactome.xlsx'

KEGG_KGML_URL = 'http://rest.kegg.jp/get/{}/kgml'
RDF_REACTOME = 'ftp://ftp.ebi.ac.uk/pub/databases/RDF/reactome/r61/reactome-biopax.tar.bz2'
RDF_WIKIPATHWAYS = 'http://data.wikipathways.org/current/rdf/wikipathways-20180810-rdf-wp.zip'
