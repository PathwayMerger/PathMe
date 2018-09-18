# -*- coding: utf-8 -*-

"""This module contains all the constants used in PathMe repo."""

import os

from bio2bel.utils import get_connection


def get_data_dir(module_name):
    """Ensures the appropriate data directory exists for the given module, then returns the file path.

    :param str module_name: The name of the module. Ex: 'pathme'
    :return: The module's data directory
    :rtype: str
    """
    module_name = module_name.lower()
    data_dir = os.path.join(PATHME_DIR, module_name)
    os.makedirs(data_dir, exist_ok=True)
    return data_dir


MODULE_NAME = 'pathme'
PATHME_DIR = os.environ.get('PATHME_DIRECTORY', os.path.join(os.path.expanduser('~'), '.pathme'))
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
    'phosphorylation': 'Ph',
    'glycosylation': 'Glyco',
    'ubiquitination': 'Ub',
    'methylation': 'Me'
}
KEGG_CITATION = '10592173'

KEGG_WIKIPATHWAYS_MAPPINGS = 'https://github.com/ComPath/curation/raw/master/mappings/kegg_wikipathways.xlsx'
KEGG_REACTOME_MAPPINGS = 'https://github.com/ComPath/curation/raw/master/mappings/kegg_reactome.xlsx'
WIKIPATHWAYS_REACTOME_MAPPINGS = 'https://github.com/ComPath/curation/raw/master/mappings/wikipathways_reactome.xlsx'

KEGG_KGML_URL = 'http://rest.kegg.jp/get/{}/kgml'
RDF_REACTOME = 'ftp://ftp.ebi.ac.uk/pub/databases/RDF/reactome/r61/reactome-biopax.tar.bz2'
RDF_WIKIPATHWAYS = 'http://data.wikipathways.org/current/rdf/wikipathways-20180910-rdf-wp.zip'

KEGG_STATS_COLUMN_NAMES = {
    'nodes': 'BEL Nodes',
    'entities': 'XML Entities',
    'edges': 'BEL Edges',
    'interactions': 'XML Interactions',
    'Protein': 'BEL Proteins',
    'gene': 'XML Genes',
    'Composite': 'BEL Composites',
    'RNA': 'BEL RNA Entities',
    'Complex': 'BEL Complexes',
    'group': 'XML Complexes',
    'Abundance': 'BEL Compounds',
    'compound entity': 'XML Compounds',
    'BiologicalProcess': 'BEL Biological Processes',
    'map': 'XML Biological Processes',
    'Reaction': 'BEL Reactions',
    'ortholog': 'XML Orthologs',
    'increases': 'BEL Increase Relations',
    'activation': 'XML Activation Relations',
    'expression': 'XML Expression Relations',
    'reversible': 'XML Reversible Reactions',
    'irreversible': 'XML Irreversible Reactions',
    'phosphorylation': 'XML Phosphorylation Relations',
    'glycosylation': 'XML Glycosylation Relations',
    'ubiquitination': 'XML Ubiquitination Relations',
    'methylation': 'XML Methylation Relations',
    'decreases': 'BEL Decrease Relations',
    'inhibition': 'XML Inhibition Relations',
    'repression': 'XML Repression Relations',
    'dephosphorylation': 'XML Dephosphorylation Relations',
    'association': 'BEL Association Relations',
    'binding/association': 'XML Binding/Association Relations',
    'indirect effect': 'XML Association Relations',
    'hasComponent': 'BEL Component Edges',
    'hasVariant': 'BEL Variant Edges',
    'hasReactant': 'BEL Reactants Edges',
    'hasProduct': 'BEL Products Edges',
    'compound': 'XML Compound Relations',
    'dissociation': 'XML Dissociation Relations',
    'hidden compound': 'XML Hidden Compound Relations',
    'missing interaction': 'XML Missing Interaction Relations',
    'state change': 'XML State Change Relations',
    'brite': 'XML Brite Hierarchy'
}
