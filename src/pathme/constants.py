# -*- coding: utf-8 -*-

"""This module contains all the constants used in PathMe repo."""

import os

from bio2bel.utils import get_connection


def get_data_dir():
    """Ensures the appropriate PathMe data directory exists for the given module, then returns the file path.

    :return: The module's data directory
    :rtype: str
    """
    os.makedirs(PATHME_DIR, exist_ok=True)
    return PATHME_DIR


def ensure_pathme_folders():
    """Ensure data folders are created."""
    os.makedirs(KEGG_DIR, exist_ok=True)
    os.makedirs(REACTOME_DIR, exist_ok=True)
    os.makedirs(WIKIPATHWAYS_DIR, exist_ok=True)
    os.makedirs(KEGG_CACHE, exist_ok=True)

    os.makedirs(KEGG_BEL, exist_ok=True)
    os.makedirs(REACTOME_BEL, exist_ok=True)
    os.makedirs(WIKIPATHWAYS_BEL, exist_ok=True)

    os.makedirs(KEGG_FILES, exist_ok=True)
    os.makedirs(REACTOME_FILES, exist_ok=True)
    os.makedirs(WIKIPATHWAYS_FILES, exist_ok=True)


MODULE_NAME = 'pathme'
PATHME_DIR = os.environ.get('PATHME_DIRECTORY', os.path.join(os.path.expanduser('~'), '.pathme'))
DATA_DIR = get_data_dir()
DEFAULT_CACHE_CONNECTION = get_connection(MODULE_NAME)

KEGG = 'kegg'
INTERPRO = 'interpro'
PFAM = 'pfam'

REACTOME = 'reactome'
WIKIPATHWAYS = 'wikipathways'

KEGG_DIR = os.path.join(DATA_DIR, KEGG)
REACTOME_DIR = os.path.join(DATA_DIR, REACTOME)
WIKIPATHWAYS_DIR = os.path.join(DATA_DIR, WIKIPATHWAYS)

KEGG_BEL = os.path.join(KEGG_DIR, 'bel')
REACTOME_BEL = os.path.join(REACTOME_DIR, 'bel')
WIKIPATHWAYS_BEL = os.path.join(WIKIPATHWAYS_DIR, 'bel')

KEGG_FILES = os.path.join(KEGG_DIR, 'xml')
REACTOME_FILES = os.path.join(REACTOME_DIR, 'rdf')
WIKIPATHWAYS_FILES = os.path.join(WIKIPATHWAYS_DIR, 'rdf')

KEGG_CACHE = os.path.join(DATA_DIR, KEGG, 'cache')

ensure_pathme_folders()

UNKNOWN = 'unknown'

KEGG_ID = 'kegg_id'
KEGG_NAME = 'kegg_name'
KEGG_TYPE = 'kegg_type'
HGNC = 'HGNC'
HGNC_SYMBOL = 'HGNC symbol'
ENSEMBL = 'ENSEMBL'
EXPASY = 'EXPASY'
ENTREZ = 'ENTREZ'
UNIPROT = 'UniProt'
CHEBI = 'ChEBI'
CHEBI_NAME = 'ChEBI name'
PUBCHEM = 'PubChem'
WIKIPEDIA = 'WIKIPEDIA'

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
RDF_WIKIPATHWAYS = 'http://data.wikipathways.org/20181010/rdf/wikipathways-20181010-rdf-wp.zip'

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

BEL_STATS_COLUMN_NAMES = {
    'nodes': 'BEL Nodes',
    'edges': 'BEL Edges',
    'Protein': 'BEL Proteins',
    'Gene': 'BEL Genes',
    'Composite': 'BEL Composites',
    'RNA': 'BEL RNA Entities',
    'Complex': 'BEL Complexes',
    'Abundance': 'BEL Compounds',
    'BiologicalProcess': 'BEL Biological Processes',
    'Reaction': 'BEL Reactions',
    'increases': 'BEL Increase Relations',
    'decreases': 'BEL Decrease Relations',
    'association': 'BEL Association Relations',
    'hasComponent': 'BEL Component Edges',
    'hasVariant': 'BEL Variant Edges',
    'hasReactant': 'BEL Reactants Edges',
    'hasProduct': 'BEL Products Edges',
    'translatedTo': 'BEL Translation Edges'
}
