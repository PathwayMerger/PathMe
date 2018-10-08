# -*- coding: utf-8 -*-

"""This module has utilities method for parsing, handling wikipathways RDF and data."""

import logging
import os
import zipfile
from typing import Dict, List, Tuple

import networkx as nx
from bio2bel_wikipathways import Manager as WikiPathwaysManager

from ..constants import DATA_DIR, ENSEMBL, ENTREZ, EXPASY, HGNC, KEGG, UNIPROT, WIKIPATHWAYS, WIKIPEDIA, INTERPRO, PFAM
from ..utils import get_files_in_folder, check_multiple

WIKIPATHWAYS_DIR = os.path.join(DATA_DIR, WIKIPATHWAYS)

log = logging.getLogger(__name__)


def evaluate_wikipathways_metadata(metadata):
    """Evaluates metadata in wikipathways and returns the string representation.

    :param metadata:
    :rtype: str
    """
    if isinstance(metadata, set):
        return ','.join(metadata)

    return metadata


def _get_update_alias_symbol(hgnc_manager, original_identifier, original_namespace):
    """Try to get current alias symbol.

    :param bio2bel_hgnc.Manager hgnc_manager: hgnc manager
    :param str original_identifier:
    :param str original_namespace:
    :rtype: tuple
    """
    query_result = hgnc_manager.get_hgnc_from_alias_symbol(original_identifier)

    if not query_result:
        log.debug('No found HGNC Symbol for id %s in (%s)', original_identifier, original_namespace)
        return original_namespace, original_identifier, original_identifier

    return HGNC, query_result.symbol, query_result.identifier


def _validate_query(hgnc_manager, query_result, original_identifier, original_namespace):
    """Process and validate HGNC query.

    :param bio2bel_hgnc.Manager hgnc_manager: hgnc manager
    :param query_result:
    :param str original_identifier:
    :param str original_namespace:
    :rtype: tuple
    """
    # If invalid entry from HGNC, try to find updated symbol
    if not query_result and original_namespace == HGNC:
        return _get_update_alias_symbol(hgnc_manager, original_identifier, HGNC)

    # Invalid entry, proceed with invalid identifier
    if not query_result:
        log.debug('No found HGNC Symbol for id %s in (%s)', original_identifier, original_namespace)
        return original_namespace, original_identifier, original_identifier

    # Multiple entries are returned, for UniProt identifiers
    if isinstance(query_result, list):
        if len(query_result) > 1:
            log.debug('UniProt identifier with multiple HGNC:s %s', query_result)
        query_result = query_result[0]

    # Correct entry, use HGNC identifier
    return HGNC, query_result.symbol, query_result.identifier


def get_valid_gene_identifier(node_ids_dict, hgnc_manager):
    """Return protein/gene identifier for a given RDF node.

    :param dict node_ids_dict: node dictionary
    :param bio2bel_hgnc.Manager hgnc_manager: hgnc manager
    :rtype: tuple[str,str,str]
    :return: namespace, name, identifier
    """
    # Try to get hgnc symbol
    if 'bdb_hgncsymbol' in node_ids_dict:

        hgnc_symbol = check_multiple(node_ids_dict['bdb_hgncsymbol'], 'bdb_hgncsymbol')
        hgnc_entry = hgnc_manager.get_gene_by_hgnc_symbol(hgnc_symbol)

        return _validate_query(hgnc_manager, hgnc_entry, hgnc_symbol, HGNC)

    # Try to get ENTREZ id
    elif 'bdb_ncbigene' in node_ids_dict or 'ncbiprotein' in node_ids_dict['uri_id']:
        if 'bdb_ncbigene' in node_ids_dict:
            entrez_id = check_multiple(node_ids_dict['bdb_ncbigene'], 'bdb_ncbigene')
        elif 'ncbiprotein' in node_ids_dict['uri_id']:
            entrez_id = check_multiple(node_ids_dict['identifier'], 'ncbiprotein')

        hgnc_entry = hgnc_manager.get_gene_by_entrez_id(entrez_id)

        return _validate_query(hgnc_manager, hgnc_entry, entrez_id, ENTREZ)

    # Try to get UniProt id
    elif 'bdb_uniprot' in node_ids_dict:
        uniprot_id = check_multiple(node_ids_dict['bdb_uniprot'], 'bdb_uniprot')
        hgnc_entry = hgnc_manager.get_gene_by_uniprot_id(uniprot_id)

        return _validate_query(hgnc_manager, hgnc_entry, uniprot_id, UNIPROT)

    # Try to get ENSEMBL id
    elif 'bdb_ensembl' in node_ids_dict or 'ena.embl' in node_ids_dict['uri_id']:
        if 'bdb_ensembl' in node_ids_dict:
            ensembl_id = check_multiple(node_ids_dict['bdb_ensembl'], 'bdb_ensembl')

        elif 'ena.embl' in node_ids_dict['uri_id']:
            ensembl_id = check_multiple(node_ids_dict['identifier'], 'bdb_ensembl')

        hgnc_entry = hgnc_manager.get_gene_by_uniprot_id(ensembl_id)

        return _validate_query(hgnc_manager, hgnc_entry, ensembl_id, ENSEMBL)

    elif 'ec-code' in node_ids_dict['uri_id']:
        ec_number = check_multiple(node_ids_dict['name'], 'ec-code')
        return EXPASY, ec_number, ec_number

    # Only wikipathways identifier is given
    elif 'bdb_wikidata' in node_ids_dict:
        name = check_multiple(node_ids_dict['name'], 'wikidata')

        # Find out whether the name is a valid HGNC symbol
        hgnc_entry = hgnc_manager.get_gene_by_hgnc_symbol(name)

        # Correct entry, use HGNC identifier
        if hgnc_entry:
            return HGNC, hgnc_entry.symbol, hgnc_entry.identifier

        log.debug('Adding WikiPathways node %s (%s)', name, WIKIPATHWAYS)
        return WIKIPATHWAYS, name, name

    elif WIKIPEDIA.lower() in node_ids_dict['uri_id']:
        wiki_name = check_multiple(node_ids_dict['identifier'], 'wikipedia_id')
        wiki_id = check_multiple(node_ids_dict['name'], 'wikipedia_name')

        log.debug('Adding Wikipedia node %s (%s)', wiki_name, WIKIPATHWAYS)

        return WIKIPEDIA, wiki_name, wiki_id

    elif KEGG.lower() in node_ids_dict['uri_id']:
        kegg_id = check_multiple(node_ids_dict['identifier'], 'kegg_id')
        kegg_name = check_multiple(node_ids_dict['name'], 'kegg_name')

        log.debug('Adding KEGG node %s ', kegg_id)

        return KEGG, kegg_name, kegg_id

    elif INTERPRO.lower() in node_ids_dict['uri_id']:
        interpro_id = check_multiple(node_ids_dict['identifier'], 'interpro_id')
        interpro_name = check_multiple(node_ids_dict['name'], 'interpro_name')
        log.debug('Adding INTERPRO node %s ', interpro_id)

        return INTERPRO, interpro_name, interpro_id

    elif PFAM.lower() in node_ids_dict['uri_id']:
        pfam_id = check_multiple(node_ids_dict['identifier'], 'pfam_id')
        pfam_name = check_multiple(node_ids_dict['name'], 'pfam_name')
        log.debug('Adding PFAM node %s ', pfam_id)

        return PFAM, pfam_name, pfam_id

    elif 'mirbase.mature' in node_ids_dict['uri_id']:
        mirbase_id = check_multiple(node_ids_dict['identifier'], 'mirbase_id')
        mirbase_name = check_multiple(node_ids_dict['name'], 'mirbase_name')
        log.debug('Adding MIRBASE node %s ', mirbase_id)

        return PFAM, mirbase_name, mirbase_id

    elif 'chembl.compound' in node_ids_dict['uri_id']:
        chembl_id = check_multiple(node_ids_dict['identifier'], 'chembl_id')
        chembl_name = check_multiple(node_ids_dict['name'], 'chembl_name')
        log.debug('Adding MIRBASE node %s ', chembl_id)

        return PFAM, chembl_name, chembl_id

    raise Exception('Unknown identifier for node %s', node_ids_dict)


def merge_two_dicts(dict1, dict2):
    """Merge two dictionaries.

    :param dict dict1:
    :param dict dict2:
    :returns: merged_dict
    :rtype: dict
    """
    merged_dict = dict1.copy()  # start with x's keys and values
    merged_dict.update(dict2)  # modifies z with y's keys and values & returns None
    return merged_dict


def convert_to_nx(nodes: Dict[str, Dict], interactions: List[Tuple[str, str, Dict]],
                  pathway_info: Dict) -> nx.MultiDiGraph:
    """Generate a NetworkX Graph from a network data structure (dict with nodes and edges).

    :param nodes: Node id as keys and Node attributes as values
    :param interactions: list of interactions
    :param pathway_info: pathway info dictionary
    """
    graph = nx.MultiDiGraph(graph_att=pathway_info)

    for node, attributes in nodes.items():
        graph.add_node(node, attr_dict=attributes)

    for subj, obj, interaction in interactions:
        graph.add_edge(subj, obj, attr=interaction)

    return graph


def debug_pathway_info(bel_graph, pathway_path, **kwargs):
    """Debug information about the pathway graph representation.

    :param pybel.BELGraph bel_graph: bel graph
    :param str pathway_path: path of the pathway
    """
    log.debug('Pathway id: {}'.format(os.path.basename(pathway_path)))

    pathway_name = bel_graph.name
    log.debug('Pathway Name: {}'.format(pathway_name))

    bel_nodes = bel_graph.number_of_nodes()
    bel_edges = bel_graph.number_of_edges()

    log.debug('Nodes imported to BEL: %s', bel_nodes)
    log.debug('Edges imported to BEL: %s', format(bel_edges))

    if 'statistics' in kwargs:
        statistics = kwargs.get('statistics')
        log.debug('RDF Nodes statistics: %', format(statistics['RDF nodes']))
        log.debug('RDF Edges statistics: %', format(statistics['RDF interactions']))


def debug_global_statistics(global_statistics):
    """Debug pathway statistics.

    :param dict global_statistics: pathway statistics
    """
    for statistics_type, rdf_types in global_statistics.items():
        log.debug('Total statistics for %s', statistics_type)

        for rdf_type, value in rdf_types.items():
            log.debug('%s: %s', rdf_type, value)


"""Download utilities"""


def get_file_name_from_url(url):
    """Get the last part of an URL.

    :param str url: URL
    :rtype: str
    """
    return url.rsplit('/', 1)[1]


def unzip_file(file_path, export_folder):
    """Unzip file into a destination folder.

    :param str file_path: name of the file
    :param str export_folder: name of the file
    """
    zip_ref = zipfile.ZipFile(file_path, 'r')
    zip_ref.extractall(export_folder)
    zip_ref.close()


def filter_wikipathways_files(files):
    """Filter files that have not 'ttl' extension or not start with 'WP'.

    :param iter[str] files: file names
    :rtype: list
    """
    return [
        file
        for file in files
        if file.endswith('.ttl') and file.startswith('WP')
    ]


def get_wikipathways_files(path, connection=None, only_canonical=True):
    """Get WikiPathways RDF files in folder.

    :param str path: folder path
    :param Optional[str] connection: database connection
    :param Optional[bool] only_canonical: only identifiers present in WP bio2bel db
    :return:
    """
    resource_files = get_files_in_folder(path)

    # Filter files in folder that have no turtle extension or do not start with 'WP'
    resource_files = filter_wikipathways_files(resource_files)

    # Skip files not present in wikipathways bio2bel db -> stuffs from reactome and so on...
    if only_canonical:

        wikipathways_manager = WikiPathwaysManager(connection)

        wikipathways_identifiers = {
            pathway.resource_id
            for pathway in wikipathways_manager.get_all_pathways()
        }

        if not wikipathways_identifiers:
            log.warning(
                'Your WikiPathays Bio2BEL Database is empty! '
                '(please run: python3 -m bio2bel_wikipathways populate)'
            )

        resource_files = [
            file
            for file in resource_files
            if file.split('.')[0] in wikipathways_identifiers
        ]

    return resource_files
