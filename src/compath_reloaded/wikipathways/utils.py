# -*- coding: utf-8 -*-

"""This module has utilities method for parsing, handling wikipathways RDF and data."""

import json
import logging
import os
import zipfile
from typing import Dict, List, Tuple

import networkx as nx
import rdflib
from bio2bel_wikipathways import Manager as WikiPathwaysManager

from ..constants import DATA_DIR, WIKIPATHWAYS
from ..utils import get_files_in_folder

WIKIPATHWAYS_DIR = os.path.join(DATA_DIR, WIKIPATHWAYS)

log = logging.getLogger(__name__)


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


def query_result_to_dict(entries, **kwargs):
    """Export to a dictionary a SPARQL query result data structure.

    :param str rdflib.plugins.sparql.processor.SPARQLResult: SPARQL query result data structure, with all the arguments queried for all entries of a certain primary type.
    :returns: entries_dict: Dictionary with all the entries id as keys and the entries arguments as values.
    :rtype: dict
    """
    entries_dict = {}
    exclude_set = {'uri_type'}

    if not entries and 'attr_empty' in kwargs:
        attr_empty = kwargs.get('attr_empty')
        entries_dict = {
            str(label): 'unknown'
            for label in attr_empty
        }

    for entry in entries:
        attributes_dict = {
            str(label): str(entry[label])
            for label in entry.labels
            if str(label) not in exclude_set and entry[label] is not None
        }

        if 'rdf_types' in kwargs:
            attributes_dict['rdf_types'] = kwargs.get('rdf_types')[str(entry.uri_id)]

        if 'identifier' in entry.labels:
            entries_dict[str(entry.identifier)] = attributes_dict
        else:
            return attributes_dict

    return entries_dict


def convert_json(graph: rdflib.Graph):
    """Convert from rdflib importated graph object to python data structure (list of dicts of each entry).

    :param rdflib.graph graph: graph object
    :rtype: list[dict]
    """
    serialized_json = graph.serialize(format='json-ld', indent=4)
    json_wp_pathway = json.loads(serialized_json.decode("utf-8"))

    return json_wp_pathway


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
    :param pathway_path str: path of the pathway
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
    """Debug pathway statistics

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
            log.warning('Your WikiPathays Bio2BEL Database is empty!')

        resource_files = [
            file
            for file in resource_files
            if file.split('.')[0] in wikipathways_identifiers
        ]

    return resource_files
