# -*- coding: utf-8 -*-

"""Common utils."""

import collections
import itertools as itt
import logging
import os
import pickle
from typing import Dict, List, Optional, Set
from urllib.parse import urlparse
from urllib.request import urlretrieve

import click
import pandas as pd
import rdflib

from pathme.constants import UNKNOWN, BEL_STATS_COLUMN_NAMES

import pybel
from pybel import from_pickle
from pybel import union
from pybel.struct.summary import count_functions, edge_summary
from pybel_tools import summary

log = logging.getLogger(__name__)


class CallCounted:
    """Decorator to determine number of calls for a method."""

    def __init__(self, method):
        self.method = method
        self.counter = 0

    def __call__(self, *args, **kwargs):
        self.counter += 1
        return self.method(*args, **kwargs)


def get_files_in_folder(path: str) -> List[str]:
    """Return the files in a given folder.

    :param path: folder path
    :return: file names in folder
    """
    return [
        file
        for file in os.listdir(path)
        if os.path.isfile(os.path.join(path, file))
    ]


def check_multiple(element, element_name):
    """Check whether element is iterable.

    :param element: variable to check
    :param element_name: name to print
    :return:
    """
    if isinstance(element, set) or isinstance(element, list):
        log.warning('Multiple {}: {}'.format(element_name, element))
        # TODO: print the wikipathways bps that return a set because they are probably wrong.
        if len(element) != 0:
            return list(element)[0]
        else:
            log.warning('Empty list/set %s', element)

    return element


def parse_id_uri(uri):
    """Get the components of a given uri (with identifier at the last position).

    :param str uri: URI
    :returns: prefix (ex: http://rdf.wikipathways.org/...)
    :returns: prefix_namespaces: if there are many namespaces, until the penultimate (ex: .../Pathway/WP22_r97775/...)
    :returns: namespace: if there are many namespaces, the last (ex: .../Interaction/)
    :returns: identifier (ex: .../c562c/)
    :rtype: tuple[str,str,str,str]
    """
    parsed_url = urlparse(uri)
    uri_suffix = parsed_url.path.split('/')

    # Returns
    # Domain (rdf.wikipathways.org),
    # Prefix namespace (pathway/WP2118_r97625/WP),
    # namespace (Interaction),
    # identifier (id61b0d9c7) in the given example ->
    # (http://rdf.wikipathways.org/Pathway/WP2118_r97625/WP/Interaction/id61b0d9c7)
    return parsed_url.netloc, '/'.join(uri_suffix[0:-2]), uri_suffix[-2], uri_suffix[-1]


def parse_namespace_uri(uri):
    """Get the prefix and namespace of a given URI (without identifier, only with a namspace at last position).

    :param str uri: URI
    :returns: prefix (ex: http://purl.org/dc/terms/...)
    :returns: namespace (ex: .../isPartOf)
    :rtype: tuple[str,str]
    """

    # Split the uri str by '/'.
    splited_uri = uri.split('/')

    # Get the uri prefix and namespace into different variables.
    prefix = '/'.join(splited_uri[0:-1])
    namespace = splited_uri[-1]
    vocabulary = namespace.split('#')[-1]

    return prefix, namespace, vocabulary


def parse_rdf(path: str, format: Optional[str] = None) -> rdflib.Graph:
    """Import a queried pathway into a rdflib Graph object.

    :param path: RDF file path
    :param format: RDF file format, default is turtle
    """
    if format is None:
        format = 'ttl'

    pickle_path = f'{path}.pickle'
    if os.path.exists(pickle_path):
        with open(pickle_path, 'rb') as file:
            return pickle.load(file)

    if not os.path.exists(path):
        raise FileNotFoundError(
            'You have still not downloaded the database file.'
            'Please run "python3 -m pathme "database" download"'
        )

    graph = rdflib.Graph()
    graph.parse(path, format=format)

    with open(pickle_path, 'wb') as file:
        pickle.dump(graph, file)

    return graph


def entry_result_to_dict(entry, **kwargs):
    """Export to a dictionary a SPARQL query result data structure.

    :param str rdflib.plugins.sparql.processor.SPARQLResult: SPARQL query result data structure, with all the arguments queried for all entries of a certain primary type.
    :returns: entries_dict: Dictionary with all the entries id as keys and the entries arguments as values.
    :rtype: dict
    """
    attributes_dict = {
        str(label): str(entry[label])
        for label in entry.labels
        if label is not None and entry[label] is not None
    }

    if 'directed_interaction' in kwargs:
        directed_interaction = kwargs.get('directed_interaction')

        if directed_interaction[0] and directed_interaction[1] in attributes_dict.keys():
            attributes_dict['participants'] = (
                attributes_dict.pop(directed_interaction[0]), attributes_dict.pop(directed_interaction[1]))

    if 'attr_empty' in kwargs:
        attr_empty = kwargs.get('attr_empty')
        empty_dict = {
            str(attr): UNKNOWN
            for attr in attr_empty
            if attr not in attributes_dict
        }
        attributes_dict.update(empty_dict)

    return attributes_dict


def entries_dict_ids_argument(entries_dict):
    entries_dict_ids = collections.defaultdict(dict)
    for entry_id, entry_att in entries_dict.items():
        entry_identifiers = {}

        for label, value in entry_att.items():
            if 'bdb' in label:
                entry_identifiers[label] = value
            else:
                entries_dict_ids[entry_id][label] = value
        entries_dict_ids[entry_id]['identifiers'] = entry_identifiers

    return entries_dict


def query_result_to_dict(entries, **kwargs) -> Dict[str, Dict[str, Dict[str, str]]]:
    """Export to a dictionary a SPARQL query result data structure.

    :param str rdflib.plugins.sparql.processor.SPARQLResult: SPARQL query result data structure, with all the arguments queried for all entries of a certain primary type.
    :returns: entries_dict: Dictionary with all the entries id as keys and the entries arguments as values.
    :rtype: dict
    """
    entries_dict = {}

    for rdf_entry in entries:
        if 'identifier' in rdf_entry.labels and rdf_entry.identifier is not None:
            id_key = str(rdf_entry.identifier)

        elif 'uri_id' in rdf_entry.labels and rdf_entry.uri_id is not None:
            id_key = rdf_entry.uri_id

        else:
            raise Exception(f'invalid data: {rdf_entry}')

        dict_rdf_entry = entry_result_to_dict(rdf_entry, **kwargs)

        if id_key not in entries_dict:
            entries_dict[id_key] = dict_rdf_entry
            if 'participants' in dict_rdf_entry:
                entries_dict[id_key]['participants'] = {dict_rdf_entry['participants']}

        else:
            entry_attributes = entries_dict[id_key]

            for label, new_value in dict_rdf_entry.items():
                if label in entry_attributes.keys():
                    value = entry_attributes[label]
                    if isinstance(value, set):
                        entries_dict[id_key][label].add(new_value)
                    elif entries_dict[id_key][label] == UNKNOWN:
                        entries_dict[id_key][label] = new_value
                    elif value != new_value:
                        entries_dict[id_key][label] = {value, new_value}
                else:
                    entries_dict[id_key][label] = new_value

    if len(entries_dict) == 1 and kwargs.get('id_dict') == False:
        return list(entries_dict.values())[0]

    elif not entries and 'attr_empty' in kwargs:

        attr_empty = kwargs.get('attr_empty')
        return {
            str(attr): UNKNOWN
            for attr in attr_empty
        }

    if kwargs.get('ids_argument') == True:
        return entries_dict_ids_argument(entries_dict)

    else:
        return entries_dict


"""Statistics functions"""


def get_entry_statitics(types_list, primary_type=None, **kwargs):
    """Get types statistics for a pathway entries type (nodes or interactions) set.

    :param str rdf_graph: primary entries type identifier (ex: DataNode or Interaction)
    :param str primary_type: primary entries type identifier (ex: DataNode or Interaction)
    """
    type_statistics = collections.defaultdict(int)

    for entry_types in types_list:
        if isinstance(entry_types, set):
            for entry_type in entry_types:
                if entry_type != primary_type:
                    type_statistics[entry_type] += 1
        else:
            type_statistics[entry_types] += 1

        if 'primary_type' in kwargs:
            if entry_types == {primary_type}:
                type_statistics['Untyped' + primary_type] += 1

            if primary_type == {'DirectedInteraction', 'Interaction'}:
                type_statistics['UntypedDirectedInteraction'] += 1

    return type_statistics, len(types_list)


def get_pathway_statitics(nodes_types, edges_types, bel_graph, **kwargs):
    rdf_nodes_statistics, rdf_total_nodes = get_entry_statitics(nodes_types)
    rdf_edges_statistics, rdf_total_edges = get_entry_statitics(edges_types)

    pathway_statistics = {
        'RDF nodes': rdf_nodes_statistics,
        'RDF interactions': rdf_edges_statistics,
        'BEL imported nodes': pybel.struct.summary.count_functions(bel_graph),
        'BEL imported edges': summary.edge_summary.count_relations(bel_graph),
        'bel_vs_rdf': {
            'RDF nodes': rdf_total_nodes,
            'RDF interactions': rdf_total_edges,
            'BEL imported nodes': bel_graph.number_of_nodes(),
            'BEL imported edges': bel_graph.number_of_edges()
        }
    }

    if 'global_statistics' in kwargs and pathway_statistics:
        global_statistics = kwargs.get('global_statistics')
        for statistics_type, rdf_types in pathway_statistics.items():
            if not rdf_types:
                continue

            for rdf_type, value in rdf_types.items():
                global_statistics[statistics_type][rdf_type] += value

        if 'all_pathways_statistics' in kwargs and pathway_statistics:
            all_pathways_statistics = kwargs.get('all_pathways_statistics')
            all_pathways_statistics[bel_graph.name] = pathway_statistics
            return global_statistics, all_pathways_statistics

        return global_statistics, pathway_statistics

    return pathway_statistics


def statistics_to_df(all_pathways_statistics):
    """Build a data frame with graph statistics.

    :param dict all_pathways_statistics: pathway statistics
    :rtype: pandas.DataFrame
    """
    pathways_statistics = collections.defaultdict(list)
    rows = []

    column_types = set()
    column_primary_types_dict = collections.defaultdict(set)

    # Get pathway type statistics
    for pathway_name, statistics_primary_type_dict in all_pathways_statistics.items():
        for statistic_primary_type, statistic_dict in statistics_primary_type_dict.items():
            column_types.update(set(statistic_dict.keys()))
            column_primary_types_dict[statistic_primary_type].update(set(statistic_dict.keys()))

    for pathway_name, statistics_primary_type_dict in all_pathways_statistics.items():
        rows.append(pathway_name)
        for column_type in column_types:
            for statistic_primary_type, statistic_dict in statistics_primary_type_dict.items():
                if column_type in statistic_dict.keys():
                    pathways_statistics['"' + column_type + '" ' + statistic_primary_type].append(
                        str(statistic_dict[column_type]))
                elif column_type in column_primary_types_dict[statistic_primary_type]:
                    pathways_statistics['"' + column_type + '" ' + statistic_primary_type].append('')

    df = pd.DataFrame(data=pathways_statistics, index=rows)

    return df


def get_bel_types(path):
    """Get BEL node and edge type statistics.

    :param str path: path to pickle
    :return: count of all nodes and edges in a BEL graph
    :rtype: dict
    """
    bel_stats = {}

    bel_graph = from_pickle(path)

    bel_stats['nodes'] = bel_graph.number_of_nodes()
    bel_stats['edges'] = bel_graph.number_of_edges()

    # Get count of all BEL function types
    bel_functions_dict = count_functions(bel_graph)
    bel_stats.update(bel_functions_dict)

    # Get count of all BEL edge types
    bel_edges_dict = edge_summary.count_relations(bel_graph)
    bel_stats.update(bel_edges_dict)

    return bel_stats


def get_bel_stats(resource_folder):
    """Get all BEL node and edge type statistics.

    :param str resource_folder: path to BEL pickles folder
    :return: count of all nodes and edges in all BEL graphs from one resource
    :rtype: dict
    """
    df = pd.DataFrame()

    pickles = get_files_in_folder(resource_folder)

    for file in pickles:
        pathway_names = []
        pathway_names.append(file.strip('.pickle'))

        bel_statistics_dict = get_bel_types(os.path.join(resource_folder, file))

        all_bel_statistics = {
            BEL_STATS_COLUMN_NAMES[key]: value
            for key, value in bel_statistics_dict.items()
        }

        # Add pathway statistic rows to DataFrame
        pathway_data = pd.DataFrame(
            all_bel_statistics,
            index=pathway_names,
            columns=BEL_STATS_COLUMN_NAMES.values(),
            dtype=int
        )

        df = df.append(pathway_data.fillna(0).astype(int))

    return df


def get_genes_from_pickles(resource_folder: str, files: List[str], manager) -> Dict[str,set]:
    """Get BEL graph gene set for all pathways in resource.

    :param str resource_folder: path to resource folder
    :param list files: list of BEL graph pickles
    :param bio2bel Manager manager: Manager
    :return: BEL graph gene sets for each pathway in resource
    :rtype: dict[str,set]
    """
    pathway_genes_dict = {}

    for file_name in files:
        graph = from_pickle(os.path.join(resource_folder, file_name))

        # Get gene set for pathway
        gene_set = get_genes_in_graph(graph)
        file_name = file_name.strip('.pickle')
        file_name = manager.get_pathway_by_id(file_name)
        pathway_genes_dict[str(file_name)] = gene_set

    return pathway_genes_dict


def get_kegg_genes_from_pickles(resource_folder, files: List[str], manager) ->  Dict[str,Set]:
    """Get BEL graph gene set for all KEGG pathways.

    :param str resource_folder: path to resource folder
    :param list files: list of BEL graph pickles
    :param bio2bel Manager manager: Manager
    :return: BEL graph gene sets for each KEGG pathway
    :rtype: dict[str,set]
    """
    pathway_genes_dict = {}

    for file_name in files:

        # Flattened graphs considered for gene sets
        if file_name.endswith('_flatten.pickle'):
            graph = from_pickle(os.path.join(resource_folder, file_name))

            # Get gene set for pathway
            gene_set = get_genes_in_graph(graph)
            file_name = file_name.strip('_flatten.pickle')
            file_name = 'path:' + file_name
            file_name = manager.get_pathway_by_id(file_name)

            pathway_genes_dict[str(file_name)] = gene_set

    return pathway_genes_dict


def get_genes_in_graph(graph: pybel.BELGraph) -> Set[str]:
    """Get BEL graph gene set for a pathway.

    :param pybel.BELGraph graph: BEL Graph
    :return: BEL graph gene set
    :rtype: set
    """
    gene_set = set()

    for node, data in graph.nodes(data=True):
        if node.function in {'Protein', 'RNA', 'Gene'} and node.namespace == 'HGNC':
            gene_set.add(node.name)

    return gene_set


def jaccard_similarity(database_gene_set, bel_genes_set):
    """Get Jaccard similarity for gene sets in database and BEL graphs.

    :param dict database_gene_set: gene sets for each pathway in database
    :param dict bel_genes_set: gene sets for each BEL graph
    :return: similarity index
    :rtype: int
    """
    jaccard_similarities = []
    count = 0
    count_no_similarity = 0

    for (database_key, database_value), (bel_key, bel_value) in itt.product(database_gene_set.items(),
                                                                            bel_genes_set.items()):

        if database_key != bel_key:
            continue

        intersection = len(set.intersection(database_value, bel_value))
        union = len(database_value.union(bel_value))
        jaccard_index = intersection / union
        jaccard_similarities.append(jaccard_index)

        if jaccard_index == 1.0:
            count += 1

        elif jaccard_index == 0.0:
            count_no_similarity += 1

    print('Jaccard index for gene sets in database vs gene sets in BEL:')
    print('{} of {} gene sets in the database and BEL graphs have a similarity of 100%.'.format(count, len(
        jaccard_similarities)))
    print('{} of {} gene sets in the database and BEL graphs have a similarity of 0%.'.format(count_no_similarity, len(
        jaccard_similarities)))

    return jaccard_similarities


"""Downloader"""


def make_downloader(url, path, export_path, decompress_file):
    """Make a function that downloads the data for you, or uses a cached version at the given path.

    :param str url: The URL of some data
    :param str export_path: folder where decompressed file will be exported
    :param method decompress_file: method to decompress file
    :return: A function that downloads the data and returns the path of the data
    :rtype: (bool -> str)
    """
    log.info(url)

    def download_data(force_download=False):
        """Download the data.

        :param bool force_download: If true, overwrites a previously cached file
        :rtype: str
        """
        if os.path.exists(path) and not force_download:
            log.info('using cached data at %s', path)
        else:
            log.info('downloading %s to %s', url, path)
            urlretrieve(url, path)

        return path

    data = download_data()

    log.info('unzipping file %s, da')
    decompress_file(data, export_path)


def summarize_helper(graphs):
    """Print in console summary of graphs.

    :param iter[graphs] graphs: BEL Graphs
    """
    click.echo('joining graphs')
    graph = union(graphs)

    click.echo('generating summary')
    summary_str = graph.summary_str()

    click.echo(summary_str)
