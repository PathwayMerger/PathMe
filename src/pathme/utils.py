# -*- coding: utf-8 -*-

"""Common utils."""

import logging
import os
import pickle
from collections import defaultdict
from typing import Dict, List, Optional
from urllib.request import urlretrieve

import click
import pandas as pd
import pybel
import rdflib
from pybel import union
from pybel_tools import summary
from pathme.constants import UNKNOWN

from .constants import DATA_DIR

log = logging.getLogger(__name__)


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


def parse_id_uri(uri):
    """Get the components of a given uri (with identifier at the last position).

    :param str uri: URI
    :returns: prefix (ex: http://rdf.wikipathways.org/...)
    :returns: prefix_namespaces: if there are many namespaces, until the penultimate (ex: .../Pathway/WP22_r97775/...)
    :returns: namespace: if there are many namespaces, the last (ex: .../Interaction/)
    :returns: identifier (ex: .../c562c/)
    :rtype: tuple[str,str,str,str]
    """
    # Split the uri str by '/'.
    splitted_uri = uri.split('/')

    # Get the uri components into different variables.
    prefix = '/'.join(splitted_uri[0:3])
    prefix_namespaces = '/'.join(splitted_uri[3:-2])
    namespace = splitted_uri[-2]
    identifier = splitted_uri[-1]

    return prefix, prefix_namespaces, namespace, identifier


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
        if label and entry[label] is not None
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


def query_result_to_dict(entries, **kwargs) -> Dict[str, Dict[str, Dict[str, str]]]:
    """Export to a dictionary a SPARQL query result data structure.

    :param str rdflib.plugins.sparql.processor.SPARQLResult: SPARQL query result data structure, with all the arguments queried for all entries of a certain primary type.
    :returns: entries_dict: Dictionary with all the entries id as keys and the entries arguments as values.
    :rtype: dict
    """
    entries_dict = {}

    for rdf_entry in entries:
        dict_rdf_entry = entry_result_to_dict(rdf_entry, **kwargs)

        if 'identifier' in rdf_entry.labels:
            id_key = str(rdf_entry.identifier)

        elif 'uri_id' in rdf_entry.labels:
            id_key = rdf_entry.uri_id

        else:
            raise Exception

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

    return entries_dict


"""Statistics functions"""


def get_entry_statitics(types_list, primary_type=None, **kwargs):
    """Get types statistics for a pathway entries type (nodes or interactions) set.

    :param str rdf_graph: primary entries type identifier (ex: DataNode or Interaction)
    :param str primary_type: primary entries type identifier (ex: DataNode or Interaction)
    """
    type_statistics = defaultdict(int)

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
    pathways_statistics = defaultdict(list)
    rows = []

    column_types = set()
    column_primary_types_dict = defaultdict(set)

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
