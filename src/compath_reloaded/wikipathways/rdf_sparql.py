# -*- coding: utf-8 -*-

"""This module contains the methods that run SPARQL queries to create the WikiPathways Graphs."""

import os
from collections import defaultdict
from typing import Dict, List, Union, Tuple

import pybel
import tqdm
from pybel import BELGraph
from pybel_tools import summary
from rdflib.namespace import Namespace, RDFS, RDF, DCTERMS, DC

from compath_reloaded.utils import parse_rdf, parse_namespace_uri
from compath_reloaded.wikipathways.convert_to_bel import convert_to_bel
from compath_reloaded.wikipathways.utils import debug_pathway_info, debug_global_statistics, query_result_to_dict

"""SPARQL string queries"""

#: SPARQL prefixes.
PREFIXES = {
    'wp': Namespace('http://vocabularies.wikipathways.org/wp#'),
    'rdfs': RDFS,
    'rdf': RDF,
    'dcterms': DCTERMS,
    'dc': DC
}

#: SPARQL query to get all the subtypes for a specific primary {type} (DataNode or Interaction) in a pathway network.
GET_ENTRIES_SUBTYPES_SPARQL = """
        SELECT DISTINCT ?uri_id ?uri_type
        WHERE {{
           ?pathway a wp:Pathway .
           ?uri_id dcterms:isPartOf ?pathway .
           ?uri_id a wp:{rdf_type} .
           ?uri_id rdf:type ?uri_type .
        }}
        """

#: SPARQL query to get all data nodes in a pathway network with some arguments.
GET_ALL_DATA_NODES_INFO_SPARQL = """
    SELECT DISTINCT ?uri_id ?uri_type ?identifier ?name
        WHERE {
           ?pathway a wp:Pathway .
           ?uri_id dcterms:isPartOf ?pathway .
           ?uri_id a wp:DataNode .
           ?uri_id rdf:type ?uri_type .
           ?uri_id dc:identifier ?uri_id .
           ?uri_id dcterms:identifier ?identifier .
           ?uri_id rdfs:label ?name .
        }
        """

#: SPARQL query to get all directed interactions in a pathway network with source and target.
GET_ALL_DIRECTED_INTERACTIONS_SPARQL = """
        SELECT DISTINCT ?source ?target ?uri_id ?uri_type
        WHERE {
           ?pathway a wp:Pathway .
           ?uri_id dcterms:isPartOf ?pathway .
           ?uri_id a wp:DirectedInteraction .
           ?uri_id rdf:type ?uri_type .
           ?uri_id wp:source ?sourceUri .
           ?uri_id wp:target ?targetUri .
           ?sourceEntry dc:identifier ?sourceUri .
           ?targetEntry dc:identifier ?targetUri .
           ?sourceEntry dcterms:identifier ?source .
           ?targetEntry dcterms:identifier ?target .
        }
        """

#: SPARQL query to get all interactions in a pathway network.
GET_PATHWAY_INFO_SPARQL = """
        SELECT DISTINCT ?title ?pw_identifier ?description ?pathway_id
        WHERE {
       ?pathway_id a wp:Pathway .
       ?pathway_id dc:title ?title .
       ?pathway_id dcterms:description ?description .
       ?pathway_id dcterms:identifier ?pw_identifier .
       }
       """

"""Queries managers"""


def _get_subtypes(query_results) -> Dict[str, set]:
    """Get all subtypes for a primary type (node or interactions).

    :param rdflib.Graph query_results: Results of the query
    :returns: Subtypes as set of values and a pathway id as keys
    """
    entries_types = defaultdict(set)

    for entry in query_results:
        entry_uri_id = str(entry.uri_id)

        _, _, entry_uri_type = parse_namespace_uri(str(entry.uri_type))

        entries_types[entry_uri_id].add(entry_uri_type)

    return dict(entries_types)


def _get_nodes(rdf_graph) -> Dict[str, Dict[str, Dict[str, str]]]:
    """Get all nodes from a RDF pathway network.

    :param rdflib.Graph rdf_graph: RDF graph object
    :returns: Nodes a dictionary with id as keys and their metadata as values
    """
    query_result_nodes = rdf_graph.query(GET_ALL_DATA_NODES_INFO_SPARQL, initNs=PREFIXES)

    entries_types = _get_subtypes(query_result_nodes)
    nodes = query_result_to_dict(query_result_nodes, rdf_types=entries_types)

    return nodes


def _get_interactions(rdf_graph) -> List[Dict[str, Union[set, Dict[str, str]]]]:
    """Get all interactions from a RDF pathway network.

    :param rdflib.Graph rdf_graph: RDF graph object
    :rtype: list
    :returns: Interactions as a list of dictionaries, where the participants are in an entry and the interaction metadata in other
    """
    query_result_directed_interactions = rdf_graph.query(GET_ALL_DIRECTED_INTERACTIONS_SPARQL, initNs=PREFIXES)
    entries_types = _get_subtypes(query_result_directed_interactions)

    directed_interactions = {}

    for interaction in query_result_directed_interactions:
        entry_uri_id = str(interaction.uri_id)
        if entry_uri_id not in directed_interactions.keys():
            directed_interactions[entry_uri_id] = {
                'participants': set(),
                'metadata': {'uri_id': entry_uri_id, 'rdf_types': entries_types[entry_uri_id]}
            }

        directed_interactions[entry_uri_id]['participants'].add((str(interaction.source), str(interaction.target)))

    return list(directed_interactions.values())


def _get_pathway_metadata(rdf_graph) -> Dict[str, str]:
    """Get information from a pathway network.

    :param rdflib.Graph rdf_graph: RDF graph object
    :returns: Metadata of a pathway as a dictionary, if empty 'unknown' will be assigned by default
    """
    return query_result_to_dict(
        rdf_graph.query(GET_PATHWAY_INFO_SPARQL, initNs=PREFIXES),
        attr_empty=['title', 'identifier', 'description', 'pathway_id'],
    )


def rdf_pathway_to_bel(graph) -> BELGraph:
    """Convert RDF graph to BELGraph

    :param graph: RDF graph
    """
    nodes = _get_nodes(graph)
    interactions = _get_interactions(graph)

    pathway_metadata = _get_pathway_metadata(graph)

    return convert_to_bel(nodes, interactions, pathway_metadata)


"""Statistics functions"""


def get_rdf_statistics(rdf_graph, primary_type) -> Dict[str, int]:
    """Get types statistics for a pathway entries type (nodes or interactions) set.

    :param str rdf_graph: primary entries type identifier (ex: DataNode or Interaction)
    :param str primary_type: primary entries type identifier (ex: DataNode or Interaction)
    """
    query_results = rdf_graph.query(GET_ENTRIES_SUBTYPES_SPARQL.format(rdf_type=primary_type),
                                    initNs=PREFIXES)

    entries_types = _get_subtypes(query_results)
    type_statistics = defaultdict(int)

    for entry_types in entries_types.values():

        for entry_type in entry_types:
            type_statistics[entry_type] += 1

        if entry_types == {primary_type}:
            type_statistics['Untyped' + primary_type] += 1

        if primary_type == {'DirectedInteraction', 'Interaction'}:
            type_statistics['UntypedDirectedInteraction'] += 1

    return type_statistics


def get_wikipathways_statistics(resource_files, resource_folder) -> Tuple[
    Dict[str, Dict[str, int]], Dict[str, Dict[str, Dict[str, int]]]]:
    """Load WikiPathways RDF to BELGraph

    :param graph: RDF file path
    """
    global_statistics = defaultdict(lambda: defaultdict(int))
    all_pathways_statistics = {}

    for rdf_file in tqdm.tqdm(resource_files, desc='Parsing WikiPathways'):
        # Parse pathway rdf_file
        pathway_path = os.path.join(resource_folder, rdf_file)

        rdf_graph = parse_rdf(pathway_path)
        bel_graph = rdf_pathway_to_bel(rdf_graph)

        rdf_nodes_statistics = get_rdf_statistics(rdf_graph, 'DataNode')
        rdf_interaction_statistics = get_rdf_statistics(rdf_graph, 'Interaction')

        rdf_total_nodes = rdf_nodes_statistics.pop('DataNode')

        if 'Interaction' in rdf_interaction_statistics.keys():
            rdf_total_interaction = rdf_interaction_statistics.pop('Interaction')

        pathway_statistics = {
            'RDF nodes': rdf_nodes_statistics,
            'RDF interactions': rdf_interaction_statistics,
            'BEL imported nodes': pybel.struct.summary.count_functions(bel_graph),
            'BEL imported edges': summary.edge_summary.count_relations(bel_graph),
            'bel_vs_rdf': {
                'RDF nodes': rdf_total_nodes,
                'RDF interactions': rdf_total_interaction,
                'BEL imported nodes': bel_graph.number_of_nodes(),
                'BEL imported edges': bel_graph.number_of_edges()
            }
        }

        debug_pathway_info(bel_graph, pathway_path, statistics=pathway_statistics)

        all_pathways_statistics[bel_graph.name] = pathway_statistics

        if pathway_statistics:
            for statistics_type, rdf_types in pathway_statistics.items():
                if not rdf_types:
                    continue

                for rdf_type, value in rdf_types.items():
                    global_statistics[statistics_type][rdf_type] += value

    debug_global_statistics(global_statistics)

    return global_statistics, all_pathways_statistics


def pathway_to_bel(file_path):
    """Convert WikiPathways RDF file to BEL.

    :param str file_path: path to the file
    :rtype: pybel.BELGraph
    """
    rdf_graph = parse_rdf(file_path)
    return rdf_pathway_to_bel(rdf_graph)


def wikipathways_to_pybel(resource_files, resource_folder):
    """Load WikiPathways graphs in PyBEL database.

    :param iter[str] resource_files: iterator with file names
    :param str resource_folder: path folder
    """
    for rdf_file in tqdm.tqdm(resource_files, desc='Loading WikiPathways graphs into PyBEL database'):
        # Parse pathway rdf_file and log stats
        pathway_path = os.path.join(resource_folder, rdf_file)
        bel_graph = pathway_to_bel(pathway_path)

        debug_pathway_info(bel_graph, pathway_path)

        pybel.to_database(bel_graph)