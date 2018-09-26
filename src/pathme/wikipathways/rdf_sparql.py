# -*- coding: utf-8 -*-

"""This module contains the methods that run SPARQL queries to create the WikiPathways Graphs."""

import os
from collections import defaultdict
from typing import Dict, List, Tuple

import pybel
import tqdm
from pybel import BELGraph
from rdflib.namespace import Namespace, RDFS, RDF, DCTERMS, DC

from pathme.utils import get_pathway_statitics, query_result_to_dict, parse_rdf
from pathme.wikipathways.convert_to_bel import convert_to_bel
from pathme.wikipathways.utils import debug_pathway_info

"""SPARQL string queries"""

#: SPARQL prefixes.
PREFIXES = {
    'wp': Namespace('http://vocabularies.wikipathways.org/wp#'),
    'rdfs': RDFS,
    'rdf': RDF,
    'dcterms': DCTERMS,
    'dc': DC,
    'hgnc': Namespace('http://identifiers.org/hgnc.symbol/')
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
GET_ALL_DATA_NODES_SPARQL = """
    SELECT DISTINCT ?uri_id (STRAFTER(STR(?uri_type), str(wp:)) AS ?node_types) ?identifier (STRAFTER(STR(?hgnc_uri), str(hgnc:)) AS ?hgnc_symbol) ?hgnc_uri ?name
    WHERE {
       ?pathway a wp:Pathway .
       ?uri_id dcterms:isPartOf ?pathway .
       ?uri_id a wp:DataNode .
       ?uri_id rdf:type ?uri_type .
       ?uri_id dcterms:identifier ?identifier .
       optional {?uri_id wp:bdbHgncSymbol ?hgnc_uri .}
       ?uri_id rdfs:label ?name .
    }
    """

#: SPARQL query to get all data nodes in a pathway network with some arguments.
GET_ALL_COMPLEXES_SPARQL = """
    SELECT DISTINCT ?uri_id (STRAFTER(STR(?uri_type), str(wp:)) AS ?node_types) ?participants ?name
    WHERE {
       ?pathway a wp:Pathway .
       ?uri_id dcterms:isPartOf ?pathway .
       ?uri_id a wp:Complex .
       ?uri_id rdf:type ?uri_type .
       ?uri_id wp:participants ?participants .
       ?participants dcterms:identifier ?participants .
    }
    """

# TODO: Check interaction complexes.
#: SPARQL query to get all directed interactions in a pathway network with source and target.
GET_ALL_DIRECTED_INTERACTIONS_SPARQL = """
    SELECT DISTINCT ?source ?target ?uri_id (STRAFTER(STR(?uri_type), str(wp:)) AS ?interaction_types)
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
        SELECT DISTINCT ?title ?identifier ?description ?pathway_id
        WHERE {
       ?pathway_id a wp:Pathway .
       ?pathway_id dc:title ?title .
       ?pathway_id dcterms:description ?description .
       ?pathway_id dcterms:identifier ?identifier .
       }
       """

"""Queries managers"""


def _get_pathway_metadata(rdf_graph) -> Dict[str, Dict[str, Dict[str, str]]]:
    """Get information from a pathway network.

    :param rdflib.Graph rdf_graph: RDF graph object
    :returns: Metadata of a pathway as a dictionary, if empty 'unknown' will be assigned by default
    """
    return query_result_to_dict(
        rdf_graph.query(GET_PATHWAY_INFO_SPARQL, initNs=PREFIXES),
        attr_empty=['title', 'identifier', 'description', 'pathway_id'],
        id_dict=False
    )


def _get_nodes(rdf_graph) -> Dict[str, Dict[str, Dict[str, str]]]:
    """Get all nodes from a RDF pathway network.

    :param rdflib.Graph rdf_graph: RDF graph object
    :returns: Nodes dict with nodes ids as keys and their metadata as values
    """
    return query_result_to_dict(
        rdf_graph.query(GET_ALL_DATA_NODES_SPARQL, initNs=PREFIXES)
    )


def _get_complexes(rdf_graph) -> Dict[str, Dict[str, Dict[str, str]]]:
    """Get all complexes from a pathway RDF network.

    :param rdflib.Graph rdf_graph: RDF graph object
    :returns: Nodes dict with nodes ids as keys and their metadata as values
    """
    return query_result_to_dict(
        rdf_graph.query(GET_ALL_COMPLEXES_SPARQL, initNs=PREFIXES)
    )


def _get_interactions(rdf_graph) -> List[dict]:
    """Get all interactions from a RDF pathway network.

    :param rdflib.Graph rdf_graph: RDF graph object
    :returns: Interactions as a list of dictionaries, where the participants are in an entry and the interaction metadata in other
    """
    return list(query_result_to_dict(
        rdf_graph.query(GET_ALL_DIRECTED_INTERACTIONS_SPARQL, initNs=PREFIXES),
        directed_interaction=('source', 'target')
    ).values())


def _get_pathway_components(graph) -> Tuple[
    Dict[str, Dict[str, Dict[str, str]]], Dict[str, Dict[str, Dict[str, str]]], List[dict]]:
    """Get all components in data structures from a RDF pathway network.

    :param graph: RDF graph object
    :returns: Returns at once the retrievals of each component type (nodes, complexes, interactions) functions.
    """
    return _get_nodes(graph), _get_complexes(graph), _get_interactions(graph)


"""Statistics functions"""


def get_wp_statistics(resource_files, resource_folder) -> Tuple[
    Dict[str, Dict[str, int]], Dict[str, Dict[str, Dict[str, int]]]]:
    """Load WikiPathways RDF to BELGraph

    :param iter[str] resource_files: RDF file path
    :param str resource_folder: RDF file path
    """
    global_statistics = defaultdict(lambda: defaultdict(int))
    all_pathways_statistics = {}

    for rdf_file in tqdm.tqdm(resource_files, desc='Parsing WikiPathways'):
        # Parse pathway rdf_file
        pathway_path = os.path.join(resource_folder, rdf_file)
        rdf_graph = parse_rdf(pathway_path)

        pathway_metadata = _get_pathway_metadata(rdf_graph)

        nodes, complexes, interactions = _get_pathway_components(rdf_graph)
        bel_graph = convert_to_bel(nodes, complexes, interactions, pathway_metadata)

        nodes.update(complexes)

        nodes_types = [
            node['node_types']
            for node in nodes.values()
        ]
        edges_types = [
            interaction['interaction_types']
            for interaction in interactions
        ]

        global_statistics, all_pathways_statistics = get_pathway_statitics(
            nodes_types, edges_types, bel_graph, global_statistics=global_statistics,
            all_pathways_statistics=all_pathways_statistics
        )

    return global_statistics, all_pathways_statistics


"""Conversion functions"""


def rdf_wikipathways_to_bel(rdf_graph) -> BELGraph:
    """Convert RDF graph to BELGraph

    :param rdf_graph: RDF graph
    """
    nodes, complexes, interactions = _get_pathway_components(rdf_graph)
    metadata = _get_pathway_metadata(rdf_graph)
    return convert_to_bel(nodes, complexes, interactions, metadata)


def wikipathways_to_bel(file_path):
    """Convert WikiPathways RDF file to BEL.

    :param str file_path: path to the file
    :rtype: pybel.BELGraph
    """
    rdf_graph = parse_rdf(file_path)
    return rdf_wikipathways_to_bel(rdf_graph)


def wikipathways_to_pybel(resource_files, resource_folder):
    """Load WikiPathways graphs in PyBEL database.

    :param iter[str] resource_files: iterator with file names
    :param str resource_folder: path folder
    """
    for rdf_file in tqdm.tqdm(resource_files, desc='Loading WikiPathways graphs into PyBEL database'):
        # Parse pathway rdf_file and log stats
        pathway_path = os.path.join(resource_folder, rdf_file)
        bel_graph = wikipathways_to_bel(pathway_path)

        debug_pathway_info(bel_graph, pathway_path)

        pybel.to_database(bel_graph)
