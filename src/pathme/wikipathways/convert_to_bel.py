# -*- coding: utf-8 -*-

"""This module contains the methods to convert a WikiPathways RDF network into a BELGraph."""

import logging
from typing import Dict, List, Tuple

from bio2bel_hgnc import Manager
from pathme.constants import HGNC
from pathme.utils import parse_id_uri, check_multiple
from pathme.wikipathways.utils import evaluate_wikipathways_metadata, get_valid_gene_identifier
from pybel import BELGraph
from pybel.dsl import abundance, activity, BaseEntity, bioprocess, complex_abundance, gene, protein, reaction, rna

log = logging.getLogger(__name__)

__all__ = [
    'convert_to_bel',
]


def convert_to_bel(nodes: Dict[str, Dict], complexes: Dict[str, Dict], interactions: Dict[str, Dict],
                   pathway_info, hgnc_manager: Manager) -> BELGraph:
    """Convert  RDF graph info to BEL."""
    graph = BELGraph(
        name=pathway_info['title'],
        version='1.0.0',
        description=evaluate_wikipathways_metadata(pathway_info['description']),
        pathway_id=pathway_info['pathway_id'],
        authors="Sarah Mubeen, Daniel Domingo-Fernández & Josep Marín-Llaó",
        contact='daniel.domingo.fernandez@scai.fraunhofer.de',
    )

    nodes = nodes_to_bel(nodes, hgnc_manager)
    nodes.update(complexes_to_bel(complexes, nodes, graph))

    for interaction_key, interaction in interactions.items():
        participants = interaction['participants']
        add_edges(graph, participants, nodes, interactions, interaction)

    return graph


def nodes_to_bel(nodes: Dict[str, Dict], hgnc_manager: Manager) -> Dict[str, BaseEntity]:
    """Convert node to Bel"""
    return {
        node_id: node_to_bel(node_att, hgnc_manager)
        for node_id, node_att in nodes.items()
    }


def node_to_bel(node: Dict, hgnc_manager: Manager) -> BaseEntity:
    """Create a BEL node."""
    node_types = node['node_types']
    uri_id = node['uri_id']

    # Get identifier from if exists else use uri_id as identifier
    if 'identifier' in node:
        identifier = node['identifier']
    else:
        identifier = uri_id

    identifier = check_multiple(identifier, 'identifier')

    uri_id = check_multiple(uri_id, 'uri_id')
    _, _, namespace, _ = parse_id_uri(uri_id)

    name = check_multiple(node['name'], 'name')

    # Get dictinoary of multiple identifiers
    if 'identifiers' in node:
        node_ids_dict = node['identifiers']
    else:
        node_ids_dict = node

    if 'Protein' in node_types:
        namespace, name, identifier = get_valid_gene_identifier(node_ids_dict, hgnc_manager)
        return protein(namespace=namespace, name=name, identifier=identifier)

    elif 'Rna' in node_types:
        namespace, name, identifier = get_valid_gene_identifier(node_ids_dict, hgnc_manager)
        return rna(namespace=namespace, name=name, identifier=identifier)

    elif 'GeneProduct' in node_types:
        namespace, name, identifier = get_valid_gene_identifier(node_ids_dict, hgnc_manager)
        return gene(namespace=HGNC, name=name, identifier=identifier)

    elif 'Metabolite' in node_types:
        # Parse URI to get namespace
        _, _, namespace, _ = parse_id_uri(uri_id)
        return abundance(namespace=namespace, name=name, identifier=identifier)

    elif 'Pathway' in node_types:
        # Parse URI to get namespace
        _, _, namespace, _ = parse_id_uri(uri_id)
        return bioprocess(namespace=namespace, name=name, identifier=identifier)

    elif 'DataNode' in node_types:
        # Parse URI to get namespace
        _, _, namespace, _ = parse_id_uri(uri_id)
        return abundance(namespace=namespace, name=name, identifier=identifier)

    else:
        log.debug('Unknown %s', node_types)


def complexes_to_bel(complexes: Dict[str, Dict], nodes: Dict[str, BaseEntity], graph: BELGraph) -> Dict[
    str, BaseEntity]:
    """Convert node to Bel"""
    return {
        complex_id: complex_to_bel(complex, nodes, graph)
        for complex_id, complex in complexes.items()
    }


def complex_to_bel(complex, nodes, graph: BELGraph):
    members = {
        nodes[member_id]
        for member_id in complex['participants']
        if member_id in nodes
    }

    _, _, _, identifier = parse_id_uri(complex['uri_id'])

    complex_bel_node = complex_abundance(members=members, identifier=identifier)
    graph.add_node_from_data(complex_bel_node)

    return complex_bel_node

def get_reaction_node(participants, nodes, interactions):
    reactants = set()
    products = set()

    for source, target in participants:
        source = get_node(source, nodes, interactions)
        if source:
            reactants.add(source)

        target = get_node(target, nodes, interactions)
        if target:
            products.add(target)

    return reaction(reactants=reactants, products=products)

def get_node(node, nodes, interactions):
    if node not in nodes:
        if '/Interaction/' in str(node):
            _, _, _, identifier = parse_id_uri(node)

            if identifier in interactions:
                return get_reaction_node(interactions[identifier]['participants'], nodes, interactions)

        log.debug('No valid id for node %s', node)
        return None
    else:
        return nodes[node]


def add_edges(graph: BELGraph, participants, nodes, interactions: Dict, att: Dict):
    """Add edges to BELGraph."""
    uri_id = att['uri_id']
    edge_types = att['interaction_types']
    _, _, namespace, interaction_id = parse_id_uri(uri_id)

    if 'Conversion' in edge_types:
        graph.add_node_from_data(get_reaction_node(participants, nodes, interactions))

    else:
        for source, target in participants:
            u = get_node(source, nodes, interactions)
            v = get_node(target, nodes, interactions)

            if u and v:
                add_simple_edge(graph, u, v, edge_types, uri_id)


def add_simple_edge(graph: BELGraph, u, v, edge_types, uri_id):
    """Add simple edge to graph.

    :param graph: BEL Graph
    :param u: source
    :param v: target
    :param edge_types: edge type dict
    :param uri_id: citation URI
    """
    if 'Stimulation' in edge_types:
        graph.add_increases(u, v, citation=uri_id, evidence='', object_modifier=activity(), annotations={})

    elif 'Inhibition' in edge_types:
        graph.add_decreases(u, v, citation=uri_id, evidence='', object_modifier=activity(), annotations={})

    elif 'Catalysis' in edge_types:
        graph.add_increases(u, v, citation=uri_id, evidence='', object_modifier=activity(), annotations={})

    elif 'DirectedInteraction' in edge_types:
        graph.add_association(u, v, citation=uri_id, evidence='', annotations={'EdgeTypes': edge_types})

    elif 'Interaction' in edge_types:
        log.debug('No interaction subtype for %s', str(uri_id))

    elif 'TranscriptionTranslation' in edge_types:
        graph.add_translation(u, v)

    else:
        log.debug('No handled edge type %s', str(uri_id))
