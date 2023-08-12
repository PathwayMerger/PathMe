# -*- coding: utf-8 -*-

"""This module contains the methods to convert a WikiPathways RDF network into a BELGraph."""

import logging
from typing import Any, Dict, Iterable, Mapping, Optional, Tuple

import pybel
from pybel import BELGraph
from pybel.dsl import BaseEntity, abundance, activity, bioprocess, complex_abundance, gene, protein, rna
from .utils import check_multiple, evaluate_wikipathways_metadata, get_valid_gene_identifier, QueryResult
from ..constants import ACTIVITY_ALLOWED_MODIFIERS, HGNC
from ..utils import add_bel_metadata, parse_id_uri

__all__ = [
    'convert_to_bel',
]

logger = logging.getLogger(__name__)


def _check_empty_complex(complex_dict: Dict[str, Dict], nodes: Dict[str, BaseEntity]):
    members = {
        nodes[member_id]
        for member_id in complex_dict['participants']
        if member_id in nodes
    }
    if not members:
        return False

    return True


def convert_to_bel(
    nodes: QueryResult,
    complexes: QueryResult,
    interactions: QueryResult,
    pathways: QueryResult,
) -> BELGraph:
    """Convert  RDF graph info to BEL."""
    graph = BELGraph(
        name=pathways['title'],
        version='1.0.0',
        description=evaluate_wikipathways_metadata(pathways['description']),
        authors="Sarah Mubeen, Daniel Domingo-Fernández & Josep Marín-Llaó",
        contact='daniel.domingo.fernandez@scai.fraunhofer.de',
    )

    add_bel_metadata(graph)

    pathway_id = graph.graph['pathway_id'] = pathways['pathway_id']

    nodes = {
        node_id: node_to_bel(node, pathway_id)
        for node_id, node in nodes.items()
    }
    nodes.update(complexes_to_bel(complexes, nodes, graph))

    for interaction in interactions.values():
        participants = interaction['participants']
        add_edges(graph, participants, nodes, interactions, interaction)

    return graph


def node_to_bel(node: Dict, pathway_id) -> BaseEntity:
    """Create a BEL node."""
    node_types = node['node_types']
    uri_id = node['uri_id']

    # Get identifier from if exists else use uri_id as identifier
    if 'identifier' in node:
        identifier = node['identifier']
    else:
        identifier = uri_id

    identifier = check_multiple(identifier, 'identifier', pathway_id)

    uri_id = check_multiple(uri_id, 'uri_id', pathway_id)
    _, _, namespace, _ = parse_id_uri(uri_id)

    name = check_multiple(node['name'], 'name', pathway_id)

    # Get dictinoary of multiple identifiers
    if 'identifiers' in node:
        node_ids_dict = node['identifiers']
    else:
        node_ids_dict = node

    if any(node_type in node_types for node_type in ('Protein', 'Rna', 'GeneProduct')):
        namespace, name, identifier = get_valid_gene_identifier(node_ids_dict, pathway_id)
        if 'Protein' in node_types:
            return protein(namespace=namespace.lower(), name=name, identifier=identifier)
        elif 'Rna' in node_types:
            return rna(namespace=namespace.lower(), name=name, identifier=identifier)
        else:  # 'GeneProduct' in node_types
            return gene(namespace=HGNC, name=name, identifier=identifier)

    elif 'Metabolite' in node_types:
        # Parse URI to get namespace
        _, _, namespace, _ = parse_id_uri(uri_id)
        return abundance(namespace=namespace.lower(), name=name, identifier=identifier)

    elif '/wikipathways/WP' in str(uri_id) and {'DataNode'} == node_types:
        # Check the uri_id if is a Pathway
        _, _, namespace, _ = parse_id_uri(uri_id)
        return bioprocess(namespace=namespace.lower(), name=name, identifier=identifier)

    elif 'DataNode' in node_types:
        # Parse URI to get namespace
        _, _, namespace, _ = parse_id_uri(uri_id)
        return abundance(namespace=namespace.lower(), name=name, identifier=identifier)

    else:
        logger.debug('Unknown %s [pathway=%s]', node_types, pathway_id)


def complexes_to_bel(
    complexes: Dict[str, Dict],
    nodes: Dict[str, BaseEntity],
    graph: BELGraph,
) -> Dict[str, BaseEntity]:
    """Convert node to BEL."""
    return {
        complex_id: complex_to_bel(complex_dict, nodes, graph)
        for complex_id, complex_dict in complexes.items()
        if _check_empty_complex(complex_dict, nodes)
    }


def complex_to_bel(complex_dict, nodes, graph: BELGraph):
    """Convert complex abundance to BEL."""
    members = list({
        nodes[member_id]
        for member_id in complex_dict['participants']
        if member_id in nodes
    })

    _, _, _, identifier = parse_id_uri(complex_dict['uri_id'])

    complex_bel_node = complex_abundance(members=members, identifier=identifier)
    graph.add_node_from_data(complex_bel_node)

    return complex_bel_node


def get_reaction_node(
    participants: Iterable[Tuple[str, str]],
    nodes,
    interactions,
) -> pybel.dsl.Reaction:
    reactants = set()
    products = set()

    for source, target in participants:
        source = get_node(source, nodes, interactions)
        if source:
            reactants.add(source)
        else:
            logging.debug(f'Could not find source for reaction: {source}')

        target = get_node(target, nodes, interactions)
        if target:
            products.add(target)
        else:
            logging.debug(f'Could not find target for reaction: {target}')

    return pybel.dsl.Reaction(reactants=reactants, products=products)


def get_node(
    node: str,
    nodes: Mapping[str, BaseEntity],
    interactions: Mapping[str, Any],
) -> Optional[BaseEntity]:
    if node in nodes:
        return nodes[node]

    if '/Interaction/' in str(node):
        _, _, _, identifier = parse_id_uri(node)
        if identifier in interactions:
            return get_reaction_node(interactions[identifier]['participants'], nodes, interactions)

    logger.debug('No valid id for node %s', node)


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
            if u is None:
                logger.debug(f'Source is none: {source}')
            if v is None:
                logger.debug(f'Target is none: {target}')


def add_simple_edge(graph: BELGraph, u: BaseEntity, v: BaseEntity, edge_types, uri_id):
    """Add simple edge to graph.

    :param graph: BEL Graph
    :param u: source
    :param v: target
    :param edge_types: edge type dict
    :param uri_id: citation URI
    """
    if 'Stimulation' in edge_types:
        graph.add_increases(
            u, v,
            citation=uri_id, evidence='Extracted from WikiPathways',
            object_modifier=activity() if isinstance(v, ACTIVITY_ALLOWED_MODIFIERS) else None,
            annotations={},
        )

    elif 'Inhibition' in edge_types:
        graph.add_decreases(
            u, v,
            citation=uri_id, evidence='Extracted from WikiPathways',
            object_modifier=activity() if isinstance(v, ACTIVITY_ALLOWED_MODIFIERS) else None,
            annotations={},
        )

    elif 'Catalysis' in edge_types:
        graph.add_increases(
            u, v,
            citation=uri_id, evidence='Extracted from WikiPathways',
            object_modifier=activity() if isinstance(v, ACTIVITY_ALLOWED_MODIFIERS) else None,
            annotations={},
        )

    elif 'DirectedInteraction' in edge_types:
        graph.add_regulates(
            u, v,
            citation=uri_id,
            evidence='Extracted from WikiPathways',
            annotations={
                'EdgeTypes': edge_types,
            },
        )

    elif 'Interaction' in edge_types:
        logger.debug('No interaction subtype for %s', uri_id)

    elif 'TranscriptionTranslation' in edge_types:
        graph.add_translation(u, v)

    else:
        logger.debug('No handled edge type %s', uri_id)
