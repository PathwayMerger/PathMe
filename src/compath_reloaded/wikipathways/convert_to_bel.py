# -*- coding: utf-8 -*-

"""This module contains the methods to convert a WikiPathways RDF network into a BELGraph."""
import logging
from typing import Dict, List, Tuple

from pybel import BELGraph
from pybel.dsl.edges import activity
from pybel.dsl.nodes import BaseEntity, abundance, bioprocess, complex_abundance, gene, protein, reaction, rna

from compath_reloaded.utils import parse_id_uri

log = logging.getLogger(__name__)


def convert_to_bel(nodes: Dict[str, Dict], interactions: List[Tuple[str, str, Dict]], pathway_info) -> BELGraph:
    graph = BELGraph(
        name=pathway_info['title'],
        version='1.0.0',
        description=pathway_info['description'],
        pathway_id=pathway_info['pathway_id'],
        authors="Sarah Mubeen, Daniel Domingo-Fernández & Josep Marín-Llaó",
        contact='daniel.domingo.fernandez@scai.fraunhofer.de',
    )

    nodes = nodes_to_bel(nodes)

    for interaction in interactions:
        participants = interaction['participants']
        interaction_info = interaction['metadata']

        add_edges(graph, participants, nodes, interaction_info)

    return graph


def nodes_to_bel(nodes: Dict[str, Dict]) -> Dict[str, BaseEntity]:
    """Convert node to Bel"""
    return {
        node_id: node_to_bel(node_att)
        for node_id, node_att in nodes.items()
    }


def node_to_bel(node: Dict) -> BaseEntity:
    """Create a BEL node."""
    node_types = node.pop('rdf_types')
    uri_id = node.pop('uri_id')
    _, _, node['namespace'], _ = parse_id_uri(uri_id)

    if 'Protein' in node_types:
        return protein(**node)

    elif 'Pathway' in node_types:
        return bioprocess(**node)

    elif 'Rna' in node_types:
        return rna(**node)

    elif 'Metabolite' in node_types:
        return abundance(**node)

    elif 'GeneProduct' in node_types:
        return gene(**node)

    elif 'DataNode' in node_types:
        return abundance(**node)

    else:
        log.warning('Unknown %s', node_types)


def add_edges(graph: BELGraph, participants, nodes, att: Dict):
    """Add edges to BELGraph."""
    uri_id = att['uri_id']
    edge_types = att['rdf_types']
    _, _, namespace, interaction_id = parse_id_uri(uri_id)

    members = set()

    if 'Complex' and 'ComplexBinding' in edge_types:
        for source, target in participants:
            members.update({source, target})

        members = {
            nodes[member_id]
            for member_id in members
        }

        complex_abundance(members=members, identifier=interaction_id, namespace=namespace)


    elif 'Conversion' in edge_types:
        reactants = set()
        products = set()

        for source, target in participants:
            reactants.add(target)
            products.add(target)

        reactants = {
            nodes[source_id]
            for source_id in reactants
        }

        products = {
            nodes[product_id]
            for product_id in products
        }

        reaction_node = reaction(reactants=reactants, products=products)
        graph.add_node_from_data(reaction_node)

    else:
        for source, target in participants:
            u = nodes[source]
            v = nodes[target]

            add_simple_edge(graph, u, v, edge_types, uri_id)


def add_simple_edge(graph: BELGraph, u, v, edge_types, uri_id):
    if 'Stimulation' in edge_types:
        graph.add_increases(u, v, citation=uri_id, evidence='', object_modifier=activity())

    elif 'Inhibition' in edge_types:
        graph.add_decreases(u, v, citation=uri_id, evidence='', object_modifier=activity())

    elif 'Catalysis' in edge_types:
        graph.add_increases(u, v, citation=uri_id, evidence='', object_modifier=activity())

    elif 'TranscriptionTranslation' in edge_types:
        graph.add_translation(u, v)

    elif 'DirectedInteraction' in edge_types:
        graph.add_association(u, v, citation=uri_id, evidence='', annotations={
            'EdgeTypes': {
                t: True
                for t in edge_types
            }
        })

    elif 'Interaction' in edge_types:
        pass

    else:
        pass
