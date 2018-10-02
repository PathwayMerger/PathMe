# -*- coding: utf-8 -*-

"""This module contains the methods to convert a WikiPathways RDF network into a BELGraph."""

import logging
from typing import Dict, List, Tuple

from bio2bel_hgnc import Manager
from pybel import BELGraph
from pybel.dsl import abundance, activity, BaseEntity, bioprocess, complex_abundance, gene, protein, reaction, rna

from pathme.constants import HGNC
from pathme.utils import parse_id_uri
from pathme.wikipathways.utils import evaluate_wikipathways_metadata, get_valid_gene_identifier

log = logging.getLogger(__name__)

__all__ = [
    'convert_to_bel',
]


def convert_to_bel(nodes: Dict[str, Dict], complexes: Dict[str, Dict], interactions: List[Tuple[str, str, Dict]],
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

    for interaction in interactions:
        participants = interaction.pop('participants')
        add_edges(graph, participants, nodes, interaction)

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

    if 'identifier' in node.keys():
        identifier = node['identifier']

    else:
        identifier = uri_id

    _, _, namespace, _ = parse_id_uri(uri_id)

    if isinstance(node['name'], set):
        print('{}'.format(node['name']))
        # TODO: print the wikipathways bps that return a set because they are probably wrong.
        name = list(node['name'])[0]

    if 'Protein' in node_types:
        namespace, name, identifier = get_valid_gene_identifier(node, hgnc_manager)
        return protein(namespace=namespace, name=name, identifier=identifier)

    elif 'Rna' in node_types:
        namespace, name, identifier = get_valid_gene_identifier(node, hgnc_manager)
        return rna(namespace=namespace, name=name, identifier=identifier)

    elif 'GeneProduct' in node_types:
        namespace, name, identifier = get_valid_gene_identifier(node, hgnc_manager)
        return gene(namespace=HGNC, name=name, identifier=identifier)

    elif 'Metabolite' in node_types:
        # FIX node[name]
        return abundance(namespace=namespace, name=node['name'], identifier=identifier)

    elif 'Pathway' in node_types:
        return bioprocess(namespace=namespace, name=node['name'], identifier=identifier)


    elif 'DataNode' in node_types:
        return abundance(namespace=namespace, name=node['name'], identifier=identifier)

    else:
        log.warning('Unknown %s', node_types)


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
    }

    _, _, _, identifier = parse_id_uri(complex['uri_id'])

    complex_bel_node = complex_abundance(members=members, identifier=identifier, namespace='wp_complex')
    graph.add_node_from_data(complex_bel_node)

    return complex_bel_node


def add_edges(graph: BELGraph, participants, nodes, att: Dict):
    """Add edges to BELGraph."""
    uri_id = att['uri_id']
    edge_types = att['interaction_types']
    _, _, namespace, interaction_id = parse_id_uri(uri_id)

    if 'Conversion' in edge_types:
        reactants = set()
        products = set()

        for source, target in participants:
            reactants.add(source)
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
        graph.add_association(u, v, citation=uri_id, evidence='', annotations={'EdgeTypes': edge_types})

    elif 'Interaction' in edge_types:
        pass

    else:
        pass
