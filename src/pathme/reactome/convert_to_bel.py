# -*- coding: utf-8 -*-

"""This module contains the methods to convert a Reactome RDF network into a BELGraph."""

import logging
from typing import Dict, List, Tuple

from bio2bel_chebi import Manager as ChebiManager
from bio2bel_hgnc import Manager as HgncManager
from pybel import BELGraph
from pybel.dsl import (
    abundance,
    activity,
    composite_abundance,
    complex_abundance,
    gene,
    rna,
    protein,
    reaction,
    bioprocess,
    BaseEntity,
    NamedComplexAbundance
)

from pathme.constants import ACTIVITY_ALLOWED_MODIFIERS, UNKNOWN, REACTOME_CITATION
from pathme.reactome.utils import get_valid_node_parameters, process_multiple_proteins
from pathme.utils import add_bel_metadata, parse_id_uri

log = logging.getLogger(__name__)

__all__ = [
    'convert_to_bel',
]


def convert_to_bel(nodes: Dict[str, Dict], interactions: List[Tuple[str, str, Dict]], pathway_info: Dict,
                   hgnc_manager: HgncManager, chebi_manager: ChebiManager) -> BELGraph:
    """Convert RDF graph dictionary into BEL graph."""
    uri_id = pathway_info['uri_reactome_id']

    if uri_id != UNKNOWN:
        _, _, namespace, identifier = parse_id_uri(uri_id)
    else:
        identifier = UNKNOWN

    description = pathway_info['comment']
    if isinstance(description, (set, list)):
        description = '\n'.join(description)

    """Convert graph-like dictionaries to BELGraph."""
    graph = BELGraph(
        name=pathway_info['display_name'],
        version='1.0.0',
        description=description,
        authors="Josep Marín-Llaó, Daniel Domingo-Fernández & Sarah Mubeen",
        contact='daniel.domingo.fernandez@scai.fraunhofer.de',
    )

    add_bel_metadata(graph)

    graph.graph['pathway_id'] = identifier

    nodes = nodes_to_bel(nodes, graph, hgnc_manager, chebi_manager)

    for interaction in interactions:
        participants = interaction['participants']
        interaction_metadata = interaction['metadata']

        add_edges(graph, participants, nodes, interaction_metadata)

    return graph


def nodes_to_bel(nodes: Dict[str, Dict], graph: BELGraph, hgnc_manager: HgncManager, chebi_manager: ChebiManager) -> \
        Dict[str, BaseEntity]:
    """Convert dictionary values to BEL nodes."""
    return {
        node_id: node_to_bel(node_att, graph, hgnc_manager, chebi_manager)
        for node_id, node_att in nodes.items()
    }


def node_to_bel(node: Dict, graph, hgnc_manager: HgncManager, chebi_manager: ChebiManager) -> BaseEntity:
    """Convert node dictionary to BEL node object."""
    node_types = node['entity_type']

    identifier, name, namespace = get_valid_node_parameters(node, hgnc_manager, chebi_manager)
    members = set()

    if namespace == 'hgnc_multiple_entry':
        return composite_abundance(process_multiple_proteins(identifier))

    elif 'Protein' in node_types:
        return protein(namespace=namespace.upper(), name=name, identifier=identifier)

    elif 'Dna' in node_types:
        return gene(namespace=namespace.upper(), name=name, identifier=identifier)

    elif 'Rna' in node_types:
        return rna(namespace=namespace.upper(), name=name, identifier=identifier)

    elif 'SmallMolecule' in node_types:
        return abundance(namespace=namespace.upper(), name=name, identifier=identifier)

    elif 'PhysicalEntity' in node_types:
        return abundance(namespace=namespace.upper(), name=name, identifier=identifier)

    elif 'Complex' in node_types:
        complex_components = node.get('complex_components')

        if complex_components:
            for component in complex_components:
                bel_node = node_to_bel(component, graph, hgnc_manager, chebi_manager)

                members.add(bel_node)

        if members:
            return complex_abundance(
                name=node.get('display_name'),
                members=members,
                identifier=identifier,
                namespace=namespace.upper()
            )
        else:
            return NamedComplexAbundance(
                name=node.get('display_name'),
                identifier=identifier,
                namespace=namespace.upper()
            )

    elif 'Pathway' in node_types:
        bioprocess_node = bioprocess(identifier=identifier, name=name, namespace=namespace.upper())
        graph.add_node_from_data(bioprocess_node)
        return bioprocess_node
    else:
        log.warning('Entity type not recognized', node_types)


def add_edges(graph: BELGraph, participants, nodes, att: Dict):
    """Add edges into the graph."""
    edge_types = att['interaction_type']

    if isinstance(participants, dict):

        reactants = {
            nodes[source_id]
            for source_id in participants['reactants']
        }

        products = {
            nodes[product_id]
            for product_id in participants['products']
        }

        reaction_node = reaction(reactants=reactants, products=products)
        graph.add_node_from_data(reaction_node)

    elif isinstance(participants, tuple):
        u = nodes[participants[0]]
        v = nodes[participants[1]]
        add_simple_edge(graph, u, v, edge_types)


def add_simple_edge(graph: BELGraph, u, v, edge_types):
    """Add a simple edge into the graph."""
    if 'ACTIVATION' in edge_types:
        graph.add_increases(
            u, v,
            citation=REACTOME_CITATION, evidence='Extracted from Reactome',
            object_modifier=activity() if isinstance(v, ACTIVITY_ALLOWED_MODIFIERS) else None,
            annotations={},
        )

    elif 'INHIBITION' in edge_types:
        graph.add_decreases(
            u, v,
            citation=REACTOME_CITATION, evidence='Extracted from Reactome',
            object_modifier=activity() if isinstance(v, ACTIVITY_ALLOWED_MODIFIERS) else None,
            annotations={},
        )
    else:
        log.warning('edge type %s', edge_types)
