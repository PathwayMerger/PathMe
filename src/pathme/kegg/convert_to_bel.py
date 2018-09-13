# -*- coding: utf-8 -*-

"""This module contains the methods to convert a KEGG RDF network into a BELGraph."""
import logging
from collections import defaultdict
from itertools import product

from pathme.constants import CHEBI, HGNC, KEGG_CITATION, KEGG_MODIFICATIONS, KEGG
from pathme.kegg.kegg_xml_parser import (
    get_all_reactions,
    get_all_relationships,
    get_entity_nodes,
    get_complex_components,
    get_reaction_pathway_edges,
    import_xml_etree
)
from pybel import BELGraph
from pybel.dsl.edges import activity
from pybel.dsl.node_classes import CentralDogma
from pybel.dsl.nodes import abundance, bioprocess, complex_abundance, composite_abundance, protein, pmod, reaction
from pybel.struct.summary import count_functions, edge_summary

log = logging.getLogger(__name__)

"""Populate empty BEL graph with KEGG pathway entities and interactions"""


def kegg_to_bel(path, hgnc_manager, chebi_manager, flatten=False):
    """Convert KGML file to a BELGraph.

    :param str path: path to KGML file
    :param bio2bel_hgnc.Manager hgnc_manager: HGNC manager
    :param bio2bel_chebi.Manager chebi_manager: ChEBI manager
    :param bool flatten: flat nodes
    :rtype: BELGraph
    """
    # Load xml
    xml_tree = import_xml_etree(path)
    root = xml_tree.getroot()

    graph = BELGraph(
        name=root.attrib['title'],
        version='1.0.0',
        description=root.attrib['link'],
        pathway_id=root.attrib['name'],
        authors="Daniel Domingo-Fernández, Josep Marín-Llaó and Sarah Mubeen",
        contact='daniel.domingo.fernandez@scai.fraunhofer.de'
    )

    # Parse file and get entities and interactions
    genes_dict, compounds_dict, maps_dict, orthologs_dict = get_entity_nodes(xml_tree, hgnc_manager, chebi_manager)
    relations_list = get_all_relationships(xml_tree)

    # Get compounds and reactions
    substrates_dict, products_dict = get_all_reactions(xml_tree, compounds_dict)
    reactions_dict = get_reaction_pathway_edges(xml_tree, substrates_dict, products_dict)

    # Get complexes
    complex_ids, flattened_complexes = get_complex_components(xml_tree, genes_dict, flattened=flatten)

    # Add nodes and edges to graph
    nodes = xml_entities_to_bel(graph, genes_dict, compounds_dict, maps_dict, flattened=flatten)

    nodes = xml_complexes_to_bel(
        graph=graph,
        node_dict=nodes,
        complex_ids=complex_ids,
        flatten_complexes=flattened_complexes if flatten else None
    )

    add_edges(graph, relations_list, nodes)
    add_reaction_edges(graph, reactions_dict, nodes)

    return graph


"""Get all entities from XML tree and convert to BEL nodes"""


def xml_entities_to_bel(graph, genes_dict, compounds_dict, maps_dict, flattened=False):
    """Convert gene and compound entities in XML to BEL nodes.

    :param graph: BELGraph
    :param dict genes_dict: dictionary of genes in XML
    :param dict compounds_dict: dictionary of compounds in XML
    :param dict maps_dict: dictionary of pathway maps in XML
    :param bool flattened: True to flatten to list of similar genes grouped together
    :return: dictionary of BEL nodes
    :rtype: dict
    """
    if flattened:
        node_dict = {
            node_id: flatten_gene_to_bel_node(graph, node_att)
            for node_id, node_att in genes_dict.items()
        }
    else:
        node_dict = {
            node_id: gene_to_bel_node(graph, node_att)
            for node_id, node_att in genes_dict.items()
        }

    for node_id, node_att in compounds_dict.items():
        node_dict[node_id] = compound_to_bel(graph, node_att)

    for node_id, node_att in maps_dict.items():
        node_dict[node_id] = map_to_bel_node(graph, node_att)

    return node_dict


def xml_complexes_to_bel(graph, node_dict, complex_ids, flatten_complexes=None):
    """ Convert complexes in XML to BEL nodes where each complex is made up of proteins
    and/or composites (i.e. groups of related proteins)

    :param dict node_dict: dictionary of BEL nodes
    :param dict complex_ids: dictionary of complex IDs and component IDs
    :param Optional[dict]: dictionary of complex IDs and flattened list of all components
    :return: dictionary of BEL nodes
    :rtype: dict
    """
    member_dict = defaultdict(list)

    if flatten_complexes is not None:
        for node_id, node_att in flatten_complexes.items():
            node_dict[node_id] = flatten_complex_to_bel_node(graph, node_att)

    # For all complexes, add BEL node component info
    else:
        for complex_id, member_ids in complex_ids.items():
            for member in member_ids:
                member_dict[complex_id].append(node_dict[member])

        for complex_id, bel_members in member_dict.items():
            node_dict[complex_id] = complexes_to_bel_node(graph, bel_members)

    return node_dict


def complexes_to_bel_node(graph, members):
    complex_node = complex_abundance(members=members)
    graph.add_node_from_data(complex_node)

    return complex_node


def gene_to_bel_node(graph, node):
    """Create a protein or protein composite BEL node and add to BEL Graph.

    :param graph: BELGraph
    :param list[dict] node: dictionary of node attributes
    :return: corresponding BEL node
    :rtype: pybel.dsl.BaseEntity
    """
    members = list()

    # Create a protein BEL node
    if len(node) == 1:
        for attribute in node:

            if HGNC in attribute:
                protein_node = protein(namespace=HGNC, name=attribute['HGNC symbol'], identifier=attribute[HGNC])
                graph.add_node_from_data(protein_node)
                return protein_node

            elif 'UniProt' in attribute:
                protein_node = protein(namespace='UniProt', name=attribute['UniProt'], identifier=attribute['UniProt'])
                graph.add_node_from_data(protein_node)
                return protein_node

            else:
                protein_node = protein(namespace=KEGG, name=attribute['kegg_id'], identifier=attribute['kegg_id'])
                graph.add_node_from_data(protein_node)
                return protein_node

    # Create a composite abundance BEL node
    else:
        for member in node:
            bel_node = gene_to_bel_node(graph, [member])
            members.append(bel_node)

        protein_composite = composite_abundance(members=members)
        graph.add_node_from_data(protein_composite)
        return protein_composite


def flatten_gene_to_bel_node(graph, node):
    """Create a protein or list of protein BEL nodes and add to BEL Graph.

    :param graph: BELGraph
    :param dict node: dictionary of node attributes
    :return: BEL node dictionary
    :rtype: dict
    """
    # if only 1 protein node, return corresponding BEL node
    if len(node) == 1:
        node_dict = node[0]

        if HGNC in node_dict:
            protein_node = protein(namespace=HGNC, name=node_dict['HGNC symbol'], identifier=node_dict[HGNC])
            graph.add_node_from_data(protein_node)
            return protein_node

        else:
            protein_node = protein(namespace=KEGG, name=node_dict['kegg_id'], identifier=node_dict['kegg_id'])
            graph.add_node_from_data(protein_node)
            return protein_node

    proteins_list = []
    # if multiple protein nodes, return corresponding list of BEL nodes
    for node_dict in node:

        if HGNC in node_dict:
            protein_node = protein(namespace=HGNC, name=node_dict['HGNC symbol'], identifier=node_dict[HGNC])
            graph.add_node_from_data(protein_node)
            proteins_list.append(protein_node)

        else:
            protein_node = protein(namespace=KEGG, name=node_dict['kegg_id'], identifier=node_dict['kegg_id'])
            graph.add_node_from_data(protein_node)
            proteins_list.append(protein_node)

    return proteins_list


def compound_to_bel(graph, node):
    """Create an abundance BEL node.

    :param graph: BELGraph
    :param dict node: dictionary of node attributes
    :return: BEL node dictionary
    :rtype: dict
    """
    for attribute in node:

        if CHEBI in attribute:

            identifier = attribute[CHEBI]
            name = attribute['ChEBI name']
            namespace = CHEBI

            compound = abundance(namespace=namespace, name=name, identifier=identifier)
            graph.add_node_from_data(compound)
            return compound

        else:

            identifier = attribute['PubChem']
            name = attribute['PubChem']
            namespace = 'PubChem'

            compound = abundance(namespace=namespace, name=name, identifier=identifier)
            graph.add_node_from_data(compound)
            return compound


def map_to_bel_node(graph, node):
    """Create a biological process BEL node.

    :param graph: BELGraph
    :param dict node: dictionary of node attributes
    :return: BEL node dictionary
    :rtype: dict
    """
    for attribute in node:
        name = attribute['map_name']
        identifier = attribute['kegg_id']

        bio_process = bioprocess(namespace=KEGG, name=name, identifier=identifier)
        graph.add_node_from_data(bio_process)
        return bio_process


def flatten_complex_to_bel_node(graph, node):
    """Create complex abundance BEL node.

    :param dict node: dictionary of node attributes
    :return: BEL node dictionary
    :rtype: dict
    """
    members = list()

    for attribute in node:

        if HGNC in attribute:
            protein_node = protein(namespace=HGNC, name=attribute['HGNC symbol'], identifier=attribute[HGNC])
            members.append(protein_node)

        else:
            protein_node = protein(namespace=KEGG, name=attribute['kegg_id'], identifier=attribute['kegg_id'])
            members.append(protein_node)

    complex_members = complex_abundance(members=members)
    graph.add_node_from_data(complex_members)

    return complex_members


"""Get edges between BEL nodes"""


def add_edges(graph, edges, nodes):
    """Add edges to BEL graph.

    :param graph: BELGraph
    :param list edges: list of relationships with entity IDs and interaction types
    :param dict nodes: dictionary of BEL nodes
    """
    for source, target, relation in edges:

        u = nodes[source]
        v = nodes[target]

        # If subject and object are lists, create edges between all products
        if isinstance(u, list) and isinstance(v, list):
            for pair in product(u, v):
                add_simple_edge(graph, pair[0], pair[1], relation)

        # If source is protein list and target is not, add edges between members in list and target
        elif isinstance(u, list) and not isinstance(v, list):
            for member in u:
                add_simple_edge(graph, member, v, relation)

        # If source is not a list and target is proteins list, add edges between them
        elif not isinstance(u, list) and isinstance(v, list):
            for member in v:
                add_simple_edge(graph, u, member, relation)

        # If entities are not lists, add edges between them
        else:
            add_simple_edge(graph, u, v, relation)


def add_reaction_edges(graph, reaction_dict, nodes):
    """Add edges from reactants to products and enzymes to reactions to BEL Graph.

    :param graph: BELGraph
    :param dict reaction_dict: dictionary of reaction IDs and reactant and product IDs
    :param dict nodes: dictionary of BEL nodes
    """
    for k, v in reaction_dict.items():

        # Get BEL gene node(s)
        enzyme = nodes[k]

        # Get compound nodes
        for source, target, reaction_type in v:

            reactants_list = []
            products_list = []

            # Get reactant compound node
            for source_id in source:
                substrate = nodes[source_id]
                reactants_list.append(substrate)

            # Get product compound node
            for target_id in target:
                product = nodes[target_id]
                products_list.append(product)

                # Add reaction BEL node to graph
                reaction_node = reaction(reactants=reactants_list, products=products_list)
                graph.add_node_from_data(reaction_node)

                if isinstance(enzyme, list):
                    for gene_type in enzyme:
                        add_simple_edge(graph, gene_type, reaction_node, reaction_type)
                else:
                    add_simple_edge(graph, enzyme, reaction_node, reaction_type)


def add_simple_edge(graph, u, v, relation_type):
    """Add corresponding edge type to BEL graph.

    :param graph: BELGraph
    :param u: source node
    :param v: target node
    :param list relation_type: list of entity IDs and types of relations
    """
    # Subject activity increases protein modification of object
    if relation_type in {'phosphorylation', 'glycosylation', 'ubiquitination', 'meythylation'}:

        # If the object is a gene, miRNA, RNA, or protein, add protein modification
        if isinstance(v, CentralDogma):
            v = v.with_variants(pmod(KEGG_MODIFICATIONS[relation_type]))
        graph.add_increases(u, v, citation='', evidence='', subject_modifier=activity())

    # Subject activity decreases protein modification (i.e. dephosphorylation) of object
    elif relation_type == 'dephosphorylation':

        # If the object is a gene, miRNA, RNA, or protein, add protein modification
        if isinstance(v, CentralDogma):
            v = v.with_variants(pmod('Ph'))
        graph.add_decreases(u, v, citation=KEGG_CITATION, evidence='', subject_modifier=activity())

    # Subject increases activity of object
    elif relation_type == 'activation':
        graph.add_increases(u, v, citation=KEGG_CITATION, evidence='', object_modifier=activity())

    # Catalytic activity of subject increases transformation of reactant(s) to product(s)
    elif relation_type in {'reversible', 'irreversible'}:
        graph.add_increases(u, v, citation='KEGG_CITATION', evidence='', subject_modifier=activity('cat'))

    # Subject decreases activity of object
    elif relation_type == 'inhibition':
        graph.add_decreases(u, v, citation=KEGG_CITATION, evidence='', object_modifier=activity())

    # Indirect effect and binding/association are noted to be equivalent relation types
    elif relation_type in {'indirect effect', 'binding/association'}:
        graph.add_association(u, v, citation=KEGG_CITATION, evidence='')

    # Subject increases expression of object
    elif relation_type == 'expression':

        # Expression object is converted to RNA abundance
        if isinstance(v, CentralDogma):
            v = v.get_rna()
        graph.add_increases(u, v, citation=KEGG_CITATION, evidence='')

    # Subject decreases expression of object
    elif relation_type == 'repression':

        # Repression object is converted to RNA abundance
        if isinstance(v, CentralDogma):
            v = v.get_rna()
        graph.add_decreases(u, v, citation=KEGG_CITATION, evidence='')

    elif relation_type in {'dissociation', 'hidden compound', 'missing interaction', 'state change'}:
        pass

    else:
        raise ValueError('Unexpected relation type {}'.format(relation_type))


def get_bel_types(path, hgnc_manager, chebi_manager, flatten=False):
    """Get all BEL node and edge type statistics.

    :param str path: path to KGML file
    :param bio2bel_hgnc.Manager hgnc_manager: HGNC manager
    :param bio2bel_chebi.Manager chebi_manager: ChEBI manager
    :param bool flatten: flat nodes
    :return: count of all nodes and edges in BEL graph
    :rtype: dict
    """
    bel_stats = {}

    bel_graph = kegg_to_bel(path, hgnc_manager, chebi_manager, flatten=True if flatten else False)

    bel_stats['nodes'] = bel_graph.number_of_nodes()
    bel_stats['edges'] = bel_graph.number_of_edges()

    # Get count of all BEL function types
    bel_functions_dict = count_functions(bel_graph)
    bel_stats.update(bel_functions_dict)

    # Get count of all BEL edge types
    bel_edges_dict = edge_summary.count_relations(bel_graph)
    bel_stats.update(bel_edges_dict)

    return bel_stats
