# -*- coding: utf-8 -*-

"""This module contains functions to parse KGML files."""

import itertools as itt
from collections import defaultdict
import xml.etree.ElementTree as ET

import networkx as nx
import requests
from bio2bel_kegg.constants import API_KEGG_GET
from bio2bel_kegg.parsers.description import parse_description
from bio2bel_hgnc import Manager as HgncManager
from bio2bel_chebi import Manager as ChebiManager
from ..constants import HGNC

"""Import XML"""


def import_xml_etree(filename):
    """Return XML tree from KGML file.

    :param str filename: path to KGML file
    :returns: XML Tree
    :rtype: xml.etree.ElementTree.ElementTree
    """
    try:
        tree = ET.parse(filename)
    except IOError as ioerr:
        print('File error: ' + str(ioerr))
        return None

    return tree


"""KEGG Handling functions"""


def get_entity_nodes(tree, hgnc_manager, chebi_manager):
    """Find entry elements (KEGG pathway nodes) in XML.

    :param xml.etree.ElementTree.ElementTree tree: XML tree
    :param bio2bel_hgnc.Manager hgnc_manager: HGNC Manager
    :param bio2bel_chebi.Manager chebi_manager: ChEBI Manager
    :return: genes with corresponding metadata (entry_id: [kegg_id, HGNC, UniProt])
    :return: compounds with corresponding metadata (entry_id: [compound_name, ChEBI])
    :rtype: dict
    """
    entry_dict = defaultdict(list)
    compound_dict = defaultdict(list)
    map_dict = defaultdict(list)
    ortholog_dict = defaultdict(list)

    for entry in tree.findall("entry"):

        entry_id = entry.get("id")
        kegg_ids = entry.get("name")
        kegg_type = entry.get("type")

        if kegg_type.startswith('gene'):

            for kegg_id in kegg_ids.split(' '):

                node_info = {
                    'kegg_id': kegg_id,
                    'kegg_type': kegg_type
                }

                kegg_url = API_KEGG_GET.format(kegg_id)

                node_meta_data = parse_description(requests.get(kegg_url))

                if 'DBLINKS' in node_meta_data:

                    for resource, identifier in node_meta_data['DBLINKS']:

                        if resource in {HGNC, 'UniProt'}:

                            if resource == HGNC:

                                node_info['HGNC symbol'] = hgnc_manager.get_gene_by_hgnc_id(identifier).symbol

                            node_info[resource] = identifier

                entry_dict[entry_id].append(node_info)

        elif kegg_type.startswith('compound'):

            compound_info = get_compound_info(kegg_ids, chebi_manager)
            if compound_info:
                compound_dict[entry_id].append(compound_info)

        elif kegg_type.startswith('map'):

            map_info = {'kegg_id': kegg_ids}

            for graphics in entry.iter('graphics'):
                map_name = graphics.get('name')
                map_info['map_name'] = map_name

            map_dict[entry_id].append(map_info)

        elif kegg_type.startswith('ortholog'):

            for ortholog_id in kegg_ids.split(' '):

                ortholog_info = {
                    'kegg_id': ortholog_id,
                    'kegg_type': kegg_type
                }

                ortholog_dict[entry_id].append(ortholog_info)

    return entry_dict, compound_dict, map_dict, ortholog_dict


def get_node_types(tree):
    """Find entity types in XML.

    :param xml.etree.ElementTree.ElementTree tree: XML tree
    :return: count of other entity types present in XML
    :rtype: dict
    """
    undefined_dict = defaultdict(int)

    for entry in tree.findall("entry"):

        entry_type = entry.get("type")

        if entry_type.startswith('gene'):
            continue
        elif entry_type.startswith('compound'):
            continue
        elif entry_type.startswith('group'):
            continue
        else:
            undefined_dict[entry_type] += 1

    return undefined_dict


def get_entities_in_complex(tree, entry_dict):
    """Find all entities participating in a complex in XML.

    :param xml.etree.ElementTree.ElementTree tree: XML tree
    :param dict entry_dict: dictionary of entities with metadata
    :return: set of tuples with kegg ids of entities in complex {((entity1, entity2,...,nth_entity),complex_id)}
    :rtype: set[tuple]
    """
    complexes = defaultdict(list)
    complex_dict = defaultdict(list)
    complex_id_dict = defaultdict(list)
    products = set()

    # get component ids of entities participating in complex
    for entry in (tree.findall("entry")):

        entry_id = entry.get("id")
        entry_type = entry.get("type")

        for component in entry.iter("component"):

            component_id = component.get("id")
            complexes[entry_id].append(component_id)
            complexes[entry_id].append(entry_type)

            # get KEGG ids of component ids
            for k, v in entry_dict.items():

                for node_info in v:

                    if component_id == k:
                        kegg_id = node_info.get('kegg_id', 'kegg_id')
                        complex_dict[component_id].append(kegg_id)

    # get complex ids and kegg ids of entities participating in complex
    for entity_id, components in complexes.items():

        for entity, kegg_ids in complex_dict.items():

            for component in components:

                if component == entity:
                    complex_id_dict[entity_id].append(kegg_ids)

    # get complex
    for k, v in complex_id_dict.items():

        for element in itt.product(*v):
            products.add((element, k))

    return products


def get_complex_components(tree, genes_dict, flattened=False):
    """Get IDs of complex components to construct complexes of protein composites (i.e. similar proteins)
    or get dictionary of flattened lists of all proteins involved in complexes

    :param xml.etree.ElementTree.ElementTree tree: XML tree
    :param dict genes_dict: dictionary of all genes in pathway
    :param bool flattened: True to flatten all complex participants
    :return: dictionary of complex IDs and component IDs
    :return: flattened dictionary of complex IDs and component metadata
    :rtype: dict
    """
    component_info = defaultdict(list)
    complex_ids = defaultdict(list)
    complexes = defaultdict(list)
    all_components = []
    flattened_complexes = defaultdict(list)

    for entry in (tree.findall("entry")):
        entry_id = entry.get("id")

        for component in entry.iter("component"):
            component_id = component.get("id")

            # Get complex IDs and each of their component IDs
            complex_ids[entry_id].append(component_id)

            # Get the IDs of all components participating in complexes
            if component_id not in all_components:
                all_components.append(component_id)

    # Get node info for each component
    for k, v in genes_dict.items():
        for (component_id, node_info) in itt.product(all_components,v):
            if component_id == k:
                component_info[component_id].append(node_info)

    # Flatten lists of components in complexes
    if flattened:

        for k, v in complex_ids.items():
            for comp_id in v:
                for comp_k, comp_v in component_info.items():
                    if comp_id == comp_k:
                        complexes[k].append(comp_v)

        for k, v in complexes.items():
            for component in v:
                for info in component:
                    flattened_complexes[k].append(info)


    return complex_ids, flattened_complexes

"""Get all interactions in KEGG pathways"""


def get_all_relationships(tree):
    """Find all relationships between 2 entities.

    :param xml.etree.ElementTree.ElementTree tree: XML tree
    :return: relationships list [(relation_entry1, relation_entry2, relation_subtype)]
    :rtype: list[tuple]
    """
    relations_list = []

    for relation in tree.findall("relation"):

        relation_entry1 = relation.get("entry1")
        relation_entry2 = relation.get("entry2")
        relation_type = relation.get('type')

        for subtype in relation.iter('subtype'):

            relation_subtype = subtype.get("name")
            relation_value = subtype.get("value")

            if relation_type in {'ECrel', 'PCrel'}:
                relations_list.append((relation_entry1, relation_entry2, 'binding/association'))
                relations_list.append((relation_entry1, relation_value, 'binding/association'))
                relations_list.append((relation_value, relation_entry2, 'binding/association'))

            elif relation_type.startswith('maplink'):
                relations_list.append((relation_entry1, relation_entry2, 'binding/association'))

            else:
                relations_list.append((relation_entry1, relation_entry2, relation_subtype))

    return relations_list


def get_edge_types(tree):
    """Get edge types of relations between 2 entities in XML.

    :param xml.etree.ElementTree.ElementTree tree: XML tree
    :return: count of types of relations present in XML
    :rtype: dict
    """
    edge_types_dict = defaultdict(int)

    for relation in tree.findall("relation"):
        for subtype in relation.iter('subtype'):
            relation_subtype = subtype.get("name")
            edge_types_dict[relation_subtype] += 1

    return edge_types_dict


def get_compound_info(compound_name, chebi_manager):
    """Return information from kegg compound.

    :param str compound_name: compound identifier in kegg
    :param bio2bel_chebi.Manager chebi_manager: ChEBI Manager
    :return: dictionary with compound info
    :rtype: dict
    """
    node_info = {'compound_name': compound_name}

    kegg_url = API_KEGG_GET.format(compound_name.strip('cdp:'))

    node_meta_data = parse_description(requests.get(kegg_url))

    # Adds CHEBI and PubChem identifier to node dictionary
    if 'DBLINKS' in node_meta_data:
        
        for resource, identifier in node_meta_data['DBLINKS']:

            if resource in {'ChEBI', 'PubChem'}:
                node_info[resource] = identifier
                if resource == 'ChEBI':

                    # Split multiple identifiers and get their names
                    for chebi_id in identifier.split(' '):
                        chebi_entry = chebi_manager.get_chemical_by_chebi_id(chebi_id)

                        if not chebi_entry:
                            continue

                        node_info['ChEBI name'] = chebi_entry.name

    return node_info

def get_all_reactions(tree):
    """Find all substrates and products participating in reaction.

    :param xml.etree.ElementTree.ElementTree tree: XML tree
    :return: dictionary with substrate info
    :return: dictionary with product info
    :rtype: dict
    """
    substrates_dict = defaultdict(list)
    products_dict = defaultdict(list)

    for reaction in tree.findall("reaction"):

        reaction_id = reaction.get("id")

        for substrate in reaction.iter('substrate'):
            substrates_dict[reaction_id].append(substrate.get("id"))

        for product in reaction.iter('product'):
            products_dict[reaction_id].append(product.get("id"))

    return substrates_dict, products_dict


def get_reaction_edge_types(tree):
    """Get edge types of reactions between 2 entities in XML.

    :param xml.etree.ElementTree.ElementTree tree: XML tree
    :return: count of types of reactions present in XML
    :rtype: dict
    """
    rxn_edge_types_dict = defaultdict(int)

    for reaction in tree.findall("reaction"):
        reaction_type = reaction.get("type")
        rxn_edge_types_dict[reaction_type] += 1

    return rxn_edge_types_dict


def get_reaction_pathway_edges(xml_tree, substrates_dict, products_dict):
    """Find all reaction edges.

    :param xml.etree.ElementTree.ElementTree xml_tree: xml tree
    :param dict substrates_dict: dictionary with substrate info
    :param dict products_dict: dictionary with product info
    :return: dictionary of reaction elements
    :rtype: dict
    """
    reactions_dict = defaultdict(list)

    for reaction in xml_tree.findall("reaction"):

        reaction_type = reaction.get("type")
        reaction_id = reaction.get("id")

        reaction_substrates = substrates_dict[reaction_id]
        reaction_products = products_dict[reaction_id]

        # Add edge from substrates to products with compound info
        reactions_dict[reaction_id].append((reaction_substrates, reaction_products, reaction_type))

        # If reaction is reversible, flip the reaction order and add a new edge
        if reaction_type == "reversible":
            reactions_dict[reaction_id].append((reaction_products, reaction_substrates, reaction_type))

    return reactions_dict


def get_pathway_edges(node_dict, relations_list, **kwargs):
    """Find all kegg_ids for source and target entries in all interactions and creates
        edges between them with interaction type labels.

    :param dict node_dict: dictionary of all entities and attributes
    :param list relations_list: list of all entities and their interactions
    :param Optional[set] kwargs: set of tuples with kegg ids of entities in complex {((entity1, entity2,...,nth_entity), complex_id)}
    :return: edges set {(entity_1, entity_2), relation_subtype)}
    :rtype: set[tuple]
    """
    source_dict = defaultdict(list)
    target_dict = defaultdict(list)
    edges = set()

    # get kegg ids for entities
    for k, v in node_dict.items():

        for source, target, relation in relations_list:

            for kegg_id in v:
                if not kegg_id['kegg_id'].startswith('undefined'):

                    if source == k:
                        source_dict[source].append(kegg_id['kegg_id'])

                    elif target == k:
                        target_dict[target].append(kegg_id['kegg_id'])

    # add complexes to source and target dictionaries
    if 'set_of_complexes' in kwargs:

        set_of_complexes = kwargs.get('set_of_complexes')

        for group, group_id in set_of_complexes:

            for source, target, relation in relations_list:

                if source == group_id:
                    source_dict[source].append(group)

                if target == group_id:
                    target_dict[target].append(group)

    # create edges between sources and targets
    for source, target, relation in relations_list:

        for source_key, source_value in source_dict.items():

            for target_key, target_value in target_dict.items():

                if source == source_key and target == target_key:
                    product = itt.product(source_value, target_value)

                    for element in product:
                        edge = (element[0], element[1], relation)

                        edges.add(edge)

    return edges


"""Populate empty networkx graph with KEGG pathway entities and interactions"""


def add_nodes_to_graph(graph, node_dict, **kwargs):
    """Add KEGG pathway nodes to networkx graph.

    :param networkx.MultiDiGraph graph: networkx graph where nodes are to be added
    :param dict node_dict: dictionary of all nodes in pathway
    :param Optional[set]: set of complexes
    """
    for entry_id, attrib_list in node_dict.items():

        for attributes in attrib_list:

            if attributes['kegg_id'].startswith('hsa'):
                try:
                    graph.add_node(
                        attributes['kegg_id'],
                        kegg_type=attributes['kegg_type'],
                        entry_id=entry_id,
                        HGNC=attributes['HGNC'],
                        UniProt=attributes['UniProt']
                    )
                except KeyError:
                    graph.add_node(attributes['kegg_id'],
                                   kegg_type=attributes['kegg_type'],
                                   entry_id=entry_id)
                # add complexes to graph
                if 'set_of_complexes' in kwargs:
                    set_of_complexes = kwargs.get('set_of_complexes')

                    for complexes, complex_id in set_of_complexes:
                        graph.add_node(complexes,
                                       kegg_type='complex',
                                       complex_id=complex_id)


def add_reaction_nodes_to_graph(graph, compound_dict):
    """Add KEGG compound nodes to graph.

    :param networkx.MultiDiGraph graph: networkx graph where reaction element nodes are to be added
    :param dict compound_dict: dictionary with compound info
    """
    for compound_key, compound_value in compound_dict.items():
        for compound in compound_value:
            try:
                graph.add_node(compound['compound_name'],
                               kegg_type='compound',
                               entry_id=compound_key,
                               ChEBI=compound['ChEBI'])
            except KeyError:
                graph.add_node(compound['compound_name'],
                               kegg_type='compound',
                               entry_id=compound_key)


def add_edges_to_network(graph, set_of_edges):
    """Add edges from KEGG pathway to a graph in place.

    :param set set_of_edges: tuple list of all edges and relation type (entity1, entity2, typle of relation)
    :param networkx.MultiDiGraph graph: networkx graph where edges are going to be added
    """
    for entry1, entry2, relation in set_of_edges:
        graph.add_edge(entry1, entry2, interaction=relation)


def add_reaction_edges_to_network(graph, reactions_list):
    """Add reaction edges from KEGG pathway to a graph in place.

    :param networkx.MultiDiGraph graph: networkx graph where reaction edges are to be added
    :param reactions_list: list of reactions
    """
    for substrate, product, reaction in reactions_list:

        if len(substrate) == 1 and len(product) == 1:
            graph.add_edge(substrate[0], product[0], interaction=reaction)

        elif len(substrate) == 1 and len(product) > 1:
            graph.add_edge(substrate[0], tuple(product), interaction=reaction)

        elif len(substrate) > 1 and len(product) == 1:
            graph.add_edge(tuple(substrate), product[0], interaction=reaction)

        else:
            graph.add_edge(tuple(substrate), tuple(product), interaction=reaction)

    # TODO: eliminate repeated rxn edges


def populate_graph(path):
    """Convert KGML file to a networkx graph.

    :param str path: path to KGML file
    :rtype: networkx.MultiDiGraph
    """
    graph = nx.MultiDiGraph()

    # Load xml
    xml_tree = import_xml_etree(path)

    # Initialize HgncManager
    hgnc_manager = HgncManager()

    # Initialize ChebiManager
    chebi_manager = ChebiManager()

    # Parse file and get entities and interactions
    genes_dict, compounds_dict, maps_dict, orthologs_dict = get_entity_nodes(xml_tree, hgnc_manager, chebi_manager)
    interactions = get_all_relationships(xml_tree)

    # Get compounds and reactions
    substrates_dict, products_dict = get_all_reactions(xml_tree)
    reactions = get_reaction_pathway_edges(xml_tree, substrates_dict, products_dict)

    # Get complexes and interactions
    complexes = get_entities_in_complex(xml_tree, genes_dict)
    edges = get_pathway_edges(genes_dict, interactions, set_of_complexes=complexes)

    # Add the nodes and edges to the graph
    add_nodes_to_graph(graph, genes_dict, set_of_complexes=complexes)
    add_reaction_nodes_to_graph(graph, compounds_dict)
    add_edges_to_network(graph, edges)
    add_reaction_edges_to_network(graph, reactions)

    return graph
