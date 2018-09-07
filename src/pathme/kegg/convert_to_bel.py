# -*- coding: utf-8 -*-

"""This module contains the methods to convert a KEGG RDF network into a BELGraph."""
import logging
from collections import defaultdict
from itertools import product

from bio2bel_chebi import Manager as ChebiManager
from bio2bel_hgnc import Manager as HgncManager
from pybel import BELGraph
from pybel.dsl.edges import activity
from pybel.dsl.nodes import abundance, bioprocess, complex_abundance, composite_abundance, protein, pmod, reaction

from pathme.constants import CHEBI, HGNC, KEGG_CITATION, KEGG_MODIFICATIONS, KEGG
from pathme.kegg.kegg_xml_parser import (
    get_all_reactions,
    get_all_relationships,
    get_entity_nodes,
    get_complex_components,
    get_reaction_pathway_edges,
    import_xml_etree
)

log = logging.getLogger(__name__)

"""Populate empty BEL graph with KEGG pathway entities and interactions"""


def kegg_to_bel(path):
    """Convert KGML file to a BELGraph.

    :param str path: path to KGML file
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
        contact='daniel.domingo.fernandez@scai.fraunhofer.de',
    )

    # Initialize HgncManager
    hgnc_manager = HgncManager()

    # Initialize ChebiManager
    chebi_manager = ChebiManager()

    # Parse file and get entities
    genes_dict, compounds_dict, maps_dict, orthologs_dict = get_entity_nodes(xml_tree, hgnc_manager, chebi_manager)

    # Get complexes
    complex_ids, flattened_complexes = get_complex_components(xml_tree, genes_dict, flattened=True)

    # Get interactions
    relations_list = get_all_relationships(xml_tree)

    # Get compounds and reactions
    substrates_dict, products_dict = get_all_reactions(xml_tree)
    reactions_dict = get_reaction_pathway_edges(xml_tree, substrates_dict, products_dict)

    # Add nodes and edges to graph
    nodes = xml_entities_to_bel(genes_dict, compounds_dict, maps_dict, flattened=True)
    nodes = xml_complexes_to_bel(nodes, complex_ids, flatten_complexes=flattened_complexes)
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


def xml_complexes_to_bel(node_dict, complex_ids, **kwargs):
    """ Convert complexes in XML to BEL nodes where each complex is made up of proteins
    and/or composites (i.e. groups of related proteins).

    :param dict node_dict: dictionary of BEL nodes
    :param dict complex_ids: dictionary of complex IDs and component IDs
    :param Optional[dict] kwargs: dictionary of complex IDs and flattened list of all components
    :return: dictionary of BEL nodes
    :rtype: dict
    """
    member_dict = defaultdict(list)

    if 'flatten_complexes' in kwargs:
        flatten_complexes = kwargs.get('flatten_complexes')

        for node_id, node_att in flatten_complexes.items():
            node_dict[node_id] = flatten_complex_to_bel_node(node_att)

    else:
        for node_k, node_v in node_dict.items():
            for complex_k, complex_v in complex_ids.items():
                for component in complex_v:
                    if component == node_k:
                        member_dict[complex_k].append(node_v)

        for k, v in member_dict.items():
            node_dict[k] = complex_abundance(v)

    return node_dict


def gene_to_bel_node(graph, node):
    """Create a protein or protein composite BEL node and add to BEL Graph.

    :param graph: BELGraph
    :param dict node: dictionary of node attributes
    :return: BEL node dictionary
    :rtype: dict
    """
    members = set()

    # Create a protein BEL node
    if len(node) == 1:
        for attribute in node:

            if HGNC in attribute:
                protein_node = protein(namespace=HGNC, name=attribute[HGNC], identifier=attribute['HGNC symbol'])
                graph.add_node_from_data(protein_node)
                return (protein_node)

            elif 'UniProt' in attribute:
                protein_node = protein(namespace='UniProt', name=attribute['UniProt'], identifier=attribute['UniProt'])
                graph.add_node_from_data(protein_node)
                return (protein_node)

            else:
                protein_node = protein(namespace=KEGG, name=attribute['kegg_id'], identifier=attribute['kegg_id'])
                graph.add_node_from_data(protein_node)
                return (protein_node)

    # Create a composite abundance BEL node
    else:
        for member in node:
            bel_node = gene_to_bel_node([member])
            members.add(bel_node)

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
    proteins_list = []

    for attribute in node:

        # Create a protein BEL node
        if len(node) == 1:

            if HGNC in attribute:
                protein_node = protein(namespace=HGNC, name=attribute['HGNC symbol'], identifier=attribute['HGNC'])
                graph.add_node_from_data(protein_node)
                return (protein_node)

            elif 'UniProt' in attribute:
                protein_node = protein(namespace='UniProt', name=attribute['UniProt'], identifier=attribute['UniProt'])
                graph.add_node_from_data(protein_node)
                return (protein_node)

            else:
                protein_node = protein(namespace=KEGG, name=attribute['kegg_id'], identifier=attribute['kegg_id'])
                graph.add_node_from_data(protein_node)
                return (protein_node)

        # Create a list of protein BEL nodes
        else:

            if HGNC in attribute:
                protein_node = protein(namespace=HGNC, name=attribute['HGNC symbol'], identifier=attribute['HGNC'])
                graph.add_node_from_data(protein_node)
                proteins_list.append(protein_node)

            elif 'UniProt' in attribute:
                protein_node = protein(namespace='UniProt', name=attribute['UniProt'], identifier=attribute['UniProt'])
                graph.add_node_from_data(protein_node)
                proteins_list.append(protein_node)

            else:
                protein_node = protein(namespace=KEGG, name=attribute['kegg_id'], identifier=attribute['kegg_id'])
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


def flatten_complex_to_bel_node(node):
    """Create complex abundance BEL node.

    :param dict node: dictionary of node attributes
    :return: BEL node dictionary
    :rtype: dict
    """
    members = set()
    for attributes in node:
        identifier = attributes[HGNC]
        name = attributes['HGNC symbol']
        member = protein(namespace=HGNC, name=name, identifier=identifier)
        members.add(member)

    return complex_abundance(members=members)


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

        # If entity is a list of proteins, add an edge to each protein node in list
        if not isinstance(type(u), complex_abundance) and not isinstance(type(v), complex_abundance):

            for member, component in product(u, v):
                if type(member) != str:
                    if type(component) != str:
                        add_simple_edge(graph, member, component, relation)

        # If source is a list of proteins and target is a complex, add an edge between them
        elif not isinstance(type(u), complex_abundance) and isinstance(type(v), complex_abundance):
            for member in u:
                if type(member) != str:
                    add_simple_edge(graph, member, v, relation)

        # If source is a complex and target is a list of proteins, add an edge between them
        elif isinstance(type(u), complex_abundance) and not isinstance(type(v), complex_abundance):
            for component in v:
                if type(component) != str:
                    add_simple_edge(graph, u, component, relation)
                    
        # If entities are both complexes, add an edge between them
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

                # Get enzyme-reaction edges
                for gene_type in enzyme:
                    add_simple_edge(graph, gene_type, reaction_node, reaction_type)


def add_simple_edge(graph, u, v, relation_type):
    """Add corresponding edge type to BEL graph.

    :param graph: BELGraph
    :param u: source node
    :param v: target node
    :param list relation_type: list of entity IDs and types of relations
    """
    # Subject activity increases protein modification of object
    if relation_type in {'phosphorylation', 'glycosylation', 'ubiquitination', 'meythylation'}:
        v = v.with_variants(pmod(KEGG_MODIFICATIONS[relation_type]))
        graph.add_increases(u, v, citation='', evidence='', subject_modifier=activity())

    # Subject activity decreases protein modification (i.e. dephosphorylation) of object
    elif relation_type == 'dephosphorylation':
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
        v = v.get_rna()
        graph.add_increases(u, v, citation=KEGG_CITATION, evidence='')

    # Subject decreases expression of object
    elif relation_type == 'repression':

        # Repression object is converted to RNA abundance
        v = v.get_rna()
        graph.add_decreases(u, v, citation=KEGG_CITATION, evidence='')

    elif relation_type in {'dissociation', 'hidden compound', 'missing interaction', 'state change'}:
        pass

    else:
        raise ValueError('Unexpected relation type {}'.format(relation_type))
