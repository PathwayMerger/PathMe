# -*- coding: utf-8 -*-

"""This module has utilities method for parsing and handling KEGG KGML files."""

import os
from collections import defaultdict
import pandas as pd
import requests
import tqdm
import networkx as nx
from bio2bel_kegg import Manager as KeggManager
from pathme.wikipathways.utils import get_files_in_folder
from pathme.kegg.kegg_xml_parser import import_xml_etree, populate_graph, get_node_types, get_edge_types, get_reaction_edge_types
from ..constants import DATA_DIR, KEGG, KEGG_KGML_URL


def get_kegg_pathway_ids(connection=None):
    """Return a list of all pathway identifiers stored in the KEGG database.

    :param Optional[str] connection: connection to the database
    :returns: list of all kegg_pathway_ids
    :rtype: list
    """
    kegg_manager = KeggManager(connection=connection)
    kegg_pathways_ids = [
        pathway.resource_id.replace('path:', '')
        for pathway in kegg_manager.get_all_pathways()
    ]

    if not kegg_pathways_ids:
        raise EnvironmentError('Your database is empty. Please run python3 -m bio2bel_kegg populate')

    return kegg_pathways_ids


def download_kgml_files(kegg_pathway_ids):
    """Downloads KEGG KGML files by querying the KEGG API.

    :param list kegg_pathway_ids: list of kegg ids
    """
    for kegg_id in tqdm.tqdm(kegg_pathway_ids, desc='Downloading KEGG files'):
        request = requests.get(KEGG_KGML_URL.format(kegg_id))
        with open(os.path.join(DATA_DIR, KEGG, '{}.xml'.format(kegg_id)), 'w+') as file:
            file.write(request.text)
            file.close()



def parse_kegg(path):
    """Parse a folder and returns graph objects.

    :param str path: path to folder containing XML files
    :rtype: list[networkx.MultiDiGraph]
    """
    pathways = []

    files = get_files_in_folder(path)

    for file_name in tqdm.tqdm(files, desc='Parsing KEGG files'):
        file_path = os.path.join(path, file_name)

        pathways.append((file_name, populate_graph(file_path)))

    return pathways


def get_undefined_node_statistics_in_xml(path):
    """Parse a folder and get additional entity type stats for each pathway in XML.

    :param str path: path to folder containing XML files
    :return additional entity node types in XML
    :rtype: pandas.DataFrame
    """
    additional_node_types = []

    files = get_files_in_folder(path)

    for file_name in tqdm.tqdm(files, desc='Parsing KEGG files for undefined entity types'):
        tree = import_xml_etree(os.path.join(path, file_name))

        # Get additional node types dictionary with statistics
        undefined_nodes_dict = get_node_types(tree)
        undefined_nodes_dict['filename'] = file_name.strip('.xml')
        additional_node_types.append(undefined_nodes_dict)

    # Add pathway IDs with additional node types data to DataFrame
    df = pd.DataFrame.from_dict(additional_node_types)
    df = df.set_index('filename')

    return df


def get_pathway_entity_types(pathway):
    """Get entity type stats for all nodes in networkx pathway

    :param str pathway: graph object
    :return: count of all entity types present in graph
    :rtype: dict
    """
    entity_types_in_pathway = defaultdict(int)

    # Get node type attribute from graph
    kegg_path = nx.get_node_attributes(pathway, 'kegg_type')

    for key, value in kegg_path.items():
        entity_types_in_pathway[value] += 1

    return entity_types_in_pathway


def get_node_type_statistics_in_network(pathways):
    """Get entity type stats for all nodes in graph objects.

    :param list pathways: list of pathway IDs and networkx pathways
    :return all entity node types in graphs
    :rtype: pandas.DataFrame
    """
    path_types = []

    # Get node types dictionary with statistics
    for file_name, path in pathways:
        entity_type_dict = get_pathway_entity_types(path)
        entity_type_dict['filename'] = file_name.strip('.xml')
        path_types.append(entity_type_dict)

    # Add pathway ID with node type info to dataframe
    kegg_nodes_df = pd.DataFrame.from_dict(path_types)
    kegg_nodes_df = kegg_nodes_df.set_index('filename')

    return kegg_nodes_df


def get_xml_relation_type_statistics(path):
    """Parse a folder and get interaction type statistics in XML.

    :param str path: path to folder containing XML files
    :return: relation edge types in XML
    :rtype: pandas.DataFrame
    """
    edge_types = []

    files = get_files_in_folder(path)

    for file_name in tqdm.tqdm(files, desc='Parsing KEGG files for relation types'):
        tree = import_xml_etree(os.path.join(path, file_name))
        # Get dictionary of all edge types in XML
        relation_type_dict = get_edge_types(tree)
        relation_type_dict['filename'] = file_name.strip('.xml')
        edge_types.append(relation_type_dict)

    # Add pathway id with additional node types data to dataframe
    relations_df = pd.DataFrame.from_dict(edge_types)
    relations_df = relations_df.set_index('filename')

    return relations_df


def get_xml_reaction_type_statistics(path):
    """Parse a folder and get reaction interaction type statistics in XML.

    :param str path: path to folder containing XML files
    :return: reaction edge types in XML
    :rtype: pandas.DataFrame
    """
    reaction_edge_types = []

    files = get_files_in_folder(path)

    for file_name in tqdm.tqdm(files, desc='Parsing KEGG files for relation types'):
        tree = import_xml_etree(os.path.join(path, file_name))
        # Get dictionary of all edge types in XML
        reaction_type_dict = get_reaction_edge_types(tree)
        reaction_type_dict['filename'] = file_name.strip('.xml')
        reaction_edge_types.append(reaction_type_dict)

    # Add pathway id with additional node types data to dataframe
    reactions_df = pd.DataFrame.from_dict(reaction_edge_types)
    reactions_df = reactions_df.set_index('filename')

    return reactions_df