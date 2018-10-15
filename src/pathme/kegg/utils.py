# -*- coding: utf-8 -*-

"""This module has utilities method for parsing and handling KEGG KGML files."""

import os

import pandas as pd
import requests
import tqdm
from bio2bel_kegg.manager import Manager as KeggManager

from pathme.kegg.convert_to_bel import get_bel_types
from pathme.kegg.kegg_xml_parser import import_xml_etree, get_xml_types
from pathme.wikipathways.utils import get_files_in_folder
from pathme.constants import KEGG_FILES, KEGG_KGML_URL, KEGG_STATS_COLUMN_NAMES

__all__ = [
    'download_kgml_files',
    'get_kegg_statistics',
    'get_kegg_pathway_ids'
]


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
    """Download KEGG KGML files by querying the KEGG API.

    :param list kegg_pathway_ids: list of kegg ids
    """
    for kegg_id in tqdm.tqdm(kegg_pathway_ids, desc='Downloading KEGG files'):
        request = requests.get(KEGG_KGML_URL.format(kegg_id))
        with open(os.path.join(KEGG_FILES, '{}.xml'.format(kegg_id)), 'w+') as file:
            file.write(request.text)
            file.close()


def get_kegg_statistics(path, hgnc_manager, chebi_manager, flatten=None):
    """Parse a folder and get KEGG statistics.

    :param graph: path
    :param bio2bel_hgnc.Manager hgnc_manager: HGNC manager
    :param bio2bel_chebi.Manager chebi_manager: ChEBI manager
    :param str path: path to folder containing XML files
    :return: KEGG KGML file and BEL graph statistics
    :rtype: pandas.DataFrame
    """
    df = pd.DataFrame()
    export_file_name = 'KEGG_pathway_stats_{}.csv'.format('flatten' if flatten else 'non_flatten')

    # Get list of all files in folder
    files = get_files_in_folder(path)

    for file_name in tqdm.tqdm(files, desc='Parsing KGML files and BEL graphs for entities and relation stats'):
        pathway_names = []
        file_path = os.path.join(path, file_name)
        tree = import_xml_etree(file_path)
        root = tree.getroot()
        pathway_names.append(root.attrib['title'])

        # Get dictionary of all entity and interaction types in XML
        xml_statistics_dict = get_xml_types(tree)

        # Get dictionary of all node and edge types in BEL Graph
        bel_statistics_dict = get_bel_types(file_path, hgnc_manager, chebi_manager, flatten=flatten)

        # Get dictionary with all XML and BEL graph stats
        xml_statistics_dict.update(bel_statistics_dict)

        # Update dictionary of all XML and BEL graph stats with corresponding column names
        all_kegg_statistics = {
            KEGG_STATS_COLUMN_NAMES[key]: value
            for key, value in xml_statistics_dict.items()
        }

        # Add pathway statistic rows to DataFrame
        pathway_data = pd.DataFrame(
            all_kegg_statistics,
            index=pathway_names,
            columns=KEGG_STATS_COLUMN_NAMES.values(),
            dtype=int

        )
        df = df.append(pathway_data.fillna(0).astype(int))

    df.to_csv(export_file_name, sep='\t')
    return df
