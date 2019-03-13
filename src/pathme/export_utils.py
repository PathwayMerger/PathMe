# -*- coding: utf-8 -*-

"""Export harmonized universe."""

import logging
import os
from typing import List

from pathme.constants import KEGG, REACTOME, WIKIPATHWAYS
from pathme.normalize_names import normalize_graph_names
from pathme.pybel_utils import flatten_complex_nodes
from pybel import BELGraph
from pybel import from_pickle
from pybel import union
from tqdm import tqdm

logger = logging.getLogger(__name__)


def get_all_pickles(kegg_path, reactome_path, wikipathways_path):
    """Return a list with all pickle paths."""
    kegg_pickles = get_files_in_folder(kegg_path)

    if not kegg_pickles:
        logger.warning('No KEGG files found. Please create the BEL KEGG files')

    reactome_pickles = get_files_in_folder(reactome_path)

    if not reactome_pickles:
        logger.warning('No Reactome files found. Please create the BEL Reactome files')

    wp_pickles = get_files_in_folder(wikipathways_path)

    if not wp_pickles:
        logger.warning('No WikiPathways files found. Please create the BEL WikiPathways files')

    return kegg_pickles, reactome_pickles, wp_pickles


def get_universe_graph(
        kegg_path: str, reactome_path: str, wikipathways_path: str,
        flatten: bool = True, normalize_names: bool = True) -> BELGraph:
    """Return universe graph."""
    kegg_pickles, reactome_pickles, wp_pickles = get_all_pickles(kegg_path, reactome_path, wikipathways_path)

    all_pickles = kegg_pickles + reactome_pickles + wp_pickles

    logger.info(f'A total of {len(all_pickles)} will be merged into the universe')

    iterator = tqdm(all_pickles, desc='Loading of the graph pickles')

    universe_list = []

    # Export KEGG
    for file in iterator:
        if not file.endswith('.pickle'):
            continue

        if file in kegg_pickles:
            graph = from_pickle(os.path.join(kegg_path, file))

            if flatten:
                flatten_complex_nodes(graph)

            if normalize_names:
                normalize_graph_names(graph, KEGG)

        elif file in reactome_pickles:
            graph = from_pickle(os.path.join(reactome_path, file))

            if flatten:
                flatten_complex_nodes(graph)

            if normalize_names:
                normalize_graph_names(graph, REACTOME)

        elif file in wp_pickles:
            graph = from_pickle(os.path.join(wikipathways_path, file))

            if flatten:
                flatten_complex_nodes(graph)

            if normalize_names:
                normalize_graph_names(graph, WIKIPATHWAYS)
        else:
            logger.warning(f'Unknown pickle file: {file}')
            continue

        universe_list.append(graph)

    logger.info('Merging all into a hairball...')

    return union(universe_list)


def get_files_in_folder(path: str) -> List[str]:
    """Return the files in a given folder.

    :param path: folder path
    :return: file names in folder
    """
    return [
        file
        for file in os.listdir(path)
        if os.path.isfile(os.path.join(path, file))
    ]
