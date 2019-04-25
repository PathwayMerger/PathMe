# -*- coding: utf-8 -*-

"""Export harmonized universe."""

import logging
import os
from typing import List

from tqdm import tqdm

from pathme.constants import KEGG, REACTOME, WIKIPATHWAYS
import networkx as nx
import pybel
from pathme.constants import KEGG, REACTOME, WIKIPATHWAYS, PATHME_DIR
from pathme.normalize_names import normalize_graph_names
from pathme.pybel_utils import flatten_complex_nodes
from pybel import BELGraph, from_pickle, union
from pybel import BELGraph, from_pickle
from pybel.struct.utils import update_metadata
from pybel.constants import RELATION

from tqdm import tqdm

logger = logging.getLogger(__name__)


def set_resource(elements, database):
    for element, data in elements:
        if 'database' in data:
            data['database'].add(database)
        else:
            data['database'] = {database}


def set_graph_resource(graph, database):
    set_resource(graph.nodes(data=True), database)
    # set_resource(graph.edges(data=True), database)


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


def left_full_data_join(g, h) -> None:
    """Wrapper around PyBEL's left_full_join to merge node data.

    :param pybel.BELGraph g: A BEL graph
    :param pybel.BELGraph h: A BEL graph
    """
    for node, data in h.nodes(data=True):
        if node in g:
            if 'database' in data and 'database' in g.nodes[node]:
                g.nodes[node]['database'].update(data['database'])
        else:
            g.add_node(node, **data)

    g.add_edges_from(
        (u, v, key, data)
        for u, v, key, data in h.edges(keys=True, data=True)
        if u not in g or v not in g[u] or key not in g[u][v]
    )

    update_metadata(h, g)

    g.warnings.extend(h.warnings)


def union_data(graphs):
    """Wrapper around PyBEL's union to instate left_full_data_join function.

    Assumes iterator is longer than 2, but not infinite.

    :param iter[BELGraph] graphs: An iterator over BEL graphs. Can't be infinite.
    :return: A merged graph
    :rtype: BELGraph

    Example usage:
    """
    graphs = tuple(graphs)

    n_graphs = len(graphs)

    if n_graphs == 0:
        raise ValueError('no graphs given')

    if n_graphs == 1:
        return graphs[0]

    target = graphs[0].copy()

    for graph in graphs[1:]:
        left_full_data_join(target, graph)

    return target

def get_universe_graph(
        kegg_path: str,
        reactome_path: str,
        wikipathways_path: str,
        *,
        flatten: bool = True,
        normalize_names: bool = True,
) -> BELGraph:
    """Return universe graph."""
    universe_graphs = _iterate_universe_graphs(
        kegg_path, reactome_path, wikipathways_path,
        flatten=flatten,
        normalize_names=normalize_names
    )
    logger.info('Merging all into a hairball...')
    return union(universe_graphs)


def _iterate_universe_graphs(
        kegg_path: str,
        reactome_path: str,
        wikipathways_path: str,
        *,
        flatten: bool = True,
        normalize_names: bool = True,
) -> BELGraph:
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
            graph = from_pickle(os.path.join(kegg_path, file), check_version=False)

            if flatten:
                flatten_complex_nodes(graph)

            if normalize_names:
                normalize_graph_names(graph, KEGG)

            set_graph_resource(graph, 'kegg')


        elif file in reactome_pickles:
            graph = from_pickle(os.path.join(reactome_path, file), check_version=False)

            if flatten:
                flatten_complex_nodes(graph)

            if normalize_names:
                normalize_graph_names(graph, REACTOME)

            set_graph_resource(graph, 'reactome')


        elif file in wp_pickles:
            graph = from_pickle(os.path.join(wikipathways_path, file), check_version=False)

            if flatten:
                flatten_complex_nodes(graph)

            if normalize_names:
                normalize_grget_set_databaseaph_names(graph, WIKIPATHWAYS)

            set_graph_resource(graph, 'wikipathways')

        else:
            logger.warning(f'Unknown pickle file: {file}')
            continue

        universe_list.append(graph)

    logger.info('Merging all into a hairball...')

    return union_data(universe_list)


def munge_node_attribute(node, attribute='name'):
    if node.get(attribute) == None:
        return str(node)
    else:
        return node.get(attribute)


def to_gml(graph: pybel.BELGraph, path: str = PATHME_DIR) -> None:
    """Write this graph to GML  file using :func:`networkx.write_gml`.
    """
    rv = nx.MultiDiGraph()
        yield graph

    for node in graph:
        rv.add_node(munge_node_attribute(node, 'name'), namespace=str(node.get('namespace')),
                    function=node.get('function'))

    for u, v, key, edge_data in graph.edges(data=True, keys=True):
        rv.add_edge(
            munge_node_attribute(u),
            munge_node_attribute(v),
            interaction=str(edge_data[RELATION]),
            bel=str(edge_data),
            key=str(key),
        )

    nx.write_gml(rv, path)

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
