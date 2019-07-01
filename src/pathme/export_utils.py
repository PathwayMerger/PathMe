# -*- coding: utf-8 -*-

"""Export harmonized universe."""

import logging
import os
from typing import Iterable, List, Tuple

import click
import networkx as nx
import pybel
from bio2bel_reactome import Manager as ReactomeManager
from bio2bel_reactome.models import Pathway
from pybel import BELGraph, union
from pybel import from_pickle
from pybel.constants import ANNOTATIONS, RELATION
from pybel.struct import add_annotation_value
from pybel.struct.mutation import collapse_all_variants, collapse_to_genes
from pybel_tools.analysis.spia import bel_to_spia_matrices, spia_matrices_to_excel
from tqdm import tqdm

from pathme.constants import KEGG, PATHME_DIR, REACTOME, WIKIPATHWAYS
from pathme.normalize_names import normalize_graph_names
from pathme.pybel_utils import flatten_complex_nodes

logger = logging.getLogger(__name__)


def add_annotation_key(graph):
    """Add annotation key in data (in place operation).

    :param pybel.BELGraph graph: BEL Graph
    """
    for u, v, k in graph.edges(keys=True):
        if ANNOTATIONS not in graph[u][v][k]:
            graph[u][v][k][ANNOTATIONS] = {}


def get_all_pickles(kegg_path: str, reactome_path: str, wikipathways_path: str) -> Tuple[List, List, List]:
    """Return a list with all pickle paths."""
    kegg_pickles = get_paths_in_folder(kegg_path)

    if not kegg_pickles:
        logger.warning('No KEGG files found. Please create the BEL KEGG files')

    reactome_pickles = get_paths_in_folder(reactome_path)

    if not reactome_pickles:
        logger.warning('No Reactome files found. Please create the BEL Reactome files')

    wp_pickles = get_paths_in_folder(wikipathways_path)

    if not wp_pickles:
        logger.warning('No WikiPathways files found. Please create the BEL WikiPathways files')

    return kegg_pickles, reactome_pickles, wp_pickles


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


def spia_export_helper(
        kegg_path: str,
        reactome_path: str,
        wikipathways_path: str,
        output: str
):
    """Export PathMe pickles to SPIA excel like file.

    :param kegg_path: directory to KEGG pickles
    :param reactome_path: directory to Reactome pickles
    :param wikipathways_path: directory to WikiPathways pickles
    :param output: output directory
    :param merge: merge equivalent pathways (default False)
    :return:
    """
    kegg_pickles, reactome_pickles, wp_pickles = get_all_pickles(kegg_path, reactome_path, wikipathways_path)

    all_pickles = kegg_pickles + reactome_pickles + wp_pickles

    click.echo(f'A total of {len(all_pickles)} will be exported')

    iterator = tqdm(all_pickles, desc='Exporting SPIA excel files')

    # Call Reactome manager and check that is populated
    reactome_manager = ReactomeManager()

    if not reactome_manager.is_populated():
        logger.warning('Reactome Manager is not populated')

    # Load each pickle and export it as excel file
    for file in iterator:
        if not file.endswith('.pickle'):
            continue

        if file in kegg_pickles:
            pathway_graph = from_pickle(os.path.join(kegg_path, file))
            normalize_graph_names(pathway_graph, KEGG)

        elif file in reactome_pickles:

            # Load BELGraph
            pathway_graph = from_pickle(os.path.join(reactome_path, file))

            # Check if pathway has children to build the merge graph
            pathway_id = file.strip('.pickle')

            # Look up in Bio2BEL Reactome
            pathway = reactome_manager.get_pathway_by_id(pathway_id)

            # Log if it is not present
            if not pathway:
                logger.warning(f'{pathway_id} not found in database')

            # Check if there are children and merge them on the fly
            for child in yield_all_children(pathway):

                child_file_path = os.path.join(reactome_path, f"{child.resource_id}.pickle")
                if not os.path.exists(child_file_path):
                    logger.warning(f'{child.resource_id} pickle does not exist')
                    continue

                # Load the pickle and union it
                child_graph = pybel.from_pickle(child_file_path)
                pathway_graph += child_graph

            # Normalize graph names
            normalize_graph_names(pathway_graph, REACTOME)

        elif file in wp_pickles:
            pathway_graph = from_pickle(os.path.join(wikipathways_path, file))
            normalize_graph_names(pathway_graph, WIKIPATHWAYS)

        else:
            logger.warning(f'Unknown pickle file: {file}')
            continue

        # Explode complex nodes
        flatten_complex_nodes(pathway_graph)

        # Collapse nodes
        collapse_all_variants(pathway_graph)
        collapse_to_genes(pathway_graph)

        spia_matrices = bel_to_spia_matrices(pathway_graph)

        output_file = os.path.join(output, f"{file.strip('.pickle')}.xlsx")

        if os.path.isfile(output_file):
            continue

        # Export excel file representing the connectivity matrix of the BEL Graph
        spia_matrices_to_excel(spia_matrices, output_file)


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

            graph.annotation_list['database'] = {KEGG, REACTOME, WIKIPATHWAYS}
            add_annotation_key(graph)
            add_annotation_value(graph, 'database', KEGG)

        elif file in reactome_pickles:
            graph = from_pickle(os.path.join(reactome_path, file), check_version=False)

            if flatten:
                flatten_complex_nodes(graph)

            if normalize_names:
                normalize_graph_names(graph, REACTOME)

            graph.annotation_list['database'] = {KEGG, REACTOME, WIKIPATHWAYS}
            add_annotation_key(graph)
            add_annotation_value(graph, 'database', REACTOME)


        elif file in wp_pickles:
            graph = from_pickle(os.path.join(wikipathways_path, file), check_version=False)

            if flatten:
                flatten_complex_nodes(graph)

            if normalize_names:
                normalize_graph_names(graph, WIKIPATHWAYS)

            graph.annotation_list['database'] = {KEGG, REACTOME, WIKIPATHWAYS}
            add_annotation_key(graph)
            add_annotation_value(graph, 'database', WIKIPATHWAYS)


        else:
            logger.warning(f'Unknown pickle file: {file}')
            continue

        graph.annotation_pattern['PathwayID'] = '.*'
        add_annotation_value(graph, 'PathwayID', file.strip(".pickle"))

        yield graph


def _munge_node_attribute(node, attribute='name'):
    """Munge node attribute."""
    if node.get(attribute) == None:
        return str(node)
    else:
        return node.get(attribute)


def to_gml(graph: pybel.BELGraph, path: str = PATHME_DIR) -> None:
    """Write this graph to GML  file using :func:`networkx.write_gml`.
    """
    rv = nx.MultiDiGraph()

    for node in graph:
        rv.add_node(_munge_node_attribute(node, 'name'), namespace=str(node.get('namespace')),
                    function=node.get('function'))

    for u, v, key, edge_data in graph.edges(data=True, keys=True):
        rv.add_edge(
            _munge_node_attribute(u),
            _munge_node_attribute(v),
            interaction=str(edge_data[RELATION]),
            bel=str(edge_data),
            key=str(key),
        )

    nx.write_gml(rv, path)


def get_paths_in_folder(directory: str) -> List[str]:
    """Return the files in a given folder.

    :param directory: folder path
    :return: file names in folder
    """
    return [
        path
        for path in os.listdir(directory)
        if os.path.isfile(os.path.join(directory, path))
    ]


def yield_all_children(pathway: Pathway) -> Iterable[Pathway]:
    """Transverse recursively reactome hierarchy and return all children for a given pathway."""
    if pathway.children:
        for child in pathway.children:
            yield child
            yield from yield_all_children(child)
