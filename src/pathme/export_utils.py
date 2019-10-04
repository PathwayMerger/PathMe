# -*- coding: utf-8 -*-

"""Export harmonized universe."""

import logging
import os
from typing import Iterable, List, Optional, Tuple

import networkx as nx
import pybel
from bio2bel_reactome import Manager as ReactomeManager
from bio2bel_reactome.models import Pathway
from pybel import BELGraph, from_pickle, union
from pybel.constants import ANNOTATIONS, RELATION
from pybel.struct import add_annotation_value
from pybel.struct.mutation import collapse_all_variants, collapse_to_genes
from pybel_tools.analysis.spia import bel_to_spia_matrices, spia_matrices_to_excel
from tqdm import tqdm

from pathme.constants import KEGG, PATHME_DIR, REACTOME, WIKIPATHWAYS
from pathme.normalize_names import normalize_graph_names
from pathme.pybel_utils import flatten_complex_nodes
from .constants import KEGG_BEL, REACTOME_BEL, WIKIPATHWAYS_BEL

logger = logging.getLogger(__name__)


def add_annotation_key(graph: BELGraph):
    """Add annotation key in data (in place operation)."""
    for u, v, k in graph.edges(keys=True):
        if ANNOTATIONS not in graph[u][v][k]:
            graph[u][v][k][ANNOTATIONS] = {}


def get_all_pickles(
    *,
    kegg_path: Optional[str] = None,
    reactome_path: Optional[str] = None,
    wikipathways_path: Optional[str] = None,
) -> Tuple[List[str], List[str], List[str]]:
    """Return a list with all pickle paths."""
    kegg_pickles = get_paths_in_folder(kegg_path or KEGG_BEL)
    if not kegg_pickles:
        logger.warning('No KEGG files found. Please create the BEL KEGG files')

    reactome_pickles = get_paths_in_folder(reactome_path or REACTOME_BEL)
    if not reactome_pickles:
        logger.warning('No Reactome files found. Please create the BEL Reactome files')

    wp_pickles = get_paths_in_folder(wikipathways_path or WIKIPATHWAYS_BEL)
    if not wp_pickles:
        logger.warning('No WikiPathways files found. Please create the BEL WikiPathways files')

    return kegg_pickles, reactome_pickles, wp_pickles


def get_universe_graph(
    *,
    kegg_path: Optional[str] = None,
    reactome_path: Optional[str] = None,
    wikipathways_path: Optional[str] = None,
    flatten: bool = True,
    normalize_names: bool = True,
) -> BELGraph:
    """Return universe graph."""
    universe_graphs = iterate_universe_graphs(
        kegg_path=kegg_path,
        reactome_path=reactome_path,
        wikipathways_path=wikipathways_path,
        flatten=flatten,
        normalize_names=normalize_names
    )
    # Just keep the graph and not the source
    universe_graphs = (graph for _, _, graph in universe_graphs)
    logger.info('Merging all into a hairball...')
    return union(universe_graphs)


def spia_export_helper(
    *,
    output: str,
    kegg_path: Optional[str] = None,
    reactome_path: Optional[str] = None,
    wikipathways_path: Optional[str] = None,
) -> None:
    """Export PathMe pickles to SPIA excel like file.

    :param output: output directory
    :param kegg_path: directory to KEGG pickles
    :param reactome_path: directory to Reactome pickles
    :param wikipathways_path: directory to WikiPathways pickles
    """
    kegg_pickles, reactome_pickles, wp_pickles = get_all_pickles(
        kegg_path=kegg_path,
        reactome_path=reactome_path,
        wikipathways_path=wikipathways_path,
    )

    paths = kegg_pickles + reactome_pickles + wp_pickles

    logger.info(f'A total of {len(paths)} will be exported')

    paths = tqdm(paths, desc='Exporting SPIA excel files')

    # Call Reactome manager and check that is populated
    reactome_manager = ReactomeManager()
    if not reactome_manager.is_populated():
        logger.warning('Reactome Manager is not populated')

    # Load each pickle and export it as excel file
    for path in paths:
        if not path.endswith('.pickle'):
            continue

        if path in kegg_pickles:
            pathway_graph = from_pickle(os.path.join(kegg_path, path))
            normalize_graph_names(pathway_graph, KEGG)

        elif path in reactome_pickles:
            # Load BELGraph
            pathway_graph = from_pickle(os.path.join(reactome_path, path))

            # Check if pathway has children to build the merge graph
            pathway_id = path.strip('.pickle')

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

        elif path in wp_pickles:
            pathway_graph = from_pickle(os.path.join(wikipathways_path, path))
            normalize_graph_names(pathway_graph, WIKIPATHWAYS)

        else:
            logger.warning(f'Unknown pickle file: {path}')
            continue

        # Explode complex nodes
        flatten_complex_nodes(pathway_graph)

        # Collapse nodes
        collapse_all_variants(pathway_graph)
        collapse_to_genes(pathway_graph)

        spia_matrices = bel_to_spia_matrices(pathway_graph)

        output_file = os.path.join(output, f"{path.strip('.pickle')}.xlsx")

        if os.path.isfile(output_file):
            continue

        # Export excel file representing the connectivity matrix of the BEL Graph
        spia_matrices_to_excel(spia_matrices, output_file)


def iterate_indra_statements(**kwargs) -> Iterable['indra.statements.Statement']:
    """Iterate over INDRA statements for the universe."""
    for _, _, graph in iterate_universe_graphs(**kwargs):
        yield from pybel.to_indra_statements(graph)


def iterate_universe_graphs(
    *,
    kegg_path: Optional[str] = None,
    reactome_path: Optional[str] = None,
    wikipathways_path: Optional[str] = None,
    flatten: bool = True,
    normalize_names: bool = True,
) -> Iterable[Tuple[str, str, BELGraph]]:
    """Return universe graph."""
    kegg_pickle_paths, reactome_pickle_paths, wp_pickle_paths = get_all_pickles(
        kegg_path=kegg_path,
        reactome_path=reactome_path,
        wikipathways_path=wikipathways_path,
    )

    n_paths = len(kegg_pickle_paths) + len(reactome_pickle_paths) + len(wp_pickle_paths)
    logger.info(f'{n_paths} graphs will be put in the universe')

    yield from _iterate_wp(wp_pickle_paths, wikipathways_path, flatten, normalize_names)
    yield from _iterate_kegg(kegg_pickle_paths, kegg_path, flatten, normalize_names)
    yield from _iterate_reactome(reactome_pickle_paths, reactome_path, flatten, normalize_names)


def _iterate_wp(wp_pickle_paths, wikipathways_path, flatten, normalize_names):
    for path in tqdm(wp_pickle_paths, desc=f'Loading WP pickles from {wikipathways_path}'):
        if not path.endswith('.pickle'):
            continue

        graph = from_pickle(os.path.join(wikipathways_path, path), check_version=False)

        if flatten:
            flatten_complex_nodes(graph)

        if normalize_names:
            normalize_graph_names(graph, WIKIPATHWAYS)

        _update_graph(graph, path, WIKIPATHWAYS)
        yield WIKIPATHWAYS, path, graph


def _iterate_kegg(kegg_pickle_paths, kegg_path, flatten, normalize_names):
    for path in tqdm(kegg_pickle_paths, desc=f'Loading KEGG pickles from {kegg_path}'):
        if not path.endswith('.pickle'):
            continue
        graph = from_pickle(os.path.join(kegg_path, path), check_version=False)

        if flatten:
            flatten_complex_nodes(graph)

        if normalize_names:
            normalize_graph_names(graph, KEGG)

        _update_graph(graph, path, KEGG)
        yield KEGG, path, graph


def _iterate_reactome(reactome_pickle_paths, reactome_path, flatten, normalize_names):
    for file in tqdm(reactome_pickle_paths, desc=f'Loading Reactome pickles from {reactome_path}'):
        if not file.endswith('.pickle'):
            continue

        graph = from_pickle(os.path.join(reactome_path, file), check_version=False)

        if flatten:
            flatten_complex_nodes(graph)

        if normalize_names:
            normalize_graph_names(graph, REACTOME)

        _update_graph(graph, file, REACTOME)
        yield REACTOME,path, graph


def _update_graph(graph, file, database):
    graph.annotation_list['database'] = {KEGG, REACTOME, WIKIPATHWAYS}
    add_annotation_key(graph)
    add_annotation_value(graph, 'database', database)
    graph.annotation_pattern['PathwayID'] = '.*'
    add_annotation_value(graph, 'PathwayID', file.strip(".pickle"))


def _munge_node_attribute(node, attribute='name'):
    """Munge node attribute."""
    if node.get(attribute) == None:
        return str(node)
    else:
        return node.get(attribute)


def to_gml(graph: pybel.BELGraph, path: str = PATHME_DIR) -> None:
    """Write this graph to GML  file using :func:`networkx.write_gml`."""
    rv = nx.MultiDiGraph()

    for node in graph:
        rv.add_node(
            _munge_node_attribute(node, 'name'),
            namespace=str(node.get('namespace')),
            function=node.get('function'),
        )

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
    """Transverse recursively the Reactome hierarchy and return all children for a given pathway."""
    if pathway.children:
        for child in pathway.children:
            yield child
            yield from yield_all_children(child)
