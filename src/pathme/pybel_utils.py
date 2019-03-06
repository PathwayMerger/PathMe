# -*- coding: utf-8 -*-

"""PyBEL generalized utils."""

import logging
from typing import Mapping, Iterable

from pybel.dsl import BaseEntity

logger = logging.getLogger(__name__)

from pybel import BELGraph

from pybel_tools.node_utils import list_abundance_cartesian_expansion, reaction_cartesian_expansion


def flatten_complex_nodes(graph: BELGraph) -> None:
    logger.info("Flat complexes and composites")
    list_abundance_cartesian_expansion(graph)
    reaction_cartesian_expansion(graph)


def multi_relabel(graph: BELGraph, mapping_dict: Mapping[BaseEntity, Iterable[BaseEntity]]) -> None:
    """Expand one victim to multiple survivor nodes, in place."""
    for victim, survivors in mapping_dict.items():
        for survivor in survivors:

            for u, _, k, d in graph.in_edges(victim, keys=True, data=True):
                graph.add_edge(u, survivor, key=k, **d)

            for _, v, k, d in graph.out_edges(victim, keys=True, data=True):
                graph.add_edge(survivor, v, key=k, **d)

    graph.remove_nodes_from(mapping_dict.keys())
