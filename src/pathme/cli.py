# -*- coding: utf-8 -*-

"""Command line interface."""

import logging
import os

import click
import networkx as nx
from pybel import to_pickle
from pybel.struct.mutation import collapse_all_variants, collapse_to_genes, remove_isolated_list_abundances
from pybel.struct.summary import count_functions

import pathme.kegg.cli
import pathme.reactome.cli
import pathme.wikipathways.cli
from pathme.constants import KEGG_BEL, REACTOME_BEL, SPIA_DIR, UNIVERSE_DIR, WIKIPATHWAYS_BEL
from pathme.export_utils import get_universe_graph, spia_export_helper

logger = logging.getLogger(__name__)

commands = dict(
    kegg=pathme.kegg.cli.main,
    wikipathways=pathme.wikipathways.cli.main,
    reactome=pathme.reactome.cli.main,
)


@click.group(help='PathMe', commands=commands)
def main():
    """Run PathMe."""
    logging.basicConfig(format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")


@main.group()
def export():
    """Export commands."""


@export.command()
@click.option('-k', '--kegg_path', help='KEGG BEL folder', default=KEGG_BEL, show_default=True)
@click.option('-r', '--reactome_path', help='Reactome BEL folder.', default=REACTOME_BEL, show_default=True)
@click.option('-w', '--wikipathways_path', help='WikiPathways BEL folder', default=WIKIPATHWAYS_BEL, show_default=True)
@click.option('-o', '--output', help='Output directory', default=SPIA_DIR, show_default=True)
def spia(kegg_path, reactome_path, wikipathways_path, output):
    """Export BEL Pickles to SPIA Excel."""
    click.echo(f'Results will be exported to {output}')
    spia_export_helper(kegg_path, reactome_path, wikipathways_path, output)


@export.command()
@click.option('-k', '--kegg-path', help='KEGG BEL folder', default=KEGG_BEL, show_default=True)
@click.option('-r', '--reactome-path', help='Reactome BEL folder.', default=REACTOME_BEL, show_default=True)
@click.option('-w', '--wikipathways-path', help='WikiPathways BEL folder', default=WIKIPATHWAYS_BEL, show_default=True)
@click.option('-o', '--output', help='Output directory', default=UNIVERSE_DIR, show_default=True)
@click.option('--no-flatten', is_flag=True, help='Do not flatten complex/reactions nodes')
@click.option('--no-normalize-names', is_flag=True, help='Do not normalize names')
def universe(kegg_path, reactome_path, wikipathways_path, output, no_flatten, no_normalize_names):
    """Export harmonized PathMe universe."""
    logging.basicConfig(level=logging.info, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    logger.setLevel(logging.INFO)

    flatten = not no_flatten
    normalize_names = not no_normalize_names

    if not flatten:
        click.echo('Complexes and Reactions will be not be flatten to single nodes')

    if not normalize_names:
        click.echo('Names will not be normalized to lower case')

    click.echo("Merging graphs to universe and harmonizing...(this might take a while)")

    # Not explode will flip the boolean coming from the cli
    universe_graph = get_universe_graph(
        kegg_path,
        reactome_path,
        wikipathways_path,
        flatten=flatten,
        normalize_names=normalize_names,
    )
    click.echo(f'Number of isolates after getting universe: {nx.number_of_isolates(universe_graph)}')

    # Remove isolated list abundances
    remove_isolated_list_abundances(universe_graph)

    if flatten:
        # TODO: Remove node list solo de Reactome
        click.echo(f'Number of isolates after flattening: {nx.number_of_isolates(universe_graph)}')

    click.echo("Merging variants and genes")
    collapse_all_variants(universe_graph)
    collapse_to_genes(universe_graph)
    click.echo(f'Number of isolates after collapsing variants and to genes: {nx.number_of_isolates(universe_graph)}')

    universe_graph.name = 'PathMe Universe'

    click.echo(f"Export BEL graph to: {os.path.join(output, 'pathme_universe_bel_graph.bel.pickle')}")
    click.echo(universe_graph.summary_str())
    click.echo(count_functions(universe_graph))

    to_pickle(universe_graph, os.path.join(output, "pathme_universe_bel_graph.bel.pickle"))


if __name__ == '__main__':
    main()
