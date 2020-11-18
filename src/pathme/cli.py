# -*- coding: utf-8 -*-

"""Command line interface."""

import logging
import os
import sys

import click
import networkx as nx
from pybel import to_pickle
from pybel.struct.mutation import collapse_all_variants, collapse_to_genes, remove_isolated_list_abundances
from pybel.struct.summary import count_functions

import pathme.kegg.cli
import pathme.reactome.cli
import pathme.wikipathways.cli
from pathme.kegg import download_kgml_files, get_kegg_pathway_ids
from .constants import CX_DIR, KEGG_BEL, PPI_DIR, REACTOME_BEL, SPIA_DIR, UNIVERSE_DIR, WIKIPATHWAYS_BEL, KEGG_FILES, \
    REACTOME_FILES, WIKIPATHWAYS_FILES
from .export_utils import export_helper, get_universe_graph, iterate_universe_graphs

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


kegg_path_option = click.option(
    '-k', '--kegg-path',
    help='KEGG BEL folder',
    default=KEGG_BEL,
    show_default=True,
)
reactome_path_option = click.option(
    '-r', '--reactome-path',
    help='Reactome BEL folder.',
    default=REACTOME_BEL,
    show_default=True,
)
wikipathways_path_option = click.option(
    '-w', '--wikipathways-path',
    help='WikiPathways BEL folder',
    default=WIKIPATHWAYS_BEL,
    show_default=True,
)

no_flatten = click.option('--no-flatten', is_flag=True, help='Do not flatten complex/reactions nodes')
no_normalize_names = click.option('--no-normalize-names', is_flag=True, help='Do not normalize names')


@export.command()
@kegg_path_option
@reactome_path_option
@wikipathways_path_option
@click.option('-o', '--output', help='Output directory', default=SPIA_DIR, show_default=True)
def spia(kegg_path, reactome_path, wikipathways_path, output):
    """Export BEL Pickles to SPIA Excel."""
    click.echo(f'Results will be exported to {output}')
    export_helper(
        kegg_path=kegg_path,
        reactome_path=reactome_path,
        wikipathways_path=wikipathways_path,
        output=output,
    )


@export.command()
@kegg_path_option
@reactome_path_option
@wikipathways_path_option
@click.option('-o', '--output', help='Output directory', default=PPI_DIR, show_default=True)
def ppi(kegg_path, reactome_path, wikipathways_path, output):
    """Export BEL Pickles to PPI-like tsv file."""
    click.echo(f'Results will be exported to {output}')
    export_helper(
        kegg_path=kegg_path,
        reactome_path=reactome_path,
        wikipathways_path=wikipathways_path,
        output=output,
        format='ppi',
    )


@export.command()
@kegg_path_option
@reactome_path_option
@wikipathways_path_option
@click.option('-o', '--output', help='Output directory', default=CX_DIR, show_default=True)
@no_flatten
@no_normalize_names
def cx(kegg_path, reactome_path, wikipathways_path, output, no_flatten, no_normalize_names):
    """Export BEL Pickles to CX."""
    try:
        from pybel_cx import to_cx_file
    except ImportError:
        click.secho('Could not import pybel_cx. Use pip install pybel-cx.')
        sys.exit(1)

    click.echo(f'Results will be exported to {output}')
    for source, path, graph in iterate_universe_graphs(
            kegg_path=kegg_path,
            reactome_path=reactome_path,
            wikipathways_path=wikipathways_path,
            flatten=(not no_flatten),
            normalize_names=(not no_normalize_names),
    ):
        with open(os.path.join(output, f"{path.strip('.pickle')}.cx.json"), 'w') as file:
            to_cx_file(graph, file)


@export.command()
@kegg_path_option
@reactome_path_option
@wikipathways_path_option
@click.option('-o', '--output', help='Output directory', default=UNIVERSE_DIR, show_default=True)
@no_flatten
@no_normalize_names
@click.option('-s', '--specie', help='Specie to geenerate universe.', default='Homo_sapiens', show_default=True)
def universe(kegg_path=KEGG_FILES,
             reactome_path=REACTOME_FILES,
             wikipathways_path=WIKIPATHWAYS_FILES,
             output=UNIVERSE_DIR,
             no_flatten=False,
             no_normalize_names=False,
             specie='Homo_sapiens'):
    """Export harmonized PathMe universe."""
    logging.basicConfig(level=logging.info, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    logger.setLevel(logging.INFO)

    flatten = not no_flatten
    normalize_names = not no_normalize_names
    specie = specie.replace(' ', '_').capitalize()

    if not flatten:
        click.echo('Complexes and Reactions will be not be flatten to single nodes')

    if not normalize_names:
        click.echo('Names will not be normalized to lower case')

    click.echo("Merging graphs to universe and harmonizing...(this might take a while)")

    # KEGG specie processing
    for _, dirs, _ in os.walk(os.path.abspath(kegg_path)):
        kegg_path = os.path.join(kegg_path, specie)
        if specie not in dirs:
            os.makedirs(kegg_path)
            kegg_ids = get_kegg_pathway_ids(connection=kegg_path, populate=True)
            if click.confirm(
                    'You are about to download KGML files from KEGG.\n'
                    'Please make sure you have read KEGG license (see: https://www.kegg.jp/kegg/rest/).'
                    ' These files cannot be distributed and their use must be exclusively with academic purposes.\n'
                    'We (PathMe developers) are not responsible for the end use of this data.\n'
            ):
                click.echo('You have read and accepted the conditions stated above.\n')
                download_kgml_files(kegg_ids, path=os.path.join(kegg_path, specie))

        break

    # Reactome specie processing
    for _, _, files in os.walk(os.path.abspath(reactome_path)):
        specie_file = '.'.join([specie, '.owl'])
        if '.'.join([specie, '.owl']) in files:
            reactome_path = os.path.join(reactome_path, specie_file)
        else:
            click.warning('Specie not found in the populated Reactome resources.')
        break

    # WikiPathways specie processing
    for _, dirs, _ in os.walk(os.path.abspath(wikipathways_path)):
        if specie in dirs:
            wikipathways_path = os.path.join(wikipathways_path, specie)
        else:
            click.warning('Specie not found in the populated WikiPathways resources.')
        break

    # Not explode will flip the boolean coming from the cli
    universe_graph = get_universe_graph(
        kegg_path=kegg_path,
        reactome_path=reactome_path,
        wikipathways_path=wikipathways_path,
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

    file_name = os.path.join(output, '_'.join([specie, 'pathme_universe_bel_graph.bel.pickle']))
    click.echo(f"Export BEL graph to: {file_name}")
    click.echo(universe_graph.summary_str())
    click.echo(count_functions(universe_graph))

    to_pickle(universe_graph, file_name)


def generate_universe(kegg_path=KEGG_FILES,
                      reactome_path=REACTOME_FILES,
                      wikipathways_path=WIKIPATHWAYS_FILES,
                      output=UNIVERSE_DIR,
                      no_flatten=False,
                      no_normalize_names=False,
                      specie='Homo_sapiens'):
    """Export harmonized PathMe universe."""

    flatten = not no_flatten
    normalize_names = not no_normalize_names
    specie = specie.replace(' ', '_').capitalize()

    if not flatten:
        click.secho('Complexes and Reactions will be not be flatten to single nodes')

    if not normalize_names:
        click.secho('Names will not be normalized to lower case')

    # KEGG specie processing
    for _, dirs, _ in os.walk(os.path.abspath(kegg_path)):
        kegg_path = os.path.join(kegg_path, specie)
        if specie not in dirs:
            kegg_ids = get_kegg_pathway_ids(populate=True)
            click.secho(
                'You are about to download KGML files from KEGG.\n'
                'Please make sure you have read KEGG license (see: https://www.kegg.jp/kegg/rest/).'
                ' These files cannot be distributed and their use must be exclusively with academic purposes.\n'
                'We (PathMe developers) are not responsible for the end use of this data.\n'
            )
            os.makedirs(kegg_path)
            download_kgml_files(kegg_ids, path=os.path.join(kegg_path, specie))
        break

    # Reactome specie processing
    for _, _, files in os.walk(os.path.abspath(reactome_path)):
        specie_file = '.'.join([specie, 'owl'])
        if specie_file in files:
            reactome_path = os.path.join(reactome_path, specie_file)
        else:
            click.secho('Specie not found in the populated Reactome resources.')
        break

    # WikiPathways specie processing
    for _, dirs, _ in os.walk(os.path.abspath(wikipathways_path)):
        if specie in dirs:
            wikipathways_path = os.path.join(wikipathways_path, specie)
        else:
            click.secho('Specie not found in the populated WikiPathways resources.')
        break

    click.secho("Merging graphs to universe and harmonizing...(this might take a while)")

    # Not explode will flip the boolean coming from the cli
    universe_graph = get_universe_graph(
        kegg_path=kegg_path,
        reactome_path=reactome_path,
        wikipathways_path=wikipathways_path,
        flatten=flatten,
        normalize_names=normalize_names,
    )
    click.secho(f'Number of isolates after getting universe: {nx.number_of_isolates(universe_graph)}')

    # Remove isolated list abundances
    remove_isolated_list_abundances(universe_graph)

    if flatten:
        # TODO: Remove node list solo de Reactome
        click.secho(f'Number of isolates after flattening: {nx.number_of_isolates(universe_graph)}')

    click.secho("Merging variants and genes")
    collapse_all_variants(universe_graph)
    collapse_to_genes(universe_graph)
    click.secho(f'Number of isolates after collapsing variants and to genes: {nx.number_of_isolates(universe_graph)}')

    universe_graph.name = 'PathMe Universe'

    file_name = os.path.join(output, '_'.join([specie, 'pathme_universe.pickle']))
    click.secho(f"Export BEL graph to: {file_name}")
    click.secho(universe_graph.summary_str())
    click.secho(count_functions(universe_graph))

    to_pickle(universe_graph, file_name)


if __name__ == '__main__':
    main()
