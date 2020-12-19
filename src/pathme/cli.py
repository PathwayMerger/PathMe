# -*- coding: utf-8 -*-

"""Command line interface."""

import logging
import os
import sys

import click
import pathme.kegg.cli
import pathme.reactome.cli
import pathme.wikipathways.cli

from .constants import CX_DIR, KEGG_BEL, PPI_DIR, REACTOME_BEL, SPIA_DIR, UNIVERSE_DIR, WIKIPATHWAYS_BEL, KEGG_FILES, \
    REACTOME_FILES, WIKIPATHWAYS_FILES
from .export_utils import export_helper, iterate_universe_graphs, generate_universe

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
def universe(
    kegg_path=KEGG_FILES,
    reactome_path=REACTOME_FILES,
    wikipathways_path=WIKIPATHWAYS_FILES,
    output=UNIVERSE_DIR,
    no_flatten=False,
    no_normalize_names=False,
    specie='Homo_sapiens'
):
    """Export harmonized PathMe universe."""
    generate_universe(kegg_path=kegg_path,
                      reactome_path=reactome_path,
                      wikipathways_path=wikipathways_path,
                      output=output,
                      no_flatten=no_flatten,
                      no_normalize_names=no_normalize_names,
                      specie=specie)


if __name__ == '__main__':
    main()
