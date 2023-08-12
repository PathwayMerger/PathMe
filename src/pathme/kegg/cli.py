# -*- coding: utf-8 -*-

"""Command line interface for KEGG that can be run with ``python -m pathme.kegg``."""

import logging
import os
from typing import Optional

import click
import time
from tqdm import tqdm

import pybel
from pyobo.cli_utils import verbose_option
from .convert_to_bel import kegg_to_bel_cache
from .utils import ensure_kgml_files
from ..constants import KEGG_BEL, KEGG_FILES
from ..export_utils import get_paths_in_folder
from ..utils import summarize_helper

logger = logging.getLogger(__name__)

__all__ = [
    'main',
]


@click.group()
def main():
    """Manage KEGG."""


KEGG_PROMPT = (
    'You are about to download KGML files from KEGG.\n'
    'Please make sure you have read KEGG license (see: https://www.kegg.jp/kegg/rest/).'
    ' These files cannot be distributed and their use must be exclusively with academic purposes.\n'
    'We (PathMe developers) are not responsible for the end use of this data.\n'
)


@main.command(help='Downloads KEGG files')
@click.option('-c', '--connection', help='bio2bel manager connection')
@click.confirmation_option(prompt=KEGG_PROMPT)
@verbose_option
def download(connection):
    """Download KEGG KGML."""
    click.echo('You have read and accepted the conditions stated above.\n')
    ensure_kgml_files(connection=connection)


@main.command()
@click.option('-f', '--flatten', is_flag=True)
@click.option('-e', '--export-folder', default=KEGG_BEL, show_default=True)
@verbose_option
def bel(flatten: bool, export_folder: str):
    """Convert KEGG to BEL."""
    t = time.time()

    if flatten:
        logger.info('Flattening mode activated')

    resource_paths = [
        path
        for path in get_paths_in_folder(KEGG_FILES)
    ]

    kegg_to_bel_cache(
        resource_files=resource_paths,
        resource_folder=KEGG_FILES,
        flatten=flatten,
        export_folder=export_folder,
    )

    logger.info('KEGG exported in %.2f seconds', time.time() - t)


@main.command()
@click.option('-e', '--export-folder', default=KEGG_BEL, show_default=True)
@verbose_option
def summarize(export_folder):
    """Summarize the KEGG export."""
    click.echo('loading KEGG graphs')
    graphs = [
        pybel.load(os.path.join(export_folder, fname))
        for fname in tqdm(get_paths_in_folder(export_folder))
    ]

    if graphs:
        summarize_helper(graphs)
    else:
        click.echo("Please export KEGG to BEL first. Run 'python3 -m pathme kegg bel' ")


if __name__ == '__main__':
    main()
