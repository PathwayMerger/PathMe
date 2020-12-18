# -*- coding: utf-8 -*-

"""Command line interface for KEGG that can be run with ``python -m pathme.kegg``."""

import logging
import os
import time

import click
from bio2bel_chebi import Manager as ChebiManager
from bio2bel_hgnc import Manager as HgncManager
from pybel import from_pickle
from tqdm import tqdm

from .convert_to_bel import kegg_to_pickles
from .utils import download_kgml_files, get_kegg_pathway_ids
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


@main.command(help='Downloads KEGG files')
@click.option('-c', '--connection', help=f"Defaults to {KEGG_FILES}")
def download(connection):
    """Download KEGG KGML."""
    kegg_ids = get_kegg_pathway_ids(connection=connection)

    if click.confirm(
            'You are about to download KGML files from KEGG.\n'
            'Please make sure you have read KEGG license (see: https://www.kegg.jp/kegg/rest/).'
            ' These files cannot be distributed and their use must be exclusively with academic purposes.\n'
            'We (PathMe developers) are not responsible for the end use of this data.\n'
    ):
        click.echo('You have read and accepted the conditions stated above.\n')
        download_kgml_files(kegg_ids)


@main.command()
@click.option('-f', '--flatten', is_flag=True, default=False)
@click.option('-e', '--export-folder', default=KEGG_BEL, show_default=True)
@click.option('-v', '--debug', is_flag=True, default=False, help='Debug mode')
def bel(flatten, export_folder, debug):
    """Convert KEGG to BEL."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    logger.setLevel(logging.INFO)

    if debug:
        click.echo("Debug mode on")
        logger.setLevel(logging.DEBUG)

    t = time.time()

    logger.info('Initiating HGNC Manager')
    hgnc_manager = HgncManager()

    if not hgnc_manager.is_populated():
        click.echo('bio2bel_hgnc was not populated. Populating now.')
        hgnc_manager.populate()

    logger.info('Initiating ChEBI Manager')
    chebi_manager = ChebiManager()

    if not chebi_manager.is_populated():
        click.echo('bio2bel_chebi was not populated. Populating now.')
        chebi_manager.populate()

    if flatten:
        logger.info('Flattening mode activated')

    resource_paths = [
        path
        for path in get_paths_in_folder(KEGG_FILES)
    ]

    kegg_to_pickles(
        resource_files=resource_paths,
        resource_folder=KEGG_FILES,
        hgnc_manager=hgnc_manager,
        chebi_manager=chebi_manager,
        flatten=flatten,
        export_folder=export_folder,
    )

    logger.info('KEGG exported in %.2f seconds', time.time() - t)


@main.command()
@click.option('-e', '--export-folder', default=KEGG_BEL, show_default=True)
def summarize(export_folder):
    """Summarize the KEGG export."""
    click.echo('loading KEGG graphs')
    graphs = [
        from_pickle(os.path.join(export_folder, fname))
        for fname in tqdm(get_paths_in_folder(export_folder))
    ]

    if graphs:
        summarize_helper(graphs)
    else:
        click.echo("Please export KEGG to BEL first. Run 'python3 -m pathme kegg bel' ")


if __name__ == '__main__':
    main()
