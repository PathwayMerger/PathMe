# -*- coding: utf-8 -*-

"""Command line interface for WikiPathways."""

import logging
import os

import click
import time
from tqdm import tqdm

from pybel import from_pickle
from .rdf_sparql import get_wp_statistics, wikipathways_to_pickles
from .utils import get_file_name_from_url, iterate_wikipathways_paths, unzip_file
from ..constants import (
    DATA_DIR, DEFAULT_CACHE_CONNECTION, RDF_WIKIPATHWAYS, WIKIPATHWAYS_BEL, WIKIPATHWAYS_DIR, WIKIPATHWAYS_FILES,
)
from ..export_utils import get_paths_in_folder
from ..utils import make_downloader, statistics_to_df, summarize_helper

logger = logging.getLogger(__name__)

__all__ = [
    'main',
]


@click.group()
def main():
    """Manage WikiPathways."""


@main.command(help='Downloads WikiPathways RDF files')
def download():
    """Download WikiPathways RDF."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    logger.setLevel(logging.INFO)

    cached_file = os.path.join(WIKIPATHWAYS_FILES, get_file_name_from_url(RDF_WIKIPATHWAYS))
    make_downloader(RDF_WIKIPATHWAYS, cached_file, WIKIPATHWAYS_FILES, unzip_file)
    logger.info('WikiPathways was downloaded')


@main.command()
@click.option('-c', '--connection', default=DEFAULT_CACHE_CONNECTION, show_default=True)
@click.option('-r', '--resource-folder')
@click.option('-d', '--export-folder', default=WIKIPATHWAYS_BEL)
@click.option('-v', '--debug', is_flag=True, default=False, help='Debug mode')
@click.option('-x', '--only-canonical', default=True, help='Parse only canonical pathways')
def bel(connection: str, resource_folder: str, export_folder: str, debug: bool, only_canonical: bool):
    """Convert WikiPathways to BEL."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")

    if debug:
        logger.setLevel(logging.DEBUG)

    t = time.time()

    if resource_folder is None:
        resource_folder = os.path.join(WIKIPATHWAYS_FILES, 'wp', 'Human')

    resource_files = iterate_wikipathways_paths(resource_folder, connection, only_canonical)

    os.makedirs(export_folder, exist_ok=True)
    wikipathways_to_pickles(resource_files, resource_folder, export_folder)

    logger.info('WikiPathways exported in %.2f seconds.', time.time() - t)


@main.command()
@click.option('-e', '--export-folder', default=WIKIPATHWAYS_BEL, show_default=True)
def summarize(export_folder):
    """Summarize the WikiPathways export."""
    click.echo('loading WikiPathways graphs')
    graphs = [
        from_pickle(os.path.join(export_folder, fname))
        for fname in tqdm(get_paths_in_folder(export_folder))
    ]

    if graphs:
        summarize_helper(graphs)
    else:
        click.echo("Please export WikiPathways to BEL first. Run 'python3 -m pathme wikipathways bel' ")


@main.command()
@click.option('-c', '--connection', default=DEFAULT_CACHE_CONNECTION, show_default=True)
@click.option('-v', '--verbose', is_flag=True)
@click.option('-x', '--only-canonical', default=True, help='Parse only canonical pathways')
@click.option('-e', '--export', default=False, help='Export to datasheet csv and xls')
def statistics(connection, verbose, only_canonical, export):
    """Generate statistics for a database."""
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    if verbose:
        logger.setLevel(logging.DEBUG)

    # TODO: Allow for an optional parameter giving the folder of the files
    resource_folder = os.path.join(WIKIPATHWAYS_DIR, 'wp', 'Human')

    resource_files = iterate_wikipathways_paths(resource_folder, connection, only_canonical)

    global_statistics, all_pathways_statistics = get_wp_statistics(resource_files, resource_folder)

    df = statistics_to_df(all_pathways_statistics)

    df.to_excel(os.path.join(DATA_DIR, 'wikipathways_statistics.xlsx'))
    df.to_csv(os.path.join(DATA_DIR, 'wikipathways_statistics.csv'))


if __name__ == '__main__':
    main()
