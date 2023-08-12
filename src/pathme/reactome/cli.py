# -*- coding: utf-8 -*-

"""Command line interface."""

import logging
import os

import click
import time
from tqdm import tqdm

from pybel import from_pickle
from .rdf_sparql import get_reactome_statistics, reactome_to_bel
from .utils import untar_file
from ..constants import DATA_DIR, RDF_REACTOME, REACTOME_BEL, REACTOME_FILES
from ..export_utils import get_paths_in_folder
from ..utils import make_downloader, statistics_to_df, summarize_helper
from ..wikipathways.utils import get_file_name_from_url

__all__ = [
    'main',
]

logger = logging.getLogger(__name__)


@click.group()
def main():
    """Manage Reactome."""


@main.command(help='Downloads Reactome RDF files')
def download():
    """Download Reactome RDF."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    logger.setLevel(logging.INFO)

    logger.info('Downloading Reactome RDF file')

    cached_file = os.path.join(REACTOME_FILES, get_file_name_from_url(RDF_REACTOME))
    make_downloader(RDF_REACTOME, cached_file, REACTOME_FILES, untar_file)
    logger.info('Reactome was downloaded')


@main.command()
@click.option('-v', '--verbose', is_flag=True)
def bel(verbose):
    """Convert Reactome to BEL."""
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    t = time.time()

    resource_file = os.path.join(REACTOME_FILES, 'Homo_sapiens.owl')

    reactome_to_bel(resource_file)

    logger.info('Reactome exported in %.2f seconds', time.time() - t)


@main.command()
@click.option('-e', '--export-folder', default=REACTOME_BEL, show_default=True)
def summarize(export_folder):
    """Summarize the Reactome export."""
    click.echo('loading Reactome graphs')
    graphs = [
        from_pickle(os.path.join(export_folder, fname))
        for fname in tqdm(get_paths_in_folder(export_folder))
    ]

    if graphs:
        summarize_helper(graphs)
    else:
        click.echo("Please export Reactome to BEL first. Run 'python3 -m pathme reactome bel' ")


@main.command()
@click.option('-v', '--verbose', is_flag=True)
@click.option('-e', '--export', default=False, help='Export to datasheet csv and xls')
def statistics(verbose, export):
    """Generate statistics for a database."""
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    if verbose:
        logger.setLevel(logging.DEBUG)

    resource_file = os.path.join(REACTOME_FILES, 'Homo_sapiens.owl')

    global_statistics, all_pathways_statistics = get_reactome_statistics(resource_file=resource_file)

    if export:
        df = statistics_to_df(all_pathways_statistics)
        df.to_excel(os.path.join(DATA_DIR, 'reactome_statistics.xlsx'))
        df.to_csv(os.path.join(DATA_DIR, 'reactome_statistics.csv'))


if __name__ == '__main__':
    main()
