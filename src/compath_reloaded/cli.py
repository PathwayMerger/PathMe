# -*- coding: utf-8 -*-

"""Command line interface."""

import logging
import os
import time

import click
from bio2bel_hgnc import Manager as HgncManager

from compath_reloaded.constants import (
    DATA_DIR,
    RDF_WIKIPATHWAYS,
    KEGG,
    WIKIPATHWAYS,
    REACTOME,
    RDF_REACTOME
)
from compath_reloaded.constants import DEFAULT_CACHE_CONNECTION
from compath_reloaded.kegg.utils import download_kgml_files, get_kegg_pathway_ids
from compath_reloaded.reactome.rdf_sparql import reactome_to_bel, get_reactome_statistics
from compath_reloaded.reactome.utils import untar_file
from compath_reloaded.utils import make_downloader, statistics_to_df
from compath_reloaded.wikipathways.rdf_sparql import get_wikipathways_statistics, wikipathways_to_pybel
from compath_reloaded.wikipathways.utils import (
    get_file_name_from_url,
    get_wikipathways_files,
    unzip_file
)

KEGG_DIR = os.path.join(DATA_DIR, KEGG)
REACTOME_DIR = os.path.join(DATA_DIR, REACTOME)
WIKIPATHWAYS_DIR = os.path.join(DATA_DIR, WIKIPATHWAYS)

log = logging.getLogger(__name__)

# Ensure data folders are created
os.makedirs(KEGG_DIR, exist_ok=True)
os.makedirs(REACTOME_DIR, exist_ok=True)
os.makedirs(WIKIPATHWAYS_DIR, exist_ok=True)


@click.group(help='ComPath Reloaded')
def main():
    """Run ComPath-Reloaded."""
    logging.basicConfig(format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")


"""KEGG"""


@main.group()
def kegg():
    """Manage KEGG."""


@kegg.command(help='Downloads KEGG files')
@click.option('-c', '--connection', help="Defaults to {}".format(os.path.join(DATA_DIR, KEGG)))
def download(connection):
    """Download KEGG KGML."""
    kegg_ids = get_kegg_pathway_ids(connection=connection)
    download_kgml_files(kegg_ids)


@kegg.command()
def populate():
    """Populate KEGG into PyBEL database."""
    # TODO
    pass


"""WikiPathways"""


@main.group()
def wikipathways():
    """Manage WikiPathways."""


@wikipathways.command(help='Downloads WikiPathways RDF files')
def download():
    """Download WikiPathways RDF."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    log.setLevel(logging.INFO)

    cached_file = os.path.join(WIKIPATHWAYS_DIR, get_file_name_from_url(RDF_WIKIPATHWAYS))
    make_downloader(RDF_WIKIPATHWAYS, cached_file, WIKIPATHWAYS, unzip_file)
    log.info('WikiPathways was downloaded')


@wikipathways.command()
@click.option('-c', '--connection', help="Defaults to {}".format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--verbose', is_flag=True)
@click.option('-x', '--only-canonical', default=True, help='Parse only canonical pathways')
@click.option('-e', '--export', default=False, help='Export to datasheet csv and xls')
def statistics(connection, verbose, only_canonical, export):
    """Generate statistics for a database."""
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    if verbose:
        log.setLevel(logging.DEBUG)

    # TODO: Allow for an optional parameter giving the folder of the files
    resource_folder = os.path.join(WIKIPATHWAYS_DIR, 'wp', 'Human')

    resource_files = get_wikipathways_files(resource_folder, connection, only_canonical)

    global_statistics, all_pathways_statistics = get_wikipathways_statistics(resource_files, resource_folder)

    df = statistics_to_df(all_pathways_statistics)

    df.to_excel(os.path.join(DATA_DIR, 'wikipathways_statistics.xlsx'))
    df.to_csv(os.path.join(DATA_DIR, 'wikipathways_statistics.csv'))


@wikipathways.command()
@click.option('-c', '--connection', help="Defaults to {}".format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--verbose', is_flag=True)
@click.option('-x', '--only-canonical', default=True, help='Parse only canonical pathways')
def load(connection, verbose, only_canonical):
    """Load WikiPathways and loads DB."""
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    if verbose:
        log.setLevel(logging.DEBUG)

    # TODO: Allow for an optional parameter giving the folder of the files
    resource_folder = os.path.join(WIKIPATHWAYS_DIR, 'wp', 'human')

    resource_files = get_wikipathways_files(resource_folder, connection, only_canonical)

    wikipathways_to_pybel(resource_files, resource_folder)


"""Reactome"""


@main.group()
def reactome():
    """Manage Reactome."""


@reactome.command(help='Downloads Reactome RDF files')
def download():
    """Download Reactome RDF."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    log.setLevel(logging.INFO)

    cached_file = os.path.join(REACTOME_DIR, get_file_name_from_url(RDF_REACTOME))
    make_downloader(RDF_REACTOME, cached_file, REACTOME, untar_file)
    log.info('Reactome was downloaded')


@reactome.command()
@click.option('-v', '--verbose', is_flag=True)
def load(verbose):
    """Populate pathways for a database."""
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    if verbose:
        log.setLevel(logging.INFO)

    t = time.time()

    log.info('Initiating HGNC Manager')
    hgnc_manager = HgncManager()

    resource_file = os.path.join(REACTOME_DIR, 'Homo_sapiens.owl')

    reactome_to_bel(resource_file, hgnc_manager)

    log.info('Reactome database laoded in %.2f seconds', time.time() - t)


@reactome.command()
@click.option('-c', '--connection', help="Defaults to {}".format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--verbose', is_flag=True)
@click.option('-x', '--only-canonical', default=True, help='Parse only canonical pathways')
@click.option('-e', '--export', default=False, help='Export to datasheet csv and xls')
def statistics(connection, verbose, only_canonical, export):
    """Generate statistics for a database."""
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    if verbose:
        log.setLevel(logging.DEBUG)

    log.info('Initiating HGNC Manager')
    hgnc_manager = HgncManager()

    resource_file = os.path.join(REACTOME_DIR, 'Homo_sapiens.owl')

    global_statistics, all_pathways_statistics = get_reactome_statistics(resource_file, hgnc_manager)

    print(global_statistics)

    if export:
        df = statistics_to_df(all_pathways_statistics)

        df.to_excel(os.path.join(DATA_DIR, 'wikipathways_statistics.xlsx'))
        df.to_csv(os.path.join(DATA_DIR, 'wikipathways_statistics.csv'))


if __name__ == '__main__':
    main()
