# -*- coding: utf-8 -*-

"""Command line interface."""

import logging
import time

import click
from bio2bel_chebi import Manager as ChebiManager
from bio2bel_hgnc import Manager as HgncManager
from tqdm import tqdm

from pathme.constants import *
from pathme.constants import DEFAULT_CACHE_CONNECTION
from pathme.kegg.convert_to_bel import kegg_to_pickles
from pathme.kegg.utils import download_kgml_files, get_kegg_pathway_ids
from pathme.reactome.rdf_sparql import get_reactome_statistics, reactome_to_bel
from pathme.reactome.utils import untar_file
from pathme.utils import CallCounted, get_files_in_folder, make_downloader, statistics_to_df, summarize_helper
from pathme.wikipathways.rdf_sparql import get_wp_statistics, wikipathways_to_pickles
from pathme.wikipathways.utils import get_file_name_from_url, get_wikipathways_files, unzip_file
from pybel import from_pickle

log = logging.getLogger(__name__)


@click.group(help='PathMe')
def main():
    """Run PathMe."""
    logging.basicConfig(format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")


"""KEGG"""


@main.group()
def kegg():
    """Manage KEGG."""


@kegg.command(help='Downloads KEGG files')
@click.option('-c', '--connection', help=f"Defaults to {KEGG_FILES}")
def download(connection):
    """Download KEGG KGML."""
    kegg_ids = get_kegg_pathway_ids(connection=connection)

    if click.confirm(
            'You are about to download KGML files from KEGG.\n'
            'Please make sure you have read KEGG license (see: https://www.kegg.jp/kegg/rest/).'
            ' These files cannot be distributed and their use must be exclusively academic.\n'
            'We (PathMe developers) are not responsible for the end use of this data.\n'
    ):
        click.echo('You have read and accepted the conditions stated above.\n')
        download_kgml_files(kegg_ids)


@kegg.command()
@click.option('-f', '--flatten', is_flag=True, default=False)
@click.option('-e', '--export-folder', default=KEGG_BEL, show_default=True)
def bel(flatten, export_folder):
    """Convert KEGG to BEL."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    log.setLevel(logging.INFO)

    t = time.time()

    log.info('Initiating HGNC Manager')
    hgnc_manager = HgncManager()

    log.info('Initiating ChEBI Manager')
    chebi_manager = ChebiManager()

    if flatten:
        log.info('Flattening mode activated')

    resource_files = [
        file
        for file in get_files_in_folder(KEGG_FILES)
    ]

    kegg_to_pickles(
        resource_files=resource_files,
        resource_folder=KEGG_FILES,
        hgnc_manager=hgnc_manager,
        chebi_manager=chebi_manager,
        flatten=flatten,
        export_folder=export_folder,
    )

    log.info('KEGG exported in %.2f seconds', time.time() - t)


@kegg.command()
@click.option('-e', '--export-folder', default=KEGG_BEL, show_default=True)
def summarize(export_folder):
    """Summarize the KEGG export."""
    click.echo('loading KEGG graphs')
    graphs = [
        from_pickle(os.path.join(export_folder, fname))
        for fname in tqdm(get_files_in_folder(export_folder))
    ]

    if graphs:
        summarize_helper(graphs)
    else:
        click.echo("Please export KEGG to BEL first. Run 'python3 -m pathme kegg bel' ")


"""WikiPathways"""


@main.group()
def wikipathways():
    """Manage WikiPathways."""


@wikipathways.command(help='Downloads WikiPathways RDF files')
def download():
    """Download WikiPathways RDF."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    log.setLevel(logging.INFO)

    cached_file = os.path.join(WIKIPATHWAYS_FILES, get_file_name_from_url(RDF_WIKIPATHWAYS))
    make_downloader(RDF_WIKIPATHWAYS, cached_file, WIKIPATHWAYS_FILES, unzip_file)
    log.info('WikiPathways was downloaded')


@wikipathways.command()
@click.option('-c', '--connection', help="Defaults to {}".format(DEFAULT_CACHE_CONNECTION))
@click.option('-d', '--debug', is_flag=True, default=False, help='Debug mode')
@click.option('-x', '--only-canonical', default=True, help='Parse only canonical pathways')
def bel(connection, debug, only_canonical):
    """Convert WikiPathways to BEL."""
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")

    if debug:
        log.setLevel(logging.DEBUG)

    logging.debug = CallCounted(logging.debug)

    log.info('Initiating HGNC Manager')
    hgnc_manager = HgncManager()

    t = time.time()

    # TODO: Allow for an optional parameter giving the folder of the files
    resource_folder = os.path.join(WIKIPATHWAYS_FILES, 'wp', 'Human')

    resource_files = get_wikipathways_files(resource_folder, connection, only_canonical)

    wikipathways_to_pickles(resource_files, resource_folder, hgnc_manager)

    log.info(
        'WikiPathways exported in %.2f seconds. A total of {} warnings regarding entities that could not be converted '
        'to standard identifiers were found.',
        time.time() - t, logging.debug.counter
    )


@wikipathways.command()
@click.option('-e', '--export-folder', default=WIKIPATHWAYS_BEL, show_default=True)
def summarize(export_folder):
    """Summarize the WikiPathways export."""
    click.echo('loading WikiPathways graphs')
    graphs = [
        from_pickle(os.path.join(export_folder, fname))
        for fname in tqdm(get_files_in_folder(export_folder))
    ]

    if graphs:
        summarize_helper(graphs)
    else:
        click.echo("Please export WikiPathways to BEL first. Run 'python3 -m pathme wikipathways bel' ")


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

    log.info('Initiating HGNC Manager')
    hgnc_manager = HgncManager()

    # TODO: Allow for an optional parameter giving the folder of the files
    resource_folder = os.path.join(WIKIPATHWAYS_DIR, 'wp', 'Human')

    resource_files = get_wikipathways_files(resource_folder, connection, only_canonical)

    global_statistics, all_pathways_statistics = get_wp_statistics(resource_files, resource_folder, hgnc_manager)

    df = statistics_to_df(all_pathways_statistics)

    df.to_excel(os.path.join(DATA_DIR, 'wikipathways_statistics.xlsx'))
    df.to_csv(os.path.join(DATA_DIR, 'wikipathways_statistics.csv'))


"""Reactome"""


@main.group()
def reactome():
    """Manage Reactome."""


@reactome.command(help='Downloads Reactome RDF files')
def download():
    """Download Reactome RDF."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    log.setLevel(logging.INFO)

    log.info('Downloading Reactome RDF file')

    cached_file = os.path.join(REACTOME_FILES, get_file_name_from_url(RDF_REACTOME))
    make_downloader(RDF_REACTOME, cached_file, REACTOME_FILES, untar_file)
    log.info('Reactome was downloaded')


@reactome.command()
@click.option('-v', '--verbose', is_flag=True)
def bel(verbose):
    """Convert Reactome to BEL."""
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    if verbose:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)

    t = time.time()

    log.info('Initiating HGNC Manager')
    hgnc_manager = HgncManager()

    resource_file = os.path.join(REACTOME_FILES, 'Homo_sapiens.owl')

    reactome_to_bel(resource_file, hgnc_manager)

    log.info('Reactome exported in %.2f seconds', time.time() - t)


@reactome.command()
@click.option('-e', '--export-folder', default=REACTOME_BEL, show_default=True)
def summarize(export_folder):
    """Summarize the Reactome export."""
    click.echo('loading Reactome graphs')
    graphs = [
        from_pickle(os.path.join(export_folder, fname))
        for fname in tqdm(get_files_in_folder(export_folder))
    ]

    if graphs:
        summarize_helper(graphs)
    else:
        click.echo("Please export Reactome to BEL first. Run 'python3 -m pathme reactome bel' ")


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

    resource_file = os.path.join(REACTOME_FILES, 'Homo_sapiens.owl')

    global_statistics, all_pathways_statistics = get_reactome_statistics(resource_file, hgnc_manager)

    if export:
        df = statistics_to_df(all_pathways_statistics)

        df.to_excel(os.path.join(DATA_DIR, 'reactome_statistics.xlsx'))
        df.to_csv(os.path.join(DATA_DIR, 'reactome_statistics.csv'))


if __name__ == '__main__':
    main()
