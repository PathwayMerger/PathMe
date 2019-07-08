# -*- coding: utf-8 -*-

"""Command line interface."""

import logging
import time

import click
import networkx as nx
from bio2bel_chebi import Manager as ChebiManager
from bio2bel_hgnc import Manager as HgncManager
from pybel import from_pickle, to_pickle
from pybel.struct.mutation import collapse_all_variants, collapse_to_genes, remove_isolated_list_abundances
from pybel.struct.summary import count_functions
from tqdm import tqdm

from pathme.constants import *
from pathme.export_utils import spia_export_helper, get_paths_in_folder, get_universe_graph
from pathme.kegg.convert_to_bel import kegg_to_pickles
from pathme.kegg.utils import download_kgml_files, get_kegg_pathway_ids
from pathme.reactome.rdf_sparql import get_reactome_statistics, reactome_to_bel
from pathme.reactome.utils import untar_file
from pathme.utils import CallCounted, make_downloader, statistics_to_df, summarize_helper
from pathme.wikipathways.rdf_sparql import get_wp_statistics, wikipathways_to_pickles
from pathme.wikipathways.utils import get_file_name_from_url, iterate_wikipathways_paths, unzip_file

logger = logging.getLogger(__name__)


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
            ' These files cannot be distributed and their use must be exclusively with academic purposes.\n'
            'We (PathMe developers) are not responsible for the end use of this data.\n'
    ):
        click.echo('You have read and accepted the conditions stated above.\n')
        download_kgml_files(kegg_ids)


@kegg.command()
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
        raise EnvironmentError('Your HGNC database is not populated. Please run python3 -m bio2bel_hgnc populate')

    logger.info('Initiating ChEBI Manager')
    chebi_manager = ChebiManager()

    if not chebi_manager.is_populated():
        raise EnvironmentError('Your CHEBI database is not populated. Please run python3 -m bio2bel_chebi populate')

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


@kegg.command()
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


"""WikiPathways"""


@main.group()
def wikipathways():
    """Manage WikiPathways."""


@wikipathways.command(help='Downloads WikiPathways RDF files')
def download():
    """Download WikiPathways RDF."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    logger.setLevel(logging.INFO)

    cached_file = os.path.join(WIKIPATHWAYS_FILES, get_file_name_from_url(RDF_WIKIPATHWAYS))
    make_downloader(RDF_WIKIPATHWAYS, cached_file, WIKIPATHWAYS_FILES, unzip_file)
    logger.info('WikiPathways was downloaded')


@wikipathways.command()
@click.option('-c', '--connection', help="Defaults to {}".format(DEFAULT_CACHE_CONNECTION))
@click.option('-r', '--resource-folder')
@click.option('-d', '--export-folder', default=WIKIPATHWAYS_BEL)
@click.option('-v', '--debug', is_flag=True, default=False, help='Debug mode')
@click.option('-x', '--only-canonical', default=True, help='Parse only canonical pathways')
def bel(connection: str, resource_folder: str, export_folder: str, debug: bool, only_canonical: bool):
    """Convert WikiPathways to BEL."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")

    if debug:
        logger.setLevel(logging.DEBUG)

    logging.debug = CallCounted(logging.debug)

    logger.info('Initiating HGNC Manager')
    hgnc_manager = HgncManager()

    if not hgnc_manager.is_populated():
        raise EnvironmentError('Your HGNC database is not populated. Please run python3 -m bio2bel_hgnc populate')

    t = time.time()

    if resource_folder is None:
        resource_folder = os.path.join(WIKIPATHWAYS_FILES, 'wp', 'Human')

    resource_files = iterate_wikipathways_paths(resource_folder, connection, only_canonical)

    os.makedirs(export_folder, exist_ok=True)
    wikipathways_to_pickles(resource_files, resource_folder, hgnc_manager, export_folder)

    logger.info(
        'WikiPathways exported in %.2f seconds. A total of %d warnings regarding entities that could not be converted '
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
        for fname in tqdm(get_paths_in_folder(export_folder))
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
        logger.setLevel(logging.DEBUG)

    logger.info('Initiating HGNC Manager')
    hgnc_manager = HgncManager()

    # TODO: Allow for an optional parameter giving the folder of the files
    resource_folder = os.path.join(WIKIPATHWAYS_DIR, 'wp', 'Human')

    resource_files = iterate_wikipathways_paths(resource_folder, connection, only_canonical)

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
    logger.setLevel(logging.INFO)

    logger.info('Downloading Reactome RDF file')

    cached_file = os.path.join(REACTOME_FILES, get_file_name_from_url(RDF_REACTOME))
    make_downloader(RDF_REACTOME, cached_file, REACTOME_FILES, untar_file)
    logger.info('Reactome was downloaded')


@reactome.command()
@click.option('-v', '--verbose', is_flag=True)
def bel(verbose):
    """Convert Reactome to BEL."""
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    t = time.time()

    logger.info('Initiating HGNC Manager')
    hgnc_manager = HgncManager()
    chebi_manager = ChebiManager()

    if not hgnc_manager.is_populated():
        raise EnvironmentError('Your HGNC database is not populated. Please run python3 -m bio2bel_hgnc populate')

    resource_file = os.path.join(REACTOME_FILES, 'Homo_sapiens.owl')

    reactome_to_bel(resource_file, hgnc_manager, chebi_manager)

    logger.info('Reactome exported in %.2f seconds', time.time() - t)


@reactome.command()
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


@reactome.command()
@click.option('-c', '--connection', help="Defaults to {}".format(DEFAULT_CACHE_CONNECTION))
@click.option('-v', '--verbose', is_flag=True)
@click.option('-x', '--only-canonical', default=True, help='Parse only canonical pathways')
@click.option('-e', '--export', default=False, help='Export to datasheet csv and xls')
def statistics(connection, verbose, only_canonical, export):
    """Generate statistics for a database."""
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    if verbose:
        logger.setLevel(logging.DEBUG)

    logger.info('Initiating HGNC Manager')
    hgnc_manager = HgncManager()
    chebi_manager = ChebiManager()

    resource_file = os.path.join(REACTOME_FILES, 'Homo_sapiens.owl')

    global_statistics, all_pathways_statistics = get_reactome_statistics(resource_file, hgnc_manager, chebi_manager)

    if export:
        df = statistics_to_df(all_pathways_statistics)

        df.to_excel(os.path.join(DATA_DIR, 'reactome_statistics.xlsx'))
        df.to_csv(os.path.join(DATA_DIR, 'reactome_statistics.csv'))


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
