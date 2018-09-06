# -*- coding: utf-8 -*-

"""This module has utilities method for parsing, handling wikipathways RDF and data."""

import logging
import os
import tarfile

from ..constants import DATA_DIR, WIKIPATHWAYS, HGNC

WIKIPATHWAYS_DIR = os.path.join(DATA_DIR, WIKIPATHWAYS)

log = logging.getLogger(__name__)

"""Download utilities"""


def get_hgnc_node_info(gene):
    """Return HGNC identifier, symbol and namespace from HGNC entry.

    :param bio2bel_hgnc.manager.models.HGNC gene:
    :rtype: tuple[str,str,str]
    """
    return gene.identifier, gene.symbol, HGNC


def untar_file(file_path, export_folder):
    """Unzip file into a destination folder.

    :param str file_path: name of the file
    :param str export_folder: name of the file
    """
    tar_ref = tarfile.open(file_path, 'r:bz2')
    tar_ref.extractall(export_folder)
    tar_ref.close()
