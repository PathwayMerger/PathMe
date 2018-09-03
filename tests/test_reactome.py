# -*- coding: utf-8 -*-

"""Tests for converting Reactome."""

import os
import unittest

from compath_reloaded.cli import REACTOME_DIR, RDF_REACTOME, REACTOME
from compath_reloaded.reactome.rdf_sparql import _get_all_entry_types
from compath_reloaded.reactome.utils import untar_file
from compath_reloaded.utils import make_downloader
from compath_reloaded.utils import parse_rdf
from compath_reloaded.wikipathways.utils import get_file_name_from_url


class TestReactome(unittest.TestCase):
    """Tests for dealing with the Reactome sub-module."""

    def setUp(self):
        """Parse two examples files."""
        cached_file = os.path.join(REACTOME_DIR, get_file_name_from_url(RDF_REACTOME))
        make_downloader(RDF_REACTOME, cached_file, REACTOME, untar_file)

        resource_file = os.path.join(REACTOME_DIR, 'Homo_sapiens.owl')

        self.rdf_graph = parse_rdf(resource_file, fmt='xml')

    def test_reactome_types(self):
        types = _get_all_entry_types(self.rdf_graph)

        self.assertIn('Complex', types)
        self.assertIn('Dna', types)
        self.assertIn('Rna', types)
        self.assertIn('Protein', types)
        self.assertIn('SmallMolecule', types)
        self.assertIn('PhysicalEntity', types)
        self.assertIn('Dna', types)
