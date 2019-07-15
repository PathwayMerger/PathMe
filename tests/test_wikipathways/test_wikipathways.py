# -*- coding: utf-8 -*-

"""Tests for converting WikiPathways."""

import tempfile

from bio2bel.testing import TemporaryConnectionMixin
from bio2bel_hgnc import Manager as HgncManager
from bio2bel_kegg.manager import Manager
from pybel import BELGraph
from pybel_tools.summary.edge_summary import count_relations
from tests.constants import WP1871, WP22, WP2359, WP2799, WP706

from pathme.kegg.kegg_xml_parser import *
from pathme.utils import parse_rdf
from pathme.wikipathways.rdf_sparql import _get_interactions, _get_nodes, wikipathways_to_bel

dir_path = os.path.dirname(os.path.realpath(__file__))
resources_path = os.path.join(dir_path, 'resources')

hgnc_test_path = os.path.join(resources_path, 'hgnc_test.json')


class WikipathwaysTest(TemporaryConnectionMixin):
    """A test case with a populated HGNC databases for WikiPathways parser."""

    @classmethod
    def setUpClass(cls):
        """Create a temporary file."""
        cls.fd, cls.path = tempfile.mkstemp()
        cls.connection = 'sqlite:///' + cls.path

        """HGNC Manager"""
        cls.manager = Manager(cls.connection)

        cls.hgnc_manager = HgncManager(engine=cls.manager.engine, session=cls.manager.session)
        cls.hgnc_manager.populate(hgnc_file_path=hgnc_test_path, use_hcop=False)

        log.info('HGNC database loaded')


class TestWikipathways(WikipathwaysTest):
    """Tests for dealing with the WikiPathways sub-module."""

    def test_get_nodes(self):
        """Test get nodes from RDF pathways files."""
        wp22_rdf_graph = parse_rdf(WP22)
        wp706_rdf_graph = parse_rdf(WP706)
        wp1871_rdf_graph = parse_rdf(WP1871)
        wp2799_rdf_graph = parse_rdf(WP2799)

        nodes_wp22 = _get_nodes(wp22_rdf_graph)
        nodes_wp706 = _get_nodes(wp706_rdf_graph)
        nodes_wp1871 = _get_nodes(wp1871_rdf_graph)
        nodes_wp2799 = _get_nodes(wp2799_rdf_graph)

        self.assertEqual(len(nodes_wp22), 17)
        self.assertEqual(len(nodes_wp706), 186)
        self.assertEqual(len(nodes_wp1871), 115)
        self.assertEqual(len(nodes_wp2799), 141)

    def test_get_interactions(self):
        """Test get nodes from RDF pathways files."""
        wp22_rdf_graph = parse_rdf(WP22)
        wp706_rdf_graph = parse_rdf(WP706)
        wp1871_rdf_graph = parse_rdf(WP1871)
        wp2799_rdf_graph = parse_rdf(WP2799)

        nodes_wp22 = _get_interactions(wp22_rdf_graph)
        nodes_wp706 = _get_interactions(wp706_rdf_graph)
        nodes_wp1871 = _get_interactions(wp1871_rdf_graph)
        nodes_wp2799 = _get_interactions(wp2799_rdf_graph)

        self.assertEqual(len(nodes_wp22), 10)
        self.assertEqual(len(nodes_wp706), 44)
        self.assertEqual(len(nodes_wp1871), 51)
        self.assertEqual(len(nodes_wp2799), 28)

    def test_wp_association_bp(self):
        """Test connect bp with protein."""
        test_graph = wikipathways_to_bel(WP2359, self.hgnc_manager)

        self.assertEqual(type(test_graph), BELGraph, msg='Error with graph type')

        self.assertEqual(test_graph.summary_dict()['Number of Nodes'], 2)
        self.assertEqual(test_graph.summary_dict()['Number of Edges'], 1)
        self.assertEqual(count_relations(test_graph)['regulates'], 1)

    def test_wp_22_to_bel(self):
        """Test 22."""
        test_graph = wikipathways_to_bel(WP22, self.hgnc_manager)

        self.assertEqual(type(test_graph), BELGraph, msg='Error with graph type')

        self.assertEqual(test_graph.summary_dict()['Number of Nodes'], 11)
        self.assertEqual(test_graph.summary_dict()['Number of Edges'], 10)
