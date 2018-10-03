# -*- coding: utf-8 -*-

"""Tests for converting WikiPathways."""

import unittest

from pybel_tools.summary.edge_summary import count_relations

from pathme.utils import parse_rdf
from pathme.wikipathways.rdf_sparql import wikipathways_to_bel, _get_nodes, _get_interactions
from pybel import BELGraph
from tests.constants import WP22, WP2359


class TestWikipathways(unittest.TestCase):
    """Tests for dealing with the WikiPathways sub-module."""

    def test_get_nodes(self):
        """Test get nodes from RDF pathways files."""
        wp22_rdf_graph = parse_rdf(WP22)
        wp706_rdf_graph = parse_rdf(WP22)
        wp1871_rdf_graph = parse_rdf(WP22)
        wp2799_rdf_graph = parse_rdf(WP22)

        nodes_wp22 = _get_nodes(wp22_rdf_graph)
        nodes_wp706 = _get_nodes(wp706_rdf_graph)
        nodes_wp1871 = _get_nodes(wp1871_rdf_graph)
        nodes_wp2799 = _get_nodes(wp2799_rdf_graph)

        self.assertEqual(len(nodes_wp22), 17)
        self.assertEqual(len(nodes_wp706), 181)
        self.assertEqual(len(nodes_wp1871), 102)
        self.assertEqual(len(nodes_wp2799), 139)

    def test_get_interactions(self):
        """Test get nodes from RDF pathways files."""
        wp22_rdf_graph = parse_rdf(WP22)
        wp706_rdf_graph = parse_rdf(WP22)
        wp1871_rdf_graph = parse_rdf(WP22)
        wp2799_rdf_graph = parse_rdf(WP22)

        nodes_wp22 = _get_interactions(wp22_rdf_graph)
        nodes_wp706 = _get_interactions(wp706_rdf_graph)
        nodes_wp1871 = _get_interactions(wp1871_rdf_graph)
        nodes_wp2799 = _get_interactions(wp2799_rdf_graph)

        self.assertEqual(len(nodes_wp22), 17)
        self.assertEqual(len(nodes_wp706), 44)
        self.assertEqual(len(nodes_wp1871), 102)
        self.assertEqual(len(nodes_wp2799), 38)

    def test_wp_association_bp(self):
        """Test connect bp with protein."""
        test_graph = wikipathways_to_bel(WP2359)

        self.assertEqual(type(test_graph), BELGraph, msg='Error with graph type')

        self.assertEqual(test_graph.summary_dict()['Number of Nodes'], 2)
        self.assertEqual(test_graph.summary_dict()['Number of Edges'], 1)
        self.assertEqual(count_relations(test_graph)['association'], 1)

    def test_wp_22_to_bel(self):
        """Test 22."""
        test_graph = wikipathways_to_bel(WP22)

        self.assertEqual(type(test_graph), BELGraph, msg='Error with graph type')

        self.assertEqual(test_graph.summary_dict()['Number of Nodes'], 11)
        self.assertEqual(test_graph.summary_dict()['Number of Edges'], 10)
