# -*- coding: utf-8 -*-

"""Tests for converting WikiPathways."""

import unittest

from pybel import BELGraph
from pybel_tools.summary.edge_summary import count_relations

from pathme.wikipathways.rdf_sparql import wikipathways_to_bel
from .constants import WP2359


class TestWikipathways(unittest.TestCase):
    """Tests for dealing with the WikiPathways sub-module."""

    def test_wikipathways_connect_bp(self):
        """Test connect bp with protein."""
        test_graph = wikipathways_to_bel(WP2359)

        self.assertEqual(type(test_graph), BELGraph, msg='Error with graph type')

        self.assertEqual(test_graph.summary_dict()['Number of Nodes'], 2)
        self.assertEqual(test_graph.summary_dict()['Number of Edges'], 1)
        self.assertEqual(count_relations(test_graph)['association'], 1)
