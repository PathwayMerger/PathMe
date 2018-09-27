# -*- coding: utf-8 -*-

"""Tests for converting WikiPathways."""

import os
import unittest

from pathme.wikipathways.utils import get_files_in_folder, merge_two_dicts, get_file_name_from_url
from tests.constants import WP_TEST_RESOURCES, WP22, WP2359


class TestUtils(unittest.TestCase):
    """Tests for utils."""

    def test_get_wikipathways_files(self):
        files = get_files_in_folder(WP_TEST_RESOURCES)

        self.assertEqual(len(files), 7)
        self.assertEqual(os.path.join(WP_TEST_RESOURCES, WP22), WP22)

    def test_merge_dict(self):
        dict_1 = {1: 'uno'}
        dict_2 = {2: 'dos'}
        merged_dict = merge_two_dicts(dict_1, dict_2)

        self.assertEqual(merged_dict, {1: 'uno', 2: 'dos'})

    def test_url(self):
        world = get_file_name_from_url('https://hello/world')

        self.assertEqual(world, 'world')
