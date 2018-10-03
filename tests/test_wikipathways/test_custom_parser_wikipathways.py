# -*- coding: utf-8 -*-

"""Tests for converting WikiPathways."""

import unittest

from pathme.wikipathways.json_rdf_parser import *


class TestWikiPathways(unittest.TestCase):
    """Tests for dealing with the WikiPathways sub-module."""

    def test_parse_id_uri(self):
        """Test parse identifier uri."""
        prefix, prefix_namespaces, namespace, identifier = parse_id_uri(
            'http://rdf.wikipathways.org/Pathway/WP22_r97775/Interaction/e3673')

        self.assertEqual(prefix, 'http://rdf.wikipathways.org')
        self.assertEqual(prefix_namespaces, 'Pathway/WP22_r97775')
        self.assertEqual(namespace, 'Interaction')
        self.assertEqual(identifier, 'e3673')

    def test_parse_namespace_uri(self):
        """Test parse namespace uri."""
        prefix, namespace, _ = parse_namespace_uri(
            'http://vocabularies.wikipathways.org/wp#participants')

        self.assertEqual(prefix, 'http://vocabularies.wikipathways.org')
        self.assertEqual(namespace, 'wp#participants')

    def test_match_entry_type(self):
        """Test match entry type."""
        entry_type = match_entry_type({'core#Collection', 'wp#Pathway'})
        self.assertEqual(entry_type, 'pathway_info')

        entry_type = match_entry_type({'wp#Interaction'})
        self.assertEqual(entry_type, 'complex')

    def test_match_attribute_label(self):
        """Test match attribute label."""
        attribute_label = match_attribute_label('source')
        self.assertEqual(attribute_label, 'source')

        attribute_label = match_attribute_label('wp#organismName')
        self.assertEqual(attribute_label, 'organismName')

        attribute_label = match_attribute_label('wp#bdbEnsembl')
        self.assertEqual(attribute_label, 'value_namespace')

    def test_match_entry(self):
        """Test match entry."""
        node_entry = {
            '@id': 'http://identifiers.org/ncbigene/5295',
            '@type': ['http://vocabularies.wikipathways.org/wp#Protein',
                      'http://vocabularies.wikipathways.org/wp#DataNode'
                      ],
            'http://purl.org/dc/elements/1.1/identifier': [{
                '@id': 'http://identifiers.org/ncbigene/5295'
            }],
            'http://purl.org/dc/elements/1.1/source': [{
                '@value': 'Entrez Gene'
            }],
            'http://purl.org/dc/terms/identifier': [{
                '@value': '5295'
            }],
            'http://purl.org/dc/terms/isPartOf': [{
                '@id': 'http://identifiers.org/wikipathways/WP22_r97775'
            }],
            'http://vocabularies.wikipathways.org/wp#bdbEnsembl': [{
                '@id': 'http://identifiers.org/ensembl/ENSG00000145675'
            }],
            'http://vocabularies.wikipathways.org/wp#bdbEntrezGene': [{
                '@id': 'http://identifiers.org/ncbigene/5295'
            }],
            'http://vocabularies.wikipathways.org/wp#bdbHgncSymbol': [{
                '@id': 'http://identifiers.org/hgnc.symbol/PIK3R1'
            }],
            'http://vocabularies.wikipathways.org/wp#bdbUniprot': [{
                '@id': 'http://identifiers.org/uniprot/E5RGI8'
            },
                {
                    '@id': 'http://identifiers.org/uniprot/H0YBC2'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/E5RJY0'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/E5RK66'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/H0YB27'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/E5RHI0'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/P27986'
                }
            ],
            'http://vocabularies.wikipathways.org/wp#isAbout': [{
                '@id': 'http://rdf.wikipathways.org/Pathway/WP22_r97775/DataNode/b5564'
            }],
            'http://www.w3.org/2000/01/rdf-schema#label': [{
                '@value': 'PIK3R1'
            }]
        }
        entry_id, entry_type = match_entry(node_entry)
        self.assertEqual(entry_type, 'nodes')
        self.assertEqual(entry_id, '5295')

        entry_pathway_version = {
            '@id': 'http://identifiers.org/wikipathways/WP22',
            'http://purl.org/pav/hasVersion': [{
                '@id': 'http://identifiers.org/wikipathways/WP22_r97775'
            }]
        }
        entry_id, entry_type = match_entry(entry_pathway_version)
        self.assertEqual(entry_type, 'pathway_info')
        self.assertEqual(entry_id, 'WP22')

        interaction_entry = {
            '@id': 'http://rdf.wikipathways.org/Pathway/WP22_r97775/WP/Interaction/a537f',
            '@type': ['http://vocabularies.wikipathways.org/wp#Interaction',
                      'http://vocabularies.wikipathways.org/wp#DirectedInteraction'
                      ],
            'http://purl.org/dc/terms/isPartOf': [{
                '@id': 'http://identifiers.org/wikipathways/WP22_r97775'
            }],
            'http://vocabularies.wikipathways.org/wp#isAbout': [{
                '@id': 'http://rdf.wikipathways.org/Pathway/WP22_r97775/Interaction/a537f'
            }],
            'http://vocabularies.wikipathways.org/wp#participants': [{
                '@id': 'http://identifiers.org/ncbigene/2885'
            },
                {
                    '@id': 'http://identifiers.org/ncbigene/1025'
                }
            ],
            'http://vocabularies.wikipathways.org/wp#source': [{
                '@id': 'http://identifiers.org/ncbigene/1025'
            }],
            'http://vocabularies.wikipathways.org/wp#target': [{
                '@id': 'http://identifiers.org/ncbigene/2885'
            }]
        }
        entry_id, entry_type = match_entry(interaction_entry)
        self.assertEqual(entry_type, 'interactions')
        self.assertEqual(entry_id, 'a537f')

    def test_match_attribute(self):
        """Test match attribute."""
        attribute_id_uri = 'http://purl.org/dc/terms/description'
        attribute_label = match_attribute(attribute_id_uri)
        self.assertEqual(attribute_label, 'description')

        attribute_id_uri = 'http://vocabularies.wikipathways.org/wp#participants'
        attribute_label = match_attribute(attribute_id_uri)
        self.assertEqual(attribute_label, 'participants')

        attribute_id_uri = 'http://vocabularies.wikipathways.org/wp#bdbEnsembl'
        attribute_label = match_attribute(attribute_id_uri)
        self.assertEqual(attribute_label, 'value_namespace')

    def set_entry_graph(self, entry, graph, parse):
        """Get the type and id from an entry and optionally parse the entry content importing it to a graph object."""
        entry_id, entry_type = match_entry(entry)
        if parse:
            parse_attributes(entry, entry_type, entry_id, graph)

        return graph, entry_id, entry_type

    def assert_get_value(self, entry_type, entry_id, graph, assert_values):
        """Assert a dict of attributes values (labels as keys and values) with the retrieval of the function get_node_attribute_value."""
        get_node_attribute_value = lambda attribute_label: get_entry_attribute_value(entry_type, entry_id,
                                                                                     attribute_label,
                                                                                     graph)
        for attribute, value in assert_values.items():
            returned_value = get_node_attribute_value(attribute)
            self.assertEqual(value, returned_value)

    def test_get_entry_attribute_value(self):
        """Test get entry attribute."""
        graph = generate_empty_pathway_graph()

        node_entry = {
            '@id': 'http://identifiers.org/ncbigene/5295',
            '@type': ['http://vocabularies.wikipathways.org/wp#Protein',
                      'http://vocabularies.wikipathways.org/wp#DataNode'
                      ],
            'http://purl.org/dc/elements/1.1/identifier': [{
                '@id': 'http://identifiers.org/ncbigene/5295'
            }],
            'http://purl.org/dc/elements/1.1/source': [{
                '@value': 'Entrez Gene'
            }],
            'http://purl.org/dc/terms/identifier': [{
                '@value': '5295'
            }],
            'http://purl.org/dc/terms/isPartOf': [{
                '@id': 'http://identifiers.org/wikipathways/WP22_r97775'
            }],
            'http://vocabularies.wikipathways.org/wp#bdbEnsembl': [{
                '@id': 'http://identifiers.org/ensembl/ENSG00000145675'
            }],
            'http://vocabularies.wikipathways.org/wp#bdbEntrezGene': [{
                '@id': 'http://identifiers.org/ncbigene/5295'
            }],
            'http://vocabularies.wikipathways.org/wp#bdbHgncSymbol': [{
                '@id': 'http://identifiers.org/hgnc.symbol/PIK3R1'
            }],
            'http://vocabularies.wikipathways.org/wp#bdbUniprot': [{
                '@id': 'http://identifiers.org/uniprot/E5RGI8'
            },
                {
                    '@id': 'http://identifiers.org/uniprot/H0YBC2'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/E5RJY0'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/E5RK66'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/H0YB27'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/E5RHI0'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/P27986'
                }
            ],
            'http://vocabularies.wikipathways.org/wp#isAbout': [{
                '@id': 'http://rdf.wikipathways.org/Pathway/WP22_r97775/DataNode/b5564'
            }],
            'http://www.w3.org/2000/01/rdf-schema#label': [{
                '@value': 'PIK3R1'
            }]
        }

        graph, entry_id, entry_type = self.set_entry_graph(node_entry, graph, True)

        attribute_values = {'identifier': '5295', 'source': 'Entrez Gene', 'isPartOf': 'WP22_r97775'}
        self.assert_get_value(entry_type, entry_id, graph, attribute_values)

    def test_set_entry_attribute(self):
        """Test get entry attribute."""
        graph = generate_empty_pathway_graph()

        node_entry = {
            '@id': 'http://identifiers.org/ncbigene/5295',
            '@type': ['http://vocabularies.wikipathways.org/wp#Protein',
                      'http://vocabularies.wikipathways.org/wp#DataNode'
                      ],
            'http://purl.org/dc/elements/1.1/identifier': [{
                '@id': 'http://identifiers.org/ncbigene/5295'
            }],
            'http://purl.org/dc/elements/1.1/source': [{
                '@value': 'Entrez Gene'
            }],
            'http://purl.org/dc/terms/identifier': [{
                '@value': '5295'
            }],
            'http://purl.org/dc/terms/isPartOf': [{
                '@id': 'http://identifiers.org/wikipathways/WP22_r97775'
            }],
            'http://vocabularies.wikipathways.org/wp#bdbEnsembl': [{
                '@id': 'http://identifiers.org/ensembl/ENSG00000145675'
            }],
            'http://vocabularies.wikipathways.org/wp#bdbEntrezGene': [{
                '@id': 'http://identifiers.org/ncbigene/5295'
            }],
            'http://vocabularies.wikipathways.org/wp#bdbHgncSymbol': [{
                '@id': 'http://identifiers.org/hgnc.symbol/PIK3R1'
            }],
            'http://vocabularies.wikipathways.org/wp#bdbUniprot': [{
                '@id': 'http://identifiers.org/uniprot/E5RGI8'
            },
                {
                    '@id': 'http://identifiers.org/uniprot/H0YBC2'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/E5RJY0'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/E5RK66'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/H0YB27'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/E5RHI0'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/P27986'
                }
            ],
            'http://vocabularies.wikipathways.org/wp#isAbout': [{
                '@id': 'http://rdf.wikipathways.org/Pathway/WP22_r97775/DataNode/b5564'
            }],
            'http://www.w3.org/2000/01/rdf-schema#label': [{
                '@value': 'PIK3R1'
            }]
        }

        graph = generate_empty_pathway_graph()
        entry_id, entry_type = match_entry(node_entry)

        set_node_attribute = lambda attribute_label, attribute_value: set_entry_attribute(entry_type, entry_id,
                                                                                          attribute_label,
                                                                                          attribute_value, graph)
        attribute_values = {'identifier': '5295', 'source': 'Entrez Gene', 'isPartOf': 'WP22_r97775'}
        for attribute, value in attribute_values.items():
            set_node_attribute(attribute, value)

        self.assert_get_value(entry_type, entry_id, graph, attribute_values)

    def test_set_interaction(self):
        """Test set interaction."""
        interaction_entry = {
            '@id': 'http://rdf.wikipathways.org/Pathway/WP22_r97775/WP/Interaction/e3673',
            '@type': ['http://vocabularies.wikipathways.org/wp#Interaction',
                      'http://vocabularies.wikipathways.org/wp#DirectedInteraction'
                      ],
            'http://purl.org/dc/terms/isPartOf': [{
                '@id': 'http://identifiers.org/wikipathways/WP22_r97775'
            }],
            'http://vocabularies.wikipathways.org/wp#isAbout': [{
                '@id': 'http://rdf.wikipathways.org/Pathway/WP22_r97775/Interaction/e3673'
            }],
            'http://vocabularies.wikipathways.org/wp#participants': [{
                '@id': 'http://identifiers.org/ncbigene/3716'
            },
                {
                    '@id': 'http://identifiers.org/ncbigene/1025'
                }
            ],
            'http://vocabularies.wikipathways.org/wp#source': [{
                '@id': 'http://identifiers.org/ncbigene/3716'
            }],
            'http://vocabularies.wikipathways.org/wp#target': [{
                '@id': 'http://identifiers.org/ncbigene/1025'
            }]
        }

        graph = generate_empty_pathway_graph()
        set_interaction(interaction_entry, graph)

        interactions = graph['interactions']
        interaction_set = {('3716', '1025', 'http://rdf.wikipathways.org/Pathway/WP22_r97775/Interaction/e3673')}
        self.assertEqual(interaction_set, interactions)

    def test_get_entry_type(self):
        """Test get entry type."""
        interaction_entry = {
            '@id': 'http://rdf.wikipathways.org/Pathway/WP22_r97775/WP/Interaction/e3673',
            '@type': ['http://vocabularies.wikipathways.org/wp#Interaction',
                      'http://vocabularies.wikipathways.org/wp#DirectedInteraction'
                      ],
            'http://purl.org/dc/terms/isPartOf': [{
                '@id': 'http://identifiers.org/wikipathways/WP22_r97775'
            }],
            'http://vocabularies.wikipathways.org/wp#isAbout': [{
                '@id': 'http://rdf.wikipathways.org/Pathway/WP22_r97775/Interaction/e3673'
            }],
            'http://vocabularies.wikipathways.org/wp#participants': [{
                '@id': 'http://identifiers.org/ncbigene/3716'
            },
                {
                    '@id': 'http://identifiers.org/ncbigene/1025'
                }
            ],
            'http://vocabularies.wikipathways.org/wp#source': [{
                '@id': 'http://identifiers.org/ncbigene/3716'
            }],
            'http://vocabularies.wikipathways.org/wp#target': [{
                '@id': 'http://identifiers.org/ncbigene/1025'
            }]
        }
        node_entry = {
            '@id': 'http://identifiers.org/ncbigene/5295',
            '@type': ['http://vocabularies.wikipathways.org/wp#Protein',
                      'http://vocabularies.wikipathways.org/wp#DataNode'
                      ],
            'http://purl.org/dc/elements/1.1/identifier': [{
                '@id': 'http://identifiers.org/ncbigene/5295'
            }],
            'http://purl.org/dc/elements/1.1/source': [{
                '@value': 'Entrez Gene'
            }],
            'http://purl.org/dc/terms/identifier': [{
                '@value': '5295'
            }],
            'http://purl.org/dc/terms/isPartOf': [{
                '@id': 'http://identifiers.org/wikipathways/WP22_r97775'
            }],
            'http://vocabularies.wikipathways.org/wp#bdbEnsembl': [{
                '@id': 'http://identifiers.org/ensembl/ENSG00000145675'
            }],
            'http://vocabularies.wikipathways.org/wp#bdbEntrezGene': [{
                '@id': 'http://identifiers.org/ncbigene/5295'
            }],
            'http://vocabularies.wikipathways.org/wp#bdbHgncSymbol': [{
                '@id': 'http://identifiers.org/hgnc.symbol/PIK3R1'
            }],
            'http://vocabularies.wikipathways.org/wp#bdbUniprot': [{
                '@id': 'http://identifiers.org/uniprot/E5RGI8'
            },
                {
                    '@id': 'http://identifiers.org/uniprot/H0YBC2'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/E5RJY0'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/E5RK66'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/H0YB27'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/E5RHI0'
                },
                {
                    '@id': 'http://identifiers.org/uniprot/P27986'
                }
            ],
            'http://vocabularies.wikipathways.org/wp#isAbout': [{
                '@id': 'http://rdf.wikipathways.org/Pathway/WP22_r97775/DataNode/b5564'
            }],
            'http://www.w3.org/2000/01/rdf-schema#label': [{
                '@value': 'PIK3R1'
            }]
        }

        interaction_types = interaction_entry['@type']
        node_types = node_entry['@type']

        interaction_type = get_entry_type(interaction_types)
        node_type = get_entry_type(node_types)

        self.assertEqual('interactions', interaction_type)
        self.assertEqual('nodes', node_type)
