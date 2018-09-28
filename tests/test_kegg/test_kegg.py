# -*- coding: utf-8 -*-

"""Tests for converting KEGG."""

from pybel.dsl.nodes import abundance, bioprocess, composite_abundance, protein
from pybel.struct.summary.node_summary import count_functions
from pybel_tools.summary.edge_summary import count_relations

from pathme.constants import *
from pathme.kegg.convert_to_bel import xml_complexes_to_bel, xml_entities_to_bel
from pathme.kegg.kegg_xml_parser import _process_kegg_api_get_entity, get_all_reactions, get_all_relationships, \
    get_complex_components, get_entity_nodes, get_reaction_pathway_edges

from tests.constants import KeggTest


class TestKegg(KeggTest):
    """Tests for dealing with the KEGG sub-module."""

    def test_get_entities_from_xml(self):
        """Test entity creation."""
        notch_genes, notch_compounds, notch_maps, notch_orthologs = get_entity_nodes(
            self.notch_tree,
            self.hgnc_manager,
            self.chebi_manager
        )
        glycolysis_genes, glycolysis_compounds, glycolysis_maps, glycolysis_orthologs = get_entity_nodes(
            self.glycolysis_tree, self.hgnc_manager, self.chebi_manager)

        self.assertEqual(len(notch_genes), 22)
        self.assertEqual(len(notch_compounds), 0)
        self.assertEqual(len(notch_maps), 2)
        self.assertEqual(len(notch_orthologs), 2)

        self.assertEqual(len(glycolysis_genes), 35)
        self.assertEqual(len(glycolysis_compounds), 31)
        self.assertEqual(len(glycolysis_maps), 7)
        self.assertEqual(len(glycolysis_orthologs), 27)

        self.assertEqual(notch_genes['3'], [
            {KEGG_ID: 'hsa:5663',
             KEGG_TYPE: 'gene',
             HGNC_SYMBOL: 'PSEN1',
             HGNC: '9508'},
            {KEGG_ID: 'hsa:5664',
             KEGG_TYPE: 'gene',
             HGNC_SYMBOL: 'PSEN2',
             HGNC: '9509'}
        ])
        self.assertEqual(glycolysis_genes['45'], [
            {KEGG_ID: 'hsa:10327',
             KEGG_TYPE: 'gene',
             HGNC_SYMBOL: 'AKR1A1',
             HGNC: '380'}
        ])
        self.assertEqual(glycolysis_compounds['83'], [
            {KEGG_ID: 'cpd:C00031', CHEBI: '4167', CHEBI_NAME: 'D-glucopyranose', PUBCHEM: '3333',
             KEGG_TYPE: 'compound'}
        ])
        self.assertEqual(notch_maps['4'], [
            {KEGG_ID: 'path:hsa04010', 'map_name': 'MAPK signaling pathway'}
        ])
        self.assertEqual(glycolysis_maps['52'], [
            {KEGG_ID: 'path:hsa00620', 'map_name': 'Pyruvate metabolism'}
        ])
        self.assertEqual(notch_orthologs['7'], [
            {KEGG_ID: 'ko:K04497', KEGG_TYPE: 'ortholog'}
        ])
        self.assertEqual(glycolysis_orthologs['73'], [
            {KEGG_ID: 'ko:K01085', KEGG_TYPE: 'ortholog'},
            {KEGG_ID: 'ko:K20866', KEGG_TYPE: 'ortholog'}
        ])

    def test_get_complex_components(self):
        """Test creation of complexes."""
        notch_genes, notch_compounds, notch_maps, notch_orthologs = get_entity_nodes(
            self.notch_tree,
            self.hgnc_manager,
            self.chebi_manager
        )
        complex_ids, flattened_complexes = get_complex_components(self.notch_tree, notch_genes, flattened=True)

        self.assertEqual(len(complex_ids), 4)
        self.assertEqual(len(flattened_complexes), 4)
        self.assertEqual(complex_ids['29'], ['5', '8'])
        self.assertEqual(flattened_complexes['29'], [
            {KEGG_ID: 'hsa:3065',
             KEGG_TYPE: 'gene',
             HGNC_SYMBOL: 'HDAC1',
             HGNC: '4852'},
            {KEGG_ID: 'hsa:3066',
             KEGG_TYPE: 'gene',
             HGNC_SYMBOL: 'HDAC2',
             HGNC: '4853'},
            {KEGG_ID: 'hsa:1487',
             KEGG_TYPE: 'gene',
             HGNC_SYMBOL: 'CTBP1',
             HGNC: '2494'},
            {KEGG_ID: 'hsa:1488',
             KEGG_TYPE: 'gene',
             HGNC_SYMBOL: 'CTBP2',
             HGNC: '2495'}
        ])

    def test_get_all_relationships(self):
        """Test relationships."""
        notch_relations = get_all_relationships(self.notch_tree)
        glycolysis_relations = get_all_relationships(self.glycolysis_tree)

        notch_relation, glycolysis_relation = None, None

        for source, target, relation in notch_relations:
            if source == '21' and target == '22':
                notch_relation = relation

        self.assertEqual(len(notch_relations), 16)
        self.assertEqual(notch_relation, 'inhibition')

        for source, target, relation in glycolysis_relations:
            if source == '51' and target == '42':
                glycolysis_relation = relation

        self.assertEqual(glycolysis_relation, 'binding/association')

    def test_get_compound(self):
        """Test compound info."""
        compound_name = 'cpd:C01172'
        compound_info = _process_kegg_api_get_entity(compound_name, 'compound', self.hgnc_manager, self.chebi_manager)

        self.assertEqual(compound_info, {
            KEGG_ID: 'cpd:C01172',
            KEGG_TYPE: 'compound',
            CHEBI: '17719',
            CHEBI_NAME: 'beta-D-glucose 6-phosphate',
            PUBCHEM: '4399'
        })
        compound_name = 'gl:G10505'
        compound_info = _process_kegg_api_get_entity(compound_name, 'compound', self.hgnc_manager, self.chebi_manager)

        self.assertEqual(compound_info, {KEGG_ID: 'gl:G10505', KEGG_TYPE: 'compound'})

    def test_get_all_reactions(self):
        """Test reactions substrates, products."""
        glycolysis_genes, glycolysis_compounds, glycolysis_maps, glycolysis_orthologs = get_entity_nodes(
            self.glycolysis_tree,
            self.hgnc_manager,
            self.chebi_manager
        )
        substrates, products = get_all_reactions(self.glycolysis_tree, glycolysis_compounds)

        self.assertEqual(len(substrates), 35)
        self.assertEqual(len(products), 35)
        self.assertEqual(substrates['62'], ['90'])
        self.assertEqual(products['49'], ['102', '136'])
        self.assertEqual(substrates['48'], ['98', '136'])
        self.assertEqual(products['48'], ['99'])

    def test_get_reaction_edges(self):
        """Test reaction pathway edges on glycolysis."""
        glycolysis_genes, glycolysis_compounds, glycolysis_maps, glycolysis_orthologs = get_entity_nodes(
            self.glycolysis_tree,
            self.hgnc_manager,
            self.chebi_manager
        )
        substrate_dict, product_dict = get_all_reactions(self.glycolysis_tree, glycolysis_compounds)
        reactions = get_reaction_pathway_edges(self.glycolysis_tree, substrate_dict, product_dict)

        returned_reaction = None

        for k, v in reactions.items():
            for substrate, product, reaction in v:

                if substrate == ['85'] and product == ['92']:
                    returned_reaction = reaction

        self.assertEqual(len(reactions), 35)
        self.assertEqual(returned_reaction, 'reversible')
        self.assertEqual(reactions['48'], [(
            ['98', '136'],
            ['99'],
            'irreversible'
        )])

    def test_get_nodes(self):
        """Test nodes."""
        glycolysis_genes, glycolysis_compounds, glycolysis_maps, glycolysis_orthologs = get_entity_nodes(
            tree=self.glycolysis_tree,
            hgnc_manager=self.hgnc_manager,
            chebi_manager=self.chebi_manager
        )

        notch_genes, notch_compounds, notch_maps, notch_orthologs = get_entity_nodes(
            tree=self.notch_tree,
            hgnc_manager=self.hgnc_manager,
            chebi_manager=self.chebi_manager
        )
        ppar_genes, ppar_compounds, ppar_maps, ppar_orthologs = get_entity_nodes(
            self.ppar_tree,
            self.hgnc_manager,
            self.chebi_manager)

        glycolysis_nodes = xml_entities_to_bel(
            graph=self.glycolysis_empty_graph,
            genes_dict=glycolysis_genes,
            compounds_dict=glycolysis_compounds,
            maps_dict=glycolysis_maps,
            flattened=False
        )
        flat_glycolysis_nodes = xml_entities_to_bel(
            graph=self.glycolysis_empty_flatten_graph,
            genes_dict=glycolysis_genes,
            compounds_dict=glycolysis_compounds,
            maps_dict=glycolysis_maps,
            flattened=True
        )
        notch_nodes = xml_entities_to_bel(
            graph=self.notch_empty_graph,
            genes_dict=notch_genes,
            compounds_dict=notch_compounds,
            maps_dict=notch_maps,
            flattened=False
        )
        flat_notch_nodes = xml_entities_to_bel(
            graph=self.notch_empty_flatten_graph,
            genes_dict=notch_genes,
            compounds_dict=notch_compounds,
            maps_dict=notch_maps,
            flattened=True
        )
        ppar_nodes = xml_entities_to_bel(
            graph=self.ppar_empty_graph,
            genes_dict=ppar_genes,
            compounds_dict=ppar_compounds,
            maps_dict=ppar_maps,
            flattened=False
        )
        flat_ppar_nodes = xml_entities_to_bel(
            graph=self.ppar_empty_flatten_graph,
            genes_dict=ppar_genes,
            compounds_dict=ppar_compounds,
            maps_dict=ppar_maps,
            flattened=True
        )

        self.assertEqual(len(glycolysis_nodes), 73)
        self.assertEqual(len(flat_glycolysis_nodes), 73)
        self.assertEqual(len(notch_nodes), 24)
        self.assertEqual(len(flat_notch_nodes), 24)

        # Test un-flattened protein nodes
        self.assertEqual(glycolysis_nodes['53'], composite_abundance([
            protein(namespace=HGNC, name='PKLR', identifier='9020'),
            protein(namespace=HGNC, name='PKM', identifier='9021')
        ]))
        self.assertEqual(notch_nodes['22'], composite_abundance([
            protein(namespace=HGNC, name='NOTCH1', identifier='7881'),
            protein(namespace=HGNC, name='NOTCH2', identifier='7882'),
            protein(namespace=HGNC, name='NOTCH3', identifier='7883'),
            protein(namespace=HGNC, name='NOTCH4', identifier='7884')
        ]))

        # Test flattened list of protein nodes
        self.assertEqual(flat_glycolysis_nodes['53'], [
            protein(namespace=HGNC, name='PKLR', identifier='9020'),
            protein(namespace=HGNC, name='PKM', identifier='9021')
        ])
        self.assertEqual(flat_notch_nodes['22'], [
            protein(namespace=HGNC, name='NOTCH1', identifier='7881'),
            protein(namespace=HGNC, name='NOTCH2', identifier='7882'),
            protein(namespace=HGNC, name='NOTCH3', identifier='7883'),
            protein(namespace=HGNC, name='NOTCH4', identifier='7884')
        ])

        # Test pathway map nodes
        self.assertEqual(flat_glycolysis_nodes['54'],
                         bioprocess(namespace=KEGG, name='Citrate cycle (TCA cycle)', identifier='path:hsa00020')
                         )
        self.assertEqual(flat_notch_nodes['4'],
                         bioprocess(namespace=KEGG, name='MAPK signaling pathway', identifier='path:hsa04010')
                         )

        # Test un-flattened compound nodes
        self.assertEqual(ppar_nodes['48'], composite_abundance([
            abundance(namespace=CHEBI, name='9(S)-HODE', identifier='	34496'),
            abundance(namespace=CHEBI, name='13(S)-HODE', identifier='34154')
        ]))

        # Test flattened compound nodes
        self.assertEqual(flat_ppar_nodes['48'], [
            abundance(namespace=CHEBI, name='9(S)-HODE', identifier='	34496'),
            abundance(namespace=CHEBI, name='13(S)-HODE', identifier='34154')
        ])
        self.assertEqual(
            flat_glycolysis_nodes['85'],
            abundance(namespace=CHEBI, name='2-phospho-D-glyceric acid', identifier='17835')
        )

    def test_complex_node(self):
        """Test complex nodes on the notch pathway."""
        notch_genes, notch_compounds, notch_maps, notch_orthologs = get_entity_nodes(
            self.notch_tree,
            self.hgnc_manager,
            self.chebi_manager
        )
        complex_ids, flattened_complexes = get_complex_components(self.notch_tree, notch_genes, flattened=False)
        flat_complex_ids, flattened_complexes = get_complex_components(self.notch_tree, notch_genes, flattened=True)

        # not flatten part
        node_dict = xml_entities_to_bel(
            graph=self.notch_empty_graph,
            genes_dict=notch_genes,
            compounds_dict=notch_compounds,
            maps_dict=notch_maps,
            flattened=False
        )
        node_dict = xml_complexes_to_bel(
            graph=self.notch_empty_graph,
            node_dict=node_dict,
            complex_ids=complex_ids,
        )

        # Flatten part
        flat_node_dict = xml_entities_to_bel(
            graph=self.notch_empty_flatten_graph,
            genes_dict=notch_genes,
            compounds_dict=notch_compounds,
            maps_dict=notch_maps,
            flattened=True
        )
        flat_node_dict = xml_complexes_to_bel(
            graph=self.notch_empty_flatten_graph,
            node_dict=flat_node_dict,
            complex_ids=flat_complex_ids,
            flatten_complexes=flattened_complexes
        )

        self.assertEqual(len(node_dict), 28)
        self.assertEqual(len(flat_node_dict), 28)

    def test_bel_nodes(self):
        """Test transforming kgml into bel nodes."""
        notch_summary_flatten_nodes = count_functions(self.notch_bel_flatten)
        notch_summary_unflatten_edges = count_relations(self.notch_bel_unflatten)
        notch_summary_unflatten_nodes = count_functions(self.notch_bel_unflatten)

        glycolysis_summary_unflatten_nodes = count_functions(self.glycolysis_bel_unflatten)
        glycolysis_summary_flatten_nodes = count_functions(self.glycolysis_bel_flatten)

        ppar_bel_unflatten_nodes = count_functions(self.ppar_bel_unflatten)
        ppar_bel_unflatten_edges = count_relations(self.ppar_bel_unflatten)
        ppar_bel_flatten_nodes = count_functions(self.ppar_bel_flatten)
        ppar_bel_flatten_edges = count_relations(self.ppar_bel_flatten)

        self.assertEqual(notch_summary_unflatten_nodes['Protein'], 48)
        self.assertEqual(notch_summary_unflatten_nodes['Composite'], 15)
        self.assertEqual(notch_summary_unflatten_nodes['Complex'], 4)
        self.assertEqual(notch_summary_unflatten_nodes['BiologicalProcess'], 2)
        self.assertEqual(notch_summary_unflatten_edges['decreases'], 6)
        self.assertEqual(notch_summary_unflatten_edges['increases'], 9)
        self.assertEqual(notch_summary_unflatten_edges['association'], 1)

        self.assertEqual(notch_summary_flatten_nodes['Protein'], 48)
        self.assertEqual(notch_summary_flatten_nodes['Composite'], 0)
        self.assertEqual(notch_summary_flatten_nodes['Complex'], 4)
        self.assertEqual(notch_summary_flatten_nodes['BiologicalProcess'], 2)

        self.assertEqual(glycolysis_summary_flatten_nodes['Protein'], 68)
        self.assertEqual(glycolysis_summary_flatten_nodes['Composite'], 0)
        self.assertEqual(glycolysis_summary_flatten_nodes['Complex'], 0)
        self.assertEqual(glycolysis_summary_flatten_nodes['Abundance'], 31)
        self.assertEqual(glycolysis_summary_flatten_nodes['BiologicalProcess'], 7)

        self.assertEqual(glycolysis_summary_unflatten_nodes['Protein'], 68)
        self.assertEqual(glycolysis_summary_unflatten_nodes['Composite'], 18)
        self.assertEqual(glycolysis_summary_unflatten_nodes['Complex'], 0)
        self.assertEqual(glycolysis_summary_flatten_nodes['Abundance'], 31)
        self.assertEqual(glycolysis_summary_unflatten_nodes['BiologicalProcess'], 7)

        self.assertEqual(self.ppar_bel_unflatten.summary_dict()['Number of Nodes'], 5)
        self.assertEqual(ppar_bel_unflatten_nodes['Abundance'], 2)
        self.assertEqual(ppar_bel_unflatten_nodes['Composite'], 1)
        # 1 protein BEL node and 1 phosphorylated protein BEL node
        self.assertEqual(ppar_bel_unflatten_nodes['Protein'], 2)
        self.assertEqual(self.ppar_bel_unflatten.summary_dict()['Number of Edges'], 5)
        self.assertEqual(ppar_bel_unflatten_edges['increases'], 2)
        self.assertEqual(ppar_bel_unflatten_edges['hasComponent'], 2)
        self.assertEqual(ppar_bel_unflatten_edges['hasVariant'], 1)

        self.assertEqual(self.ppar_bel_flatten.summary_dict()['Number of Nodes'], 4)
        self.assertEqual(ppar_bel_flatten_nodes['Abundance'], 2)
        self.assertEqual(ppar_bel_flatten_nodes['Composite'], 0)
        self.assertEqual(ppar_bel_flatten_nodes['Protein'], 2)
        self.assertEqual(self.ppar_bel_flatten.summary_dict()['Number of Edges'], 5)
        self.assertEqual(ppar_bel_flatten_edges['increases'], 4)
        self.assertEqual(ppar_bel_flatten_edges['hasVariant'], 1)
