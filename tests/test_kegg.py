# -*- coding: utf-8 -*-

"""Tests for converting KEGG."""

from pybel_tools.summary import edge_summary

from pathme.kegg.convert_to_bel import kegg_to_bel, xml_entities_to_bel, xml_complexes_to_bel
from pathme.kegg.kegg_xml_parser import *
from pybel import BELGraph
from pybel.struct.summary.node_summary import count_functions
from .constants import NOTCH_XML, GLYCOLYSIS_XML, PPAR_XML, DatabaseMixin


class TestKegg(DatabaseMixin):
    """Tests for dealing with the KEGG sub-module."""

    def setUp(self):
        """Parse two examples files."""
        self.notch_tree = import_xml_etree(NOTCH_XML)
        self.glycolysis_tree = import_xml_etree(GLYCOLYSIS_XML)

        log.info('Loading notch unflatten')
        self.notch_bel_unflatten = kegg_to_bel(NOTCH_XML, self.hgnc_manager, self.chebi_manager)

        log.info('Loading notch flatten')
        self.notch_bel_flatten = kegg_to_bel(NOTCH_XML, self.hgnc_manager, self.chebi_manager, flatten=True)

        log.info('Loading glycolysis unflatten')
        self.glycolysis_bel_unflatten = kegg_to_bel(GLYCOLYSIS_XML, self.hgnc_manager, self.chebi_manager)

        log.info('Loading glycolysis flatten')
        self.glycolysis_bel_flatten = kegg_to_bel(GLYCOLYSIS_XML, self.hgnc_manager, self.chebi_manager, flatten=True)

        log.info('Loading PPAR unflatten')
        self.ppar_bel_unflatten = kegg_to_bel(PPAR_XML, self.hgnc_manager, self.chebi_manager)

        log.info('Loading PPAR flatten')
        self.ppar_bel_flatten = kegg_to_bel(PPAR_XML, self.hgnc_manager, self.chebi_manager)

        self.glycolysis_empty_graph = BELGraph(
            name='Glycolysis',
            version='1.0.0',
            description='Glycolisis',
            pathway_id='Glycolisis',
            authors="Daniel Domingo-Fernández, Josep Marín-Llaó and Sarah Mubeen",
            contact='daniel.domingo.fernandez@scai.fraunhofer.de'
        )

        self.glycolysis_empty_flatten_graph = BELGraph(
            name='Glycolisis flatten',
            version='1.0.0',
            description='Glycolisis',
            pathway_id='Glycolisis',
            authors="Daniel Domingo-Fernández, Josep Marín-Llaó and Sarah Mubeen",
            contact='daniel.domingo.fernandez@scai.fraunhofer.de'
        )

        self.notch_empty_graph = BELGraph(
            name='Notch',
            version='1.0.0',
            description='Notch',
            pathway_id='Notch',
            authors="Daniel Domingo-Fernández, Josep Marín-Llaó and Sarah Mubeen",
            contact='daniel.domingo.fernandez@scai.fraunhofer.de'
        )

        self.notch_empty_flatten_graph = BELGraph(
            name='Notch flatten',
            version='1.0.0',
            description='Notch',
            pathway_id='Notch',
            authors="Daniel Domingo-Fernández, Josep Marín-Llaó and Sarah Mubeen",
            contact='daniel.domingo.fernandez@scai.fraunhofer.de'
        )

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
            {'kegg_id': 'hsa:5663',
             'kegg_type': 'gene',
             'HGNC symbol': 'PSEN1',
             'HGNC': '9508',
             'UniProt': 'P49768 A0A024R6A3'},
            {'kegg_id': 'hsa:5664',
             'kegg_type': 'gene',
             'HGNC symbol': 'PSEN2',
             'HGNC': '9509',
             'UniProt': 'P49810'}
        ])
        self.assertEqual(glycolysis_genes['45'], [
            {'kegg_id': 'hsa:10327',
             'kegg_type': 'gene',
             'HGNC symbol': 'AKR1A1',
             'HGNC': '380',
             'UniProt': 'P14550 V9HWI0'}
        ])
        self.assertEqual(glycolysis_compounds['83'], [
            {'compound_name': 'cpd:C00031', 'kegg_type': 'compound', 'ChEBI': '4167', 'ChEBI name': 'D-glucopyranose'}
        ])
        self.assertEqual(notch_maps['4'], [
            {'kegg_id': 'path:hsa04010', 'map_name': 'MAPK signaling pathway'}
        ])
        self.assertEqual(glycolysis_maps['52'], [
            {'kegg_id': 'path:hsa00620', 'map_name': 'Pyruvate metabolism'}
        ])
        self.assertEqual(notch_orthologs['7'], [
            {'kegg_id': 'ko:K04497', 'kegg_type': 'ortholog'}
        ])
        self.assertEqual(glycolysis_orthologs['73'], [
            {'kegg_id': 'ko:K01085', 'kegg_type': 'ortholog'},
            {'kegg_id': 'ko:K20866', 'kegg_type': 'ortholog'}
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
            {'kegg_id': 'hsa:3065',
             'kegg_type': 'gene',
             'HGNC symbol': 'HDAC1',
             'HGNC': '4852',
             'UniProt': 'Q13547 Q6IT96'},
            {'kegg_id': 'hsa:3066',
             'kegg_type': 'gene',
             'HGNC symbol': 'HDAC2',
             'HGNC': '4853',
             'UniProt': 'Q92769'},
            {'kegg_id': 'hsa:1487',
             'kegg_type': 'gene',
             'HGNC symbol': 'CTBP1',
             'HGNC': '2494',
             'UniProt': 'Q13363 X5D8Y5'},
            {'kegg_id': 'hsa:1488',
             'kegg_type': 'gene',
             'HGNC symbol': 'CTBP2',
             'HGNC': '2495',
             'UniProt': 'P56545'}
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
        compound_info = get_compound_info(compound_name, self.chebi_manager)

        self.assertEqual(compound_info['ChEBI'], '17719')
        self.assertEqual(compound_info['ChEBI name'], 'beta-D-glucose 6-phosphate')

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

    def test_get_reaction_edges(self):
        """Test reaction pathway edges on glycolysis."""
        glycolysis_genes, glycolysis_compounds, glycolysis_maps, glycolysis_orthologs = get_entity_nodes(
            self.glycolysis_tree,
            self.hgnc_manager,
            self.chebi_manager
        )
        substrate_dict, product_dict = get_all_reactions(self.glycolysis_tree, glycolysis_compounds)
        reactions = get_reaction_pathway_edges(self.glycolysis_tree, substrate_dict, product_dict)

        returned_reaction, returned_product = None, None

        for k, v in reactions.items():
            for substrate, product, reaction in v:

                if substrate == ['85'] and product == ['92']:
                    returned_reaction = reaction
                if substrate == ['136', '98']:
                    returned_product = product

        self.assertEqual(len(reactions), 35)
        self.assertEqual(returned_reaction, 'reversible')
        self.assertEqual(returned_product, ['99'])

    def test_get_pathway_edges(self):
        """Test edges."""
        notch_genes, notch_compounds, notch_maps, notch_orthologs = get_entity_nodes(
            self.notch_tree,
            self.hgnc_manager,
            self.chebi_manager
        )
        relations = get_all_relationships(self.notch_tree)
        complexes = get_entities_in_complex(self.notch_tree, notch_genes)
        edges = get_pathway_edges(notch_genes, relations, set_of_complexes=complexes)

        returned_relation, complex_relation = None, None

        for source, target, relation in edges:
            if source == 'hsa:11317' and target == 'hsa:171558':
                returned_relation = relation

            if source == ('hsa:3065', 'hsa:1487') and target == 'hsa:11317':
                complex_relation = relation

        self.assertEqual(returned_relation, 'expression')
        self.assertEqual(complex_relation, 'inhibition')

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

        path_name = None

        self.assertEqual(len(glycolysis_nodes), 73)
        self.assertEqual(len(flat_glycolysis_nodes), 73)
        self.assertEqual(len(notch_nodes), 24)

        # Test flattened list of protein nodes
        self.assertEqual(flat_glycolysis_nodes['53'], [
            {'function': 'Protein', 'namespace': 'HGNC', 'name': 'PKLR', 'identifier': '9020'},
            {'function': 'Protein', 'namespace': 'HGNC', 'name': 'PKM', 'identifier': '9021'}
        ])

        # Test pathway map nodes
        pathway_info = glycolysis_nodes['54']
        if pathway_info['identifier'] == 'path:hsa00020':
            path_name = pathway_info['name']

        self.assertEqual(path_name, 'Citrate cycle (TCA cycle)')
        self.assertEqual(
            glycolysis_nodes['85'],
            {'function': 'Abundance', 'namespace': 'ChEBI', 'identifier': '17835'}
        )

        # Test compound nodes
        self.assertEqual(glycolysis_nodes['85'], {
            'function': 'Abundance', 'namespace': 'ChEBI', 'identifier': '17835'
        })

    def test_complex_node(self):
        """Test complex nodes on the notch pathway."""
        notch_genes, notch_compounds, notch_maps, notch_orthologs = get_entity_nodes(
            self.notch_tree,
            self.hgnc_manager,
            self.chebi_manager
        )
        complex_ids, flattened_complexes = get_complex_components(self.notch_tree, notch_genes, flattened=False)
        flat_complex_ids, flattened_complexes = get_complex_components(self.notch_tree, notch_genes, flattened=True)

        node_dict = xml_entities_to_bel(
            graph=self.notch_empty_graph,
            genes_dict=notch_genes,
            compounds_dict=notch_compounds,
            maps_dict=notch_maps,
            flattened=False
        )
        node_dict = xml_complexes_to_bel(node_dict, complex_ids, flattened_complexes)

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
        notch_summary_unflatten_edges = edge_summary.count_relations(self.notch_bel_unflatten)
        notch_summary_unflatten_nodes = count_functions(self.notch_bel_unflatten)

        glycolysis_summary_unflatten_nodes = count_functions(self.glycolysis_bel_unflatten)
        glycolysis_summary_flatten_nodes = count_functions(self.glycolysis_bel_flatten)

        ppar_bel_unflatten_nodes = count_functions(self.ppar_bel_unflatten)
        ppar_bel_unflatten_edges = edge_summary.count_relations(self.ppar_bel_unflatten)
        ppar_bel_flatten_nodes = count_functions(self.ppar_bel_flatten)
        ppar_bel_flatten_edges = edge_summary.count_relations(self.ppar_bel_flatten)

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

        self.assertEqual(self.ppar_bel_flatten.summary_dict()['Number of Nodes'], 4)
        self.assertEqual(ppar_bel_unflatten_nodes['Abundance'], 2)
        self.assertEqual(ppar_bel_unflatten_nodes['Composite'], 1)
        self.assertEqual(ppar_bel_unflatten_nodes['Protein'], 1)
        self.assertEqual(self.ppar_bel_flatten.summary_dict()['Number of Edges'], 3)
        self.assertEqual(ppar_bel_unflatten_edges['increases'], 1)
        self.assertEqual(ppar_bel_unflatten_edges['hasComponent'], 2)

        self.assertEqual(self.ppar_bel_flatten.summary_dict()['Number of Nodes'], 3)
        self.assertEqual(ppar_bel_flatten_nodes['Abundance'], 2)
        self.assertEqual(ppar_bel_flatten_nodes['Composite'], 0)
        self.assertEqual(ppar_bel_flatten_nodes['Protein'], 1)
        self.assertEqual(self.ppar_bel_flatten.summary_dict()['Number of Edges'], 2)
        self.assertEqual(ppar_bel_flatten_edges['increases'], 2)
