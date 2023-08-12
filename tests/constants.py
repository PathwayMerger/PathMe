# -*- coding: utf-8 -*-

"""Tests constants."""

import logging
import os
import tempfile

from bio2bel.testing import TemporaryConnectionMixin
from bio2bel_chebi import Manager as ChebiManager
from bio2bel_hgnc import Manager as HgncManager
from bio2bel_kegg.manager import Manager
from pathme.kegg.convert_to_bel import kgml_to_bel
from pathme.kegg.kegg_xml_parser import import_xml_etree
from pybel import BELGraph

logger = logging.getLogger(__name__)

TEST_FOLDER = os.path.dirname(os.path.realpath(__file__))
KEGG_TEST_RESOURCES = os.path.join(TEST_FOLDER, 'resources', 'kegg')
WP_TEST_RESOURCES = os.path.join(TEST_FOLDER, 'resources', 'wp')
REACTOME_TEST_RESOURCES = os.path.join(TEST_FOLDER, 'resources', 'reactome')

GLYCOLYSIS_XML = os.path.join(KEGG_TEST_RESOURCES, 'hsa00010.xml')
NOTCH_XML = os.path.join(KEGG_TEST_RESOURCES, 'hsa04330.xml')
PPAR_XML = os.path.join(KEGG_TEST_RESOURCES, '03320_cpd_test.xml')

WP22 = os.path.join(WP_TEST_RESOURCES, 'WP22.ttl')
WP706 = os.path.join(WP_TEST_RESOURCES, 'WP706.ttl')
WP1871 = os.path.join(WP_TEST_RESOURCES, 'WP1871.ttl')
WP2799 = os.path.join(WP_TEST_RESOURCES, 'WP2799.ttl')
WP2359 = os.path.join(WP_TEST_RESOURCES, 'WP2359_mod.ttl')

dir_path = os.path.dirname(os.path.realpath(__file__))
resources_path = os.path.join(dir_path, 'resources')

pathways = os.path.join(resources_path, 'hsa.txt')
protein_pathway_url = os.path.join(resources_path, 'pathway_gene.txt')

hgnc_test_path = os.path.join(resources_path, 'hgnc_test.json')
chebi_test_path = os.path.join(resources_path, 'chebi_test.tsv.gz')


class KeggTest(TemporaryConnectionMixin):
    """A test case with a populated HGNC/CheBI databases for KEGG parser."""

    @classmethod
    def setUpClass(cls):
        """Create temporary file."""
        logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")
        logger.setLevel(logging.INFO)

        """Create temporary file"""

        cls.fd, cls.path = tempfile.mkstemp()
        cls.connection = 'sqlite:///' + cls.path

        # create temporary database
        cls.manager = Manager(cls.connection)

        """HGNC Manager"""

        cls.hgnc_manager = HgncManager(engine=cls.manager.engine, session=cls.manager.session)
        cls.hgnc_manager.populate(hgnc_file_path=hgnc_test_path, use_hcop=False)

        logger.info('HGNC database loaded')

        """CHEBI Manager"""

        cls.chebi_manager = ChebiManager(engine=cls.manager.engine, session=cls.manager.session)
        cls.chebi_manager._populate_compounds(url=chebi_test_path)

        logger.info('ChEBI database loaded')

        cls.notch_tree = import_xml_etree(NOTCH_XML)
        cls.glycolysis_tree = import_xml_etree(GLYCOLYSIS_XML)
        cls.ppar_tree = import_xml_etree(PPAR_XML)

        logger.info('Loading notch unflatten')
        cls.notch_bel_unflatten = kgml_to_bel(NOTCH_XML, cls.hgnc_manager, cls.chebi_manager)

        logger.info('Loading notch flatten')
        cls.notch_bel_flatten = kgml_to_bel(NOTCH_XML, cls.hgnc_manager, cls.chebi_manager, flatten=True)

        logger.info('Loading glycolysis unflatten')
        cls.glycolysis_bel_unflatten = kgml_to_bel(GLYCOLYSIS_XML, cls.hgnc_manager, cls.chebi_manager)

        logger.info('Loading glycolysis flatten')
        cls.glycolysis_bel_flatten = kgml_to_bel(GLYCOLYSIS_XML, cls.hgnc_manager, cls.chebi_manager, flatten=True)

        logger.info('Loading PPAR unflatten')
        cls.ppar_bel_unflatten = kgml_to_bel(PPAR_XML, cls.hgnc_manager, cls.chebi_manager)

        logger.info('Loading PPAR flatten')
        cls.ppar_bel_flatten = kgml_to_bel(PPAR_XML, cls.hgnc_manager, cls.chebi_manager, flatten=True)

        cls.glycolysis_empty_graph = BELGraph(
            name='Glycolysis',
            version='1.0.0',
            description='Glycolysis',
            authors="Daniel Domingo-Fernández, Josep Marín-Llaó and Sarah Mubeen",
            contact='daniel.domingo.fernandez@scai.fraunhofer.de',
        )

        cls.glycolysis_empty_graph.graph['pathway_id'] = 'path:hsa00010'

        cls.glycolysis_empty_flatten_graph = BELGraph(
            name='Glycolysis flatten',
            version='1.0.0',
            description='Glycolysis',
            authors="Daniel Domingo-Fernández, Josep Marín-Llaó and Sarah Mubeen",
            contact='daniel.domingo.fernandez@scai.fraunhofer.de',
        )

        cls.glycolysis_empty_graph.graph['pathway_id'] = 'path:hsa00010'

        cls.notch_empty_graph = BELGraph(
            name='Notch',
            version='1.0.0',
            description='Notch signaling pathway',
            authors="Daniel Domingo-Fernández, Josep Marín-Llaó and Sarah Mubeen",
            contact='daniel.domingo.fernandez@scai.fraunhofer.de',
        )

        cls.notch_empty_graph.graph['pathway_id'] = 'path:hsa04330'

        cls.notch_empty_flatten_graph = BELGraph(
            name='Notch flatten',
            version='1.0.0',
            description='Notch signaling pathway',
            authors="Daniel Domingo-Fernández, Josep Marín-Llaó and Sarah Mubeen",
            contact='daniel.domingo.fernandez@scai.fraunhofer.de',
        )

        cls.notch_empty_flatten_graph.graph['pathway_id'] = 'path:hsa04330'

        cls.ppar_empty_graph = BELGraph(
            name='PPAR',
            version='1.0.0',
            description='PPAR signaling pathway',
            authors="Daniel Domingo-Fernández, Josep Marín-Llaó and Sarah Mubeen",
            contact='daniel.domingo.fernandez@scai.fraunhofer.de',
        )

        cls.ppar_empty_graph.graph['pathway_id'] = 'path:hsa03320'

        cls.ppar_empty_flatten_graph = BELGraph(
            name='PPAR flatten',
            version='1.0.0',
            description='PPAR signaling pathway',
            authors="Daniel Domingo-Fernández, Josep Marín-Llaó and Sarah Mubeen",
            contact='daniel.domingo.fernandez@scai.fraunhofer.de',
        )

        cls.ppar_empty_flatten_graph.graph['pathway_id'] = 'path:hsa03320'

    @classmethod
    def tearDownClass(cls):
        """Close the connection in the manager and deletes the temporary database."""
        cls.manager.drop_all()
        cls.hgnc_manager.drop_all()
        cls.manager.session.close()
        cls.hgnc_manager.session.close()
        super().tearDownClass()
