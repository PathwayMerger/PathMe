# -*- coding: utf-8 -*-

"""Export harmonized universe."""

import logging
import os
from typing import List

from pathme.utils import get_files_in_folder
from pybel import BELGraph
from pybel import from_pickle
from pybel import union
from pybel.dsl import Abundance, BiologicalProcess, CentralDogma, ListAbundance, Reaction
from tqdm import tqdm

logger = logging.getLogger(__name__)

WIKIPATHWAYS_BIOL_PROCESS = {
    "lipid biosynthesis", "hsc survival", "glycolysis & gluconeogenesis",
    "triacylglyceride  synthesis", "wnt canonical signaling", "regulation of actin skeleton",
    "fatty acid metabolism", "mrna processing major splicing pathway", "senescence",
    "monocyte differentiation", "pentose phosphate pathway", "ethanolamine  phosphate",
    "hsc differentiation", "actin, stress fibers and adhesion",
    "regulation of actin cytoskeleton", "s-phase progression", "g1-s transition",
    "toll-like receptor signaling pathway", "regulation of  actin cytoskeleton",
    "proteasome degradation", "apoptosis", "bmp pathway", "ampk activation",
    "g1/s checkpoint arrest", "mapk signaling pathway",
    "chromatin remodeling and  epigenetic modifications", "wnt signaling pathway",
    "ros production", "erbb signaling pathway", "shh pathway", "inflammation",
    "dna replication", "mrna translation", "oxidative stress",
    "cell cycle checkpoint activation", "gi/go pathway", "wnt pathway",
    "g1/s transition of mitotic cell cycle", "modulation of estrogen receptor signalling",
    "dna repair", "bmp canonical signaling", "igf and insuline signaling", "unfolded protein response", "cell death",
    "p38/mapk  pathway", "glycogen metabolism", "gnrh signal pathway",
    "the intra-s-phase checkpoint mediated arrest of cell cycle progression", "tca cycle",
    "mtor protein kinase signaling pathway", "proteasome  degradation pathway", "morphine metabolism", "hsc aging",
    "gastric pepsin release", "parietal cell production", "prostaglandin pathway", "cell cycle (g1/s)  progression",
    "notch pathway", "g2/m progression", "wnt signaling", "cell adhesion", "cell cycle progression", "egfr pathway",
    "cell cycle", "angiogenesis", "g2/m-phase checkpoint", "hsc self renewal", "26s proteasome  degradation",
    "mapk signaling", "immune system up or down regulation", "m-phase progression", "insulin signaling",
    "nf kappa b pathway", "cell cycle  progression", "gi pathway",
    "cd45+ hematopoietic-    derived cell    proliferation",
    "kreb's cycle", "glycogen synthesis", "apoptosis pathway",
    "g1/s progression", "inflammasome activation", "melanin biosynthesis", "proteasomal degradation",
    "g2/m checkpoint arrest",
    "g1/s cell cycle transition", "dna damage response", "gastric histamine release"
}

WIKIPATHWAYS_METAB = {
    "2,8-dihydroxyadenine", "8,11-dihydroxy-delta-9-thc", "adp-ribosyl", "cocaethylene", "dhcer1p",
    "ecgonidine", "f2-isoprostane", "fumonisins b1", "iodine", "l-glutamate", "lactosylceramide",
    "methylecgonidine", "n-acetyl-l-aspartate", "nad+", "nadph oxidase", "neuromelanin",
    "nicotinic acid (na)", "nmn", "pip2", "sphingomyelin", "thf"
}
WIKIPATHWAYS_NAME_NORMALIZATION = {
    "Ca 2+": "ca 2+", "acetyl coa": "acetyl-coa", "acetyl-coa(mit)": "acetyl-coa",
    "h20": "h2o"
}

# Entities in Reactome that required manual curation
BLACK_LIST_REACTOME = {"5'"}
REACTOME_PROT = {
    "phospho-g2/m transition proteins", "integrin alpha5beta1, integrin alphavbeta3, cd47",
    "food proteins", "activated fgfr2", "adherens junction-associated proteins",
    "pi3k mutants,activator:pi3k", "prolyl 3-hydroxylases", "gpi-anchored proteins", "c3d, c3dg, ic3b",
    "c4s/c6s chains", "activated fgfr1 mutants and fusions", "activated fgfr3 mutants", "protein",
    "cyclin a2:cdk2 phosphorylated g2/m transition protein", "c4c, c3f", "activated raf/ksr1",
    "activated fgfr1 mutants", "g2/m transition proteins", "lman family receptors", "cyclin",
    "usp12:wdr48:wdr20,usp26", "proteins with cleaved gpi-anchors", "activated fgfr2 mutants", "c4d, ic3b",
    "c5b:c6:c7, c8, c9", "cyclin a1:cdk2 phosphorylated g2/m transition protein",
    "genetically or chemically inactive braf", "il13-downregulated proteins", "activated fgfr4 mutants",
    "rna-binding protein in rnp (ribonucleoprotein) complexes", "effector proteins", "usp3, saga complex",
    'dephosphorylated "receiver" raf/ksr1'
}


def get_all_pickles(kegg_path, reactome_path, wikipathways_path):
    """Return a list with all pickle paths."""
    kegg_pickles = get_files_in_folder(kegg_path)

    if not kegg_pickles:
        logger.warning('No KEGG files found. Please create the BEL KEGG files')

    reactome_pickles = get_files_in_folder(reactome_path)

    if not reactome_pickles:
        logger.warning('No Reactome files found. Please create the BEL Reactome files')

    wp_pickles = get_files_in_folder(wikipathways_path)

    if not wp_pickles:
        logger.warning('No WikiPathways files found. Please create the BEL WikiPathways files')

    return kegg_pickles, reactome_pickles, wp_pickles


def get_universe_graph(kegg_path: str, reactome_path: str, wikipathways_path: str) -> BELGraph:
    """Return universe graph."""
    kegg_pickles, reactome_pickles, wp_pickles = get_all_pickles(kegg_path, reactome_path, wikipathways_path)

    all_pickles = kegg_pickles + reactome_pickles + wp_pickles

    logger.info(f'A total of {len(all_pickles)} will be merged into the universe')

    iterator = tqdm(all_pickles, desc='Creating universe')

    universe_list = []

    # Export KEGG
    for file in iterator:
        if not file.endswith('.pickle'):
            continue

        if file in kegg_pickles:
            graph = from_pickle(os.path.join(kegg_path, file))

        elif file in reactome_pickles:
            graph = from_pickle(os.path.join(reactome_path, file))

        elif file in wp_pickles:
            graph = from_pickle(os.path.join(wikipathways_path, file))

        else:
            logger.warning(f'Unknown pickle file: {file}')
            continue

        universe_list.append(graph)

    return union(universe_list)


def process_reactome_multiple_genes(genes: str) -> List:
    """Process a wrong ID with multiple identifiers"""
    gene_list = []
    for counter, gene in enumerate(genes):

        # Strip the ' gene' prefix
        gene = gene.strip().strip(' gene').strip(' genes')

        # First element is always OK
        if counter == 0:
            gene_list.append(gene)

        # If the identifier starts the same than the first one, it is right
        elif gene[:2] == genes[0][:2]:
            gene_list.append(gene)

        # If the identifier is longer than 2 it is a 'valid' HGNC symbol
        elif len(gene) > 2:
            gene_list.append(gene)

        # If they start different, it might have only a number (e.g., 'ABC1, 2, 3') so it needs to be appended
        elif gene.isdigit():
            gene_list.append(genes[0][:-1] + gene)

        # If the have only one letter (e.g., HTR1A,B,D,E,F,HTR5A)
        elif len(gene) == 1:
            gene_list.append(genes[0][:-1] + gene)

    return gene_list


def munge_reactome_gene(gene):
    """Process Reactome gene"""
    if "," in gene:
        return process_reactome_multiple_genes(gene.split(","))

    elif "/" in gene:
        return process_reactome_multiple_genes(gene.split("/"))

    return gene


def calculate_database_sets(nodes, database):
    """Calculate node sets for each modality in the database"""
    gene_nodes = set()
    mirna_nodes = set()
    metabolite_nodes = set()
    bp_nodes = set()

    for node in nodes:

        if isinstance(node, ListAbundance) or isinstance(node, Reaction) or not node.name:
            continue

        # Lower case name and strip quotes or white spaces
        name = node.name.lower().strip('"').strip()

        # Dealing with Genes/miRNAs
        if isinstance(node, CentralDogma):

            ##################
            # miRNA entities #
            ##################

            if name.startswith("mir"):

                # Reactome preprocessing to flat multiple identifiers
                if database == 'reactome':
                    reactome_cell = munge_reactome_gene(name)
                    if isinstance(reactome_cell, list):
                        for name in reactome_cell:
                            mirna_nodes.add(name.replace("mir-", "mir"))
                    else:
                        mirna_nodes.add(name.strip(' genes').replace("mir-", "mir"))

                    continue

                mirna_nodes.add(name.replace("mir-", "mir"))

            ##################
            # Genes entities #
            ##################

            else:
                # Reactome preprocessing to flat multiple identifiers
                if database == 'reactome':
                    reactome_cell = munge_reactome_gene(name)
                    if isinstance(reactome_cell, list):
                        for name in reactome_cell:
                            if name in BLACK_LIST_REACTOME:  # Filter entities in black list
                                continue
                            elif name.startswith("("):  # remove redundant parentheses
                                name = name.strip("(").strip(")")

                            gene_nodes.add(name)
                    else:
                        gene_nodes.add(name)
                    continue

                # WikiPathways and KEGG do not require any processing of genes
                if name in WIKIPATHWAYS_BIOL_PROCESS:
                    bp_nodes.add(name)
                    continue
                gene_nodes.add(name)

        #######################
        # Metabolite entities #
        #######################

        elif isinstance(node, Abundance):

            if database == 'wikipathways':
                # Biological processes that are captured as abundance in BEL since they were characterized wrong in WikiPathways
                if name in WIKIPATHWAYS_BIOL_PROCESS:
                    bp_nodes.add(name)
                    continue

                elif node.namespace in {'WIKIDATA', 'WIKIPATHWAYS', 'REACTOME'} and name not in WIKIPATHWAYS_METAB:
                    bp_nodes.add(name)
                    continue

                # Fix naming in duplicate entity
                if name in WIKIPATHWAYS_NAME_NORMALIZATION:
                    name = WIKIPATHWAYS_NAME_NORMALIZATION[name]

            elif database == 'reactome':
                # Curated proteins that were coded as metabolites
                if name in REACTOME_PROT:
                    gene_nodes.add(name)
                    continue

                # Flat multiple identifiers (this is not trivial because most of ChEBI names contain commas,
                # so a clever way to fix some of the entities is to check that all identifiers contain letters)
                elif "," in name and all(
                        string.isalpha()
                        for string in name.split(",")
                ):
                    for string in name.split(","):
                        metabolite_nodes.add(name)
                    continue

            metabolite_nodes.add(name)

        #################################
        # Biological Processes entities #
        #################################

        elif isinstance(node, BiologicalProcess):
            if name.startswith('title:'):
                name = name[6:]  # KEGG normalize

            bp_nodes.add(name)

    return gene_nodes, mirna_nodes, metabolite_nodes, bp_nodes
