# -*- coding: utf-8 -*-

"""Methods to normalize names across databases."""

import logging
from collections import defaultdict
from typing import List

from networkx import relabel_nodes
from pathme.constants import REACTOME, WIKIPATHWAYS
from pathme.pybel_utils import multi_relabel
from pybel import BELGraph
from pybel.dsl import Abundance, BiologicalProcess, CentralDogma, ListAbundance, Reaction, MicroRna, Protein

logger = logging.getLogger(__name__)

# Curated list of BPs from WikiPathways
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

# WikiPathways Abundances that were categorized with the wrong function
WIKIPATHWAYS_METAB = {
    "2,8-dihydroxyadenine", "8,11-dihydroxy-delta-9-thc", "adp-ribosyl", "cocaethylene", "dhcer1p",
    "ecgonidine", "f2-isoprostane", "fumonisins b1", "iodine", "l-glutamate", "lactosylceramide",
    "methylecgonidine", "n-acetyl-l-aspartate", "nad+", "nadph oxidase", "neuromelanin",
    "nicotinic acid (na)", "nmn", "pip2", "sphingomyelin", "thf"
}

# Normalize name WikiPathways
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


def process_reactome_multiple_genes(genes: str) -> List:
    """Process a wrong ID with multiple identifiers."""
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
    """Process Reactome gene."""
    if "," in gene:
        return process_reactome_multiple_genes(gene.split(","))

    elif "/" in gene:
        return process_reactome_multiple_genes(gene.split("/"))

    return gene


def normalize_graph_names(graph: BELGraph, database: str) -> None:
    """Normalize graph names."""

    # Victim to Survivor (one to one node) mapping
    one_to_one_mapping = {}
    # Victim to Survivors (one to many nodes) mapping
    one_to_many_mapping = defaultdict(set)

    for node in graph.nodes():

        # Skip ListAbundances and Reactions since they do not have a name
        if isinstance(node, ListAbundance) or isinstance(node, Reaction) or not node.name:
            continue

        # Normalize names: Lower case name and strip quotes or white spaces
        lower_name = node.name.lower().strip('"').strip()

        # Dealing with Genes/miRNAs
        if isinstance(node, CentralDogma):

            ##################
            # miRNA entities #
            ##################

            if lower_name.startswith("mir"):

                # Reactome preprocessing to flat multiple identifiers
                if database == REACTOME:
                    reactome_cell = munge_reactome_gene(lower_name)
                    if isinstance(reactome_cell, list):
                        for lower_name in reactome_cell:
                            one_to_many_mapping[node].add(
                                MicroRna(
                                    node.namespace, name=lower_name.replace("mir-", "mir"), identifier=node.identifier
                                ))

                    one_to_one_mapping[node] = MicroRna(
                        node.namespace,
                        name=lower_name.strip(' genes').replace("mir-", "mir")  # Special case for Reactome
                    )
                    continue

                # KEGG and Reactome
                one_to_one_mapping[node] = MicroRna(
                    node.namespace, name=node.name.replace("mir-", "mir"), identifier=node.identifier
                )

            ##################
            # Genes entities #
            ##################

            else:
                # Reactome preprocessing to flat multiple identifiers
                if database == REACTOME:
                    reactome_cell = munge_reactome_gene(lower_name)
                    if isinstance(reactome_cell, list):
                        for lower_name in reactome_cell:
                            if lower_name in BLACK_LIST_REACTOME:  # Filter entities in black list
                                continue
                            elif lower_name.startswith("("):  # remove redundant parentheses
                                lower_name = lower_name.strip("(").strip(")")

                            one_to_many_mapping[node].add(
                                Protein(node.namespace, name=lower_name, identifier=node.identifier)
                            )
                    else:
                        one_to_one_mapping[node] = Protein(node.namespace, name=lower_name, identifier=node.identifier)

                    continue

                # WikiPathways and KEGG do not require any processing of genes
                elif database == WIKIPATHWAYS and lower_name in WIKIPATHWAYS_BIOL_PROCESS:
                    one_to_one_mapping[node] = BiologicalProcess(
                        node.namespace, name=lower_name, identifier=node.identifier
                    )
                    continue

                one_to_one_mapping[node] = Protein(node.namespace, name=lower_name, identifier=node.identifier)

        #######################
        # Metabolite entities #
        #######################

        elif isinstance(node, Abundance):

            if database == 'wikipathways':
                # Biological processes that are captured as abundance in BEL since they were characterized wrong in WikiPathways
                if lower_name in WIKIPATHWAYS_BIOL_PROCESS:
                    one_to_one_mapping[node] = BiologicalProcess(
                        node.namespace, name=lower_name, identifier=node.identifier
                    )
                    continue

                # Abundances to BiologicalProcesses
                elif node.namespace in {'WIKIDATA', 'WIKIPATHWAYS', 'REACTOME'} \
                        and lower_name not in WIKIPATHWAYS_METAB:
                    one_to_one_mapping[node] = BiologicalProcess(
                        node.namespace, name=lower_name, identifier=node.identifier
                    )
                    continue

                # Fix naming in duplicate entity
                if lower_name in WIKIPATHWAYS_NAME_NORMALIZATION:
                    lower_name = WIKIPATHWAYS_NAME_NORMALIZATION[lower_name]

            elif database == REACTOME:
                # Curated proteins that were coded as metabolites
                if lower_name in REACTOME_PROT:
                    one_to_one_mapping[node] = Protein(
                        node.namespace, name=lower_name, identifier=node.identifier
                    )
                    continue

                # Flat multiple identifiers (this is not trivial because most of ChEBI names contain commas,
                # so a clever way to fix some of the entities is to check that all identifiers contain letters)
                elif "," in lower_name and all(
                        string.isalpha()
                        for string in lower_name.split(",")
                ):
                    for string in lower_name.split(","):
                        one_to_many_mapping[node].add(
                            Abundance(node.namespace, name=string, identifier=node.identifier)
                        )
                    continue

            one_to_one_mapping[node] = Abundance(node.namespace, name=lower_name, identifier=node.identifier)

        #################################
        # Biological Processes entities #
        #################################

        elif isinstance(node, BiologicalProcess):
            # KEGG normalize name by removing the title prefix
            if lower_name.startswith('title:'):
                lower_name = lower_name[6:]

            one_to_one_mapping[node] = BiologicalProcess(
                node.namespace, name=lower_name, identifier=node.identifier
            )

    relabel_nodes(graph, one_to_one_mapping)
    multi_relabel(graph, one_to_many_mapping)
