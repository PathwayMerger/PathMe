# -*- coding: utf-8 -*-

"""This module contains the methods that handle the mappings stored in ComPath resources"""

import pandas as pd


def get_pathways_from_statement(mapping_statement, mapping_type):
    """Return the subject, object of the mapping.

    :param str mapping_statement: statement
    :param str mapping_type: type of relationship
    :rtype: tuple[str,str]
    """
    _pathways = mapping_statement.split(mapping_type)

    return _pathways[0].strip(), _pathways[1].strip()


def remove_star_from_pathway_name(pathway_name):
    """Remove the star that label the reference pathway in isPartOf statements.

    :param str pathway_name: pathway name
    :rtype: str
    """
    return pathway_name.replace("*", "").strip()


def parse_part_of_mapping(mapping_statement):
    """Return the pathways of a hierarchical mapping.

    :param str mapping_statement: statement
    :rtype: tuple[str,str]
    """
    pathway_1, pathway_2 = get_pathways_from_statement(mapping_statement, 'isPartOf')

    if "*" in pathway_1:
        pathway_1 = remove_star_from_pathway_name(pathway_1)
        return pathway_1, pathway_2

    pathway_2 = remove_star_from_pathway_name(pathway_2)
    return pathway_2, pathway_1


def get_mapped_pathways(dataframe):
    """Get pathways with mappings.

    :param pandas.DataFrame dataframe: data frame mappings
    :returns: pathway mappings
    :rtype: list[tuple[str][str][str]]
    """
    mappings = list()

    for _, row in dataframe.iterrows():

        equivalent_to_mappings = row['equivalentTo Mappings']

        if not pd.isnull(equivalent_to_mappings):

            for mapping_statement in equivalent_to_mappings.split("\n"):

                if mapping_statement == '':
                    continue

                reference_pathway, compared_pathway = get_pathways_from_statement(mapping_statement, "equivalentTo")

                mappings.append((reference_pathway, "equivalentTo", compared_pathway))

        is_part_of_mappings = row['isPartOf Mappings']

        if not pd.isnull(is_part_of_mappings):

            for mapping_statement in is_part_of_mappings.split('\n'):

                if mapping_statement == '':
                    continue

                reference_pathway, compared_pathway = parse_part_of_mapping(mapping_statement)

                mappings.append((reference_pathway, "isPartOf", compared_pathway))

    return mappings
