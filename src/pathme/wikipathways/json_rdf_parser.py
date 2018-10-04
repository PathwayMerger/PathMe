# -*- coding: utf-8 -*-

"""This module contains the custom parser for RDF."""
import json

import rdflib

from pathme.utils import parse_rdf
from pathme.wikipathways.utils import convert_to_nx

"""RDF CUSTOM PARSER FUNCTIONS"""

"""URI Parsers"""


def parse_id_uri(uri):
    """Get the components of a given uri (with identifier at the last position).

    :param str uri: URI
    :returns: prefix (ex: http://rdf.wikipathways.org/...)
    :returns: prefix_namespaces: if there are many namespaces, until the penultimate (ex: .../Pathway/WP22_r97775/...)
    :returns: namespace: if there are many namespaces, the last (ex: .../Interaction/)
    :returns: identifier (ex: .../c562c/)
    :rtype: tuple[str,str,str,str]
    """
    # Split the uri str by '/'.
    splitted_uri = uri.split('/')

    # Get the uri components into different variables.
    prefix = '/'.join(splitted_uri[0:3])
    prefix_namespaces = '/'.join(splitted_uri[3:-2])
    namespace = splitted_uri[-2]
    identifier = splitted_uri[-1]

    return prefix, prefix_namespaces, namespace, identifier


def parse_namespace_uri(uri):
    """Get the prefix and namespace of a given URI (without identifier, only with a namspace at last position).

    :param str uri: URI
    :returns: prefix (ex: http://purl.org/dc/terms/...)
    :returns: namespace (ex: .../isPartOf)
    :rtype: tuple[str,str]
    """
    # Split the uri str by '/'.
    splited_uri = uri.split('/')

    # Get the uri prefix and namespace into different variables.
    prefix = ('/').join(splited_uri[0:-1])
    namespace = splited_uri[-1]
    vocabulary = namespace.split('#')[-1]

    return prefix, namespace, vocabulary


"""Match Cases"""


def match_entry_type(types):
    """Get the type identifier from the type attributes of an entry.

    :param set types: set with wp vocabularies as string.
    :returns: entry_type: type identifier, also used as first key to the main graph (ex: 'nodes' -> graph[entry_type][node_id])
    :rtype: Optional[str]
    """
    # Matches the vocabularies set from wp (could be indicated one or more types) specified in RDF as URI namespaces (previous URI parsing).
    if types == {'core#Collection', 'wp#Pathway'} or types == {'wp#PublicationReference'}:
        return 'pathway_info'

    elif types == {'wp#DataNode', 'wp#Protein'} \
            or types == {'wp#DataNode', 'wp#Metabolite'} \
            or types == {'wp#DataNode', 'wp#GeneProduct'} \
            or types == {'wp#DataNode'}:
        return 'nodes'

    elif types == {'wp#Interaction', 'wp#DirectedInteraction'} \
            or types == {'wp#Catalysis', 'wp#DirectedInteraction', 'wp#Interaction'} \
            or types == {'wp#Interaction', 'wp#ComplexBinding', 'wp#Binding'}:
        return 'interactions'

    elif types == {'wp#Interaction'} or types == {'wp#Complex'}:
        return 'complex'

    raise Exception('Entry type/s not recognised %s', types)


def match_attribute_label(attribute_namespace):
    """Assign differently an attribute label depending on the namespace of the attribute label (attribute key).

    :param str attribute_namespace: the namespace of the attribute specification
    :returns: attribute_label: label to identify the attribute
    :rtype: Optional[str]
    """
    # Depending on the attribute_namespace, of the following cases:

    # Return the attribute_namespace (is the variable that is being checked, so will be one of the list).
    if attribute_namespace in {'title', 'source', 'identifier', 'description', 'isPartOf', 'hasVersion', 'page',
                               'references'}:
        return attribute_namespace

    # Return the attribute_namespace without the voacabulary resource specification (ex: wp#, split #).
    elif attribute_namespace in {'wp#isAbout', 'wp#participants', 'wp#pathwayOntologyTag', 'wp#ontologyTag',
                                 'wp#organism', 'wp#organismName', 'rdf-schema#label'}:
        return attribute_namespace.split('#')[1]

    # Return the value_namespace. The value of the atribute is also specified as a uri, with the value as the identifier. So in the following cases, the value namespace of the value uri will be later assigned as the attribute label.
    elif attribute_namespace in {'wp#bdbEnsembl', 'wp#bdbEntrezGene', 'wp#bdbHgncSymbol', 'wp#bdbPubChem', 'wp#bdbHmdb',
                                 'wp#bdbUniprot', 'wp#bdbChEBI', 'wp#bdbChemspider', 'wp#bdbWikidata'}:
        return 'value_namespace'

    else:
        raise Exception('Attribute namespace not recognised %s', attribute_namespace)


def match_entry(entry):
    """For a given entry, get the identifier of the entry and also the entry_type.

    :param dict entry:
    :returns: entry_id: for nodes would be the ncbi id and edges the interaction wp id
    :returns: entry_type:  type identifier, also used as first key to the main graph (ex: 'nodes' -> graph[entry_type][node_id])
    :rtype: Optional[tuple[str,str]]
    """
    # Get the components of the entry identifier URI. The URI identifier will be used as the entry id and the namespace for type identification, wheras the prefix to match the entry handling
    uri = entry['@id']
    prefix, _, namespace, entry_id = parse_id_uri(uri)

    # Check if the prefix is recognized (could also be treated different for specific prefix cases)
    if prefix not in {'http://identifiers.org', 'http://rdf.wikipathways.org'}:
        raise Exception('Entry not recognised: %s', prefix)

    # Get the entry type id calling the get_entry_type function, if the the entry has an attribute with type
    if '@type' in entry:
        entry_type = get_entry_type(entry['@type'])

    # Assign literally the entry type when the entry has no type attribute,
    # if the entry id namespace is recognized (like the pathway id entry), else raise an exeption
    else:
        if namespace == 'wikipathways':
            entry_type = 'pathway_info'
        else:
            raise Exception('Entry type not recognised: %s %s', prefix, namespace)

    return entry_id, entry_type


def match_attribute(uri):
    """For a given attribute @id URI, get the label to be assigned to the attribute.

    :param str uri: attribute_label as URI
    :returns: attribute_label: label to be assigned to the attribute
    :rtype: Optional[str]
    """
    # Get the prefix and namespace of the @id uri
    attribute_prefix, attribute_namespace, _ = parse_namespace_uri(uri)

    # TODO: /dc/terms attribute_prefix_namespaces

    # Check if the prefix is recognized (could also be treated different for specific prefix cases)
    if attribute_prefix not in {
        'http://www.w3.org/2000/01', 'http://identifiers.org',
        'http://vocabularies.wikipathways.org',
        'http://purl.org/dc/terms',
        'http://purl.org/dc/elements/1.1',
        'http://xmlns.com',
        'http://xmlns.com/foaf/0.1',
        'http://purl.org/pav'
    }:
        raise Exception('Invalid attribute prefix %s', attribute_prefix, attribute_namespace)

    # Get the attribute_type depending on the attribute_namespace, calling the match_attribute_label function.
    return match_attribute_label(attribute_namespace)


"""Getters and Setters"""


def get_entry_attribute_value(entry_label, node_id, attribute_label, graph):
    """For a given node (if entry type is nodes), entry_type and attribute, get the associated value in the graph object.

    :param str node_id: ncbi id
    :param str attribute_label: ex: 'source', 'title', 'description'...
    :param str entry_label: ex: 'nodes', 'title', 'description'...
    :param dict graph:
    :returns: attribute_value: associated value in the graph object
    :rtype: Optional[str]
    """
    # Check for each queried parameter if it is in the graph as key to get the value.
    if entry_label in graph:
        if entry_label == 'nodes' and node_id in graph[entry_label] and attribute_label in graph[entry_label][node_id]:
            return graph[entry_label][node_id][attribute_label]

        elif attribute_label in graph[entry_label]:
            return graph[entry_label][attribute_label]

    raise Exception('Error in get node attribute: %s value: %s', entry_label, attribute_label, node_id)


def set_entry_attribute(entry_type, node_id, attribute_label, value, graph):
    """For a given node (if the entry type is 'nodes'), entry_type and attribute, set the associated value in the graph object.

    :param str node_id: ncbi id
    :param str attribute_label: attribute label
    :param str entry_label: entry label
    :param str value: value
    :param dict graph: graph object
    """
    # Check the entry_type to handle the graph value insertion differently (A node has the node_id as an extra key).
    if entry_type == 'nodes':

        if node_id not in graph['nodes']:
            graph['nodes'][node_id] = {}

        graph['nodes'][node_id][attribute_label] = value

    elif entry_type == 'pathway_info':
        graph['pathway_info'][attribute_label] = value

    else:
        raise Exception('Error in set node attribute: %s value: %s %s', entry_type, node_id, attribute_label)


def set_interaction(entry, graph):
    """For a given interaction entry sets the source and target (and the interaction id) into the pathway graph.

    :param dict entry: entry whose type has been previous identified as 'interactions'
    :param dict graph: pathway network graph object
    """
    # Get the participants of the interaction, picking directly the node identifiers from the entry attributes with the corresponding argument URI namespace (wp#source or wp#interaction). After, for each participant, parse the URI value to obtine the nodes identifier.
    # Get the node source id
    uri_source_id = entry['http://vocabularies.wikipathways.org/wp#source'][0]['@id']
    _, _, _, source_id = parse_id_uri(uri_source_id)

    # Get the node target id
    uri_target_id = entry['http://vocabularies.wikipathways.org/wp#target'][0]['@id']
    _, _, _, target_id = parse_id_uri(uri_target_id)

    # Also get the isAbout as the identifier of the interaction. For now will be the literal URI, due to no further information (like inhibits, increments) is indicated
    uri_interaction_type = entry['http://vocabularies.wikipathways.org/wp#isAbout'][0]['@id']

    # Finally, add directy to the interactions set of the graph the three identifiers of the interaction as a tupple.
    graph['interactions'].add((source_id, target_id, uri_interaction_type))


def get_entry_type(types):
    """For a set of uris that indicate the entry's type, get the type identifier (call match_entry_type).

    :param set types: set of uris (from the entry's @type attribute) that indicate the entry's type
    :returns: entry_type: type identifier, also used as first key to the main graph (ex: 'nodes' -> graph[entry_type][node_id])
    :rtype: str
    """
    types_set = set()

    # For each type indicated in the attribute of the entry, add the namespace to a set
    for typ in types:
        _, namespace, _ = parse_namespace_uri(typ)

        types_set.add(namespace)

    # Get the type identifier in function of the namespaces of the type set (call match_entry_type).
    entry_type = match_entry_type(types_set)

    return entry_type


"""Statements Parser"""


def parse_attribute_values(entry_label, entry_id, attribute_values, attribute_label, graph):
    """For each value in attribute_values, taking into account the attribute_label type (if it is specified value_namespace thus would be the value namespace), adds a new entry to the graph (calling set_entry_attribute methode) beeing the last level of parsing. The value is added as a set if there are multiple values for the same attribute_labe or as a sigle value.

    :param str entry_label: entry label
    :param str entry_id: entry identifier
    :param list[dict] attribute_values: values
    :param str attribute_label: label
    :param dict graph: graph object
    """
    attribute_value = set()

    for attribute_raw_value in attribute_values:
        for value_label, value in attribute_raw_value.items():
            if value_label == '@id':
                _, _, value_namespace, value_identifier = parse_id_uri(value)
                attribute_value.add(value_identifier)

                if attribute_label == 'value_namespace':
                    attribute_label = value_namespace

            elif value_label == '@value':
                attribute_value.add(value)

            elif value_label == '@language':
                set_entry_attribute(entry_label, entry_id, 'language', value, graph)
            else:
                raise Exception('Error with attribute {}'.format(value_label))

    if len(attribute_value) == 1:
        attribute_value = list(attribute_value)[0]

    set_entry_attribute(entry_label, entry_id, attribute_label, attribute_value, graph)


def parse_attributes(entry, entry_type, entry_id, graph):
    """For each attribute in attributes, if is labbeled as a uri (not in {'@id', '@value', '@type'}) gets the attribute_type (calling the correspondent function) and calls the next statement (parse_attribute_values).

    :param dict attributes: attributes of the entity
    :param str entry_type: entity type
    :param str entry_id: entry identifier
    :param dict graph: graph object
    """
    for attribute_label, values in entry.items():
        if attribute_label not in {'@id', '@value', '@type'}:
            attribute_type = match_attribute(attribute_label)
            parse_attribute_values(entry_type, entry_id, values, attribute_type, graph)


def generate_empty_pathway_graph():
    """Test set entry attribute."""
    graph = {}
    graph['interactions'] = set()
    graph['nodes'] = {}
    graph['pathway_info'] = {}

    return graph


def parse_entries(entries):
    """Create the graph object which will be finally returned full, with the values of the statements parser calls (entry -> attributes -> values). First the type of the entry and the id is obtined with match_entry, and according the type retrived is haddled the entry in a particular manner. If the entry is 'interactions' type, only will be called the methode set_interaction, but if not in depth parser if the different levels will be called (parse_attributes and parse_attribute_values).

    :param list[dict] entries:
    :rtype: networkx.MultiDiGraph
    """
    # Create the graph object
    graph = generate_empty_pathway_graph()

    for entry in entries:
        entry_id, entry_type = match_entry(entry)

        if entry_type == 'interactions':
            set_interaction(entry, graph)
        elif entry_type != 'complex':
            parse_attributes(entry, entry_type, entry_id, graph)

    return graph

def convert_json(graph: rdflib.Graph):
    """Convert from rdflib importated graph object to python data structure (list of dicts of each entry).

    :param rdflib.graph graph: graph object
    :rtype: list[dict]
    """
    serialized_json = graph.serialize(format='json-ld', indent=4)
    json_wp_pathway = json.loads(serialized_json.decode("utf-8"))

    return json_wp_pathway

def parse_pathway(pathway_path):
    """After importing the indicated pathway from text file resources into a graph rdflib object(import_pathway), calls the diferent data types transformations (convert_json function) and the first statement of the parser that will return a graph data structure (parse_entries function). This retrieved graph will be converted to a networX graph (convert_to_nx function).

    :param str pathway_path: pathway identifier
    :rtype: networkx.MultiDiGraph
    """
    graph = parse_rdf(pathway_path, format='turtle')

    json_wp_pathway = convert_json(graph)
    pathway_network = parse_entries(json_wp_pathway)

    nodes = pathway_network['nodes']
    interactions = pathway_network['interactions']
    pathway_info = pathway_network['pathway_info']

    return convert_to_nx(nodes, interactions, pathway_info)
