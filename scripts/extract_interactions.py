import os
from io import StringIO

import pandas as pd
import requests
import rdflib
url = 'http://sparql.wikipathways.org/sparql'

query = """
SELECT DISTINCT
    ?pathwayIdentifier
    ?pathwayTitle
    ?interaction
    (GROUP_CONCAT(DISTINCT ?interactionType; separator=", ") AS ?interactionTypes)
    ?sourceNamespace
    ?sourceIdentifier
    ?sourceEntrezIdentifier
    ?sourceEntrezLabel
    ?targetNamespace
    ?targetIdentifier
    ?targetEntrezIdentifier
    ?targetEntrezLabel
WHERE {
   ?pathway a wp:Pathway .
   OPTIONAL { ?pathway dcterms:identifier ?pathwayIdentifier . }
   OPTIONAL { ?pathway dc:title ?pathwayTitle . }
   ?interaction dcterms:isPartOf ?pathway . 
   ?interaction a wp:DirectedInteraction .
   ?interaction a ?interactionType .
   ?interaction wp:source ?source .
   ?interaction wp:target ?target .
   OPTIONAL { 
        ?source wp:bdbEntrezGene ?sourceEntrez . 
        ?sourceEntrez dcterms:identifier ?sourceEntrezIdentifier .    
        ?sourceEntrez rdfs:label ?sourceEntrezLabel .
    }
   OPTIONAL { 
        ?target wp:bdbEntrezGene ?targetEntrez . 
        ?targetEntrez dcterms:identifier ?targetEntrezIdentifier .  
        ?targetEntrez rdfs:label ?targetEntrezLabel .  
   }
   OPTIONAL { ?source dc:source ?sourceNamespace . }
   OPTIONAL { ?target dc:source ?targetNamespace . }
   OPTIONAL { ?source dcterms:identifier ?sourceIdentifier . }
   OPTIONAL { ?source dcterms:identifier ?targetIdentifier . }
}
"""

def load_sar():
    directory = '/Users/cthoyt/dev/wikipathways-covid19/wp/Human'
    for name in os.listdir(directory):
        path = os.path.join(directory, name)
        graph = rdflib.Graph()
        graph.load(path, format='turtle')
        print(len(graph))
        res = graph.query(query)
        if not res:
            print(f'Failed on {name}')
            continue
        print(f'\nResults for {name}:\n')
        for row in res:
            print(row)

def main():
    res = requests.get(url, params=dict(query=query, format='text/csv'))
    df = pd.read_csv(StringIO(res.text), dtype=str)
    df.to_csv('test3.tsv', sep='\t', index=False)


if __name__ == '__main__':
    load_sar()
