"""
@author Sina Beier

sb59@sanger.ac.uk

Copyright is owned by the Sanger Institute (Genome Research Ltd.) and Sina Beier

This is a script for highlighting clusters and edges on a pangenome graph based on a straingraph.

"""

import argparse
import networkx as nx
import re
import matplotlib
matplotlib.use('PS')

global panG


def main():
    """
        Main function
    """
    global panG
    panG = nx.read_graphml(args.pangenome)
    subgraph(args.subgraph)

    for n in list(panG.nodes()):
        if 'weight' not in panG.nodes[n]:
            panG.remove_node(n)
            print("node removed: " + n)
        if 'overlay' not in panG.nodes[n]:
            panG.nodes[n]['overlay'] = "gray"
    for a,b in list(panG.edges()):
        if 'overlay' not in panG[a][b]:
            panG[a][b]['overlay'] = "gray"

    nx.write_graphml(panG, args.highlighted)


def subgraph(sub):
    global panG
    subG = nx.read_graphml(sub)
    for n in subG.nodes():
        l = subG.nodes[n]['cluster']
        if panG.has_node(l):
            panG.nodes[l]['overlay'] = 'green'
        else:
            panG.add_node(l, weight=1.0, strains=1, overlay='green')
    for n1, n2 in subG.edges():
        l1 = subG.nodes[n1]['cluster']
        l2 = subG.nodes[n2]['cluster']
        if panG.has_edge(l1, l2):
            panG[l1][l2]['overlay'] = 'green'
        else:
            panG.add_edge(l1, l2, weight=1.0, overlay='green')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pangenome", help="Pangenome GraphML")
    parser.add_argument("subgraph", help="Subgraph (Pattern graph or single Strain graph) to be highlighted")
    parser.add_argument("highlighted", help="Highlighted GraphML")
    args = parser.parse_args()
    main()
