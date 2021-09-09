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
    with open(args.ideogenomes) as ids:
        for line in ids:
            l = re.sub("\n", "", line)
            ideogenome(l)

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


def ideogenome(ideo):
    global panG
    ideoG = nx.read_graphml(ideo)
    for n in ideoG.nodes():
        l = ideoG.nodes[n]['cluster']
        if panG.has_node(l):
            panG.nodes[l]['overlay'] = 'orange'
        else:
            panG.add_node(l, weight=1.0, strains=1, overlay='yellow')
    for n1, n2 in ideoG.edges():
        l1 = ideoG.nodes[n1]['cluster']
        l2 = ideoG.nodes[n2]['cluster']
        if panG.has_edge(l1, l2):
            panG[l1][l2]['overlay'] = 'orange'
        else:
            panG.add_edge(l1, l2, weight=1.0, overlay='yellow')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pangenome", help="Pangenome GraphML")
    parser.add_argument("ideogenomes", help="File with list of Ideogenome GraphMLs (one per line)")
    parser.add_argument("highlighted", help="Highlighted GraphML")
    args = parser.parse_args()
    main()
