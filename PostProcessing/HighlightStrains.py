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
        if 'weight' not in panG.node[n]:
            panG.remove_node(n)
            print ("node removed: " + n)
        if 'overlay' not in panG.node[n]:
            panG.node[n]['overlay'] = panG.node[n]['color']
    for a,b in list(panG.edges()):
        if 'overlay' not in panG[a][b]:
            panG[a][b]['overlay'] = panG[a][b]['color']

    nx.write_graphml(panG, args.highlighted)


def ideogenome(ideo):
    global panG
    ideoG = nx.read_graphml(ideo)
    for n in ideoG.nodes():
        l = ideoG.nodes[n]['label']
        if panG.has_node(l):
            panG.nodes[l]['overlay'] = 'orange'
        else:
            panG.add_node(l, weight=1.0, strains=1, color='yellow')
    for n1, n2 in ideoG.edges():
        l1 = ideoG.nodes[n1]['label']
        l2 = ideoG.nodes[n2]['label']
        if panG.has_edge(l1, l2):
            panG[l1][l2]['overlay'] = 'orange'
        else:
            panG.add_edge(l1, l2, weight=1.0, color='yellow')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pangenome", help="Pangenome GraphML")
    parser.add_argument("ideogenomes", help="File with list of Ideogenome GraphMLs (one per line)")
    parser.add_argument("highlighted", help="Highlighted GraphML")
    args = parser.parse_args()
    main()
