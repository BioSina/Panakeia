"""
@author Sina Beier

sb59@sanger.ac.uk

Copyright is owned by the Sanger Institute (Genome Research Ltd.) and Sina Beier

This is a script for marking assigned clades  on a pangenome graph based on output from Florents analysis.

"""

import argparse
import networkx as nx
import re


def mark_pangenome():
    s = nx.read_graphml(args.pangenome)
    for n in s.nodes():
        c = s.node[n]['clades']
        clades = re.split(",", c)

        p = s.node[n]['specifics']
        if p != "none":
            specific = re.split(",", p)
        else:
            specific = []
        s.node[n]['highlight'] = "none"
        if args.clade in clades:
            s.node[n]['highlight'] = "contained"
        if args.clade in specific:
            s.node[n]['highlight'] = "specific"
    for u, v in s.edges():
        c = s[u][v]['clades']
        clades = re.split(",", c)
        s[u][v]['highlight'] = "none"
        if args.clade in clades:
            s[u][v]['highlight'] = "contained"

    nx.write_graphml(s, args.outfile)


def main():
    """
        Main function
    """
    mark_pangenome()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pangenome", help="Pangenome GraphML")
    parser.add_argument("outfile", help="Highlighted pangenome GraphML")
    parser.add_argument("clade", help="CladeID")

    args = parser.parse_args()
    main()
