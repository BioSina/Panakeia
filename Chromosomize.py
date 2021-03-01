"""
@author Sina Beier

sb59@sanger.ac.uk

Copyright is owned by the Sanger Institute (Genome Research Ltd.) and Sina Beier

This is a script for highlighting clusters and edges on a pangenome graph based on a straingraph.

"""

import argparse
import networkx as nx
import re
import os

global clusters
clusters = dict()

global chromosomes
chromosomes = dict()


def read_clusters(cluster):
    global clusters
    with open(cluster) as clu:
        for line in clu:
            line = re.sub('\n', "", line)
            ls = re.split('\t', line)
            preid = re.split(":", ls[0])
            left_locus = re.sub('^ ', "", preid[1])
            myid = preid[0]
            loci = ls[1:]
            loci.insert(0, left_locus)
            clusters[myid] = loci
    print("Read clusters!")


def template(template):
    global clusters
    global chromosomes
    temp = nx.read_graphml(template)
    temp = nx.MultiGraph(temp)
    edges = list(temp.edges())
    for u, v in edges:
        if temp[u][v][0]['color'] == "grey":
            temp.remove_edge(u, v)
    ccs = sorted(nx.connected_components(temp), key=len, reverse=True)
    print("Detected "+str(len(ccs))+" chromosomes.")
    c = 0
    for comp in ccs:
        c = c+1
        chromosom = "chr"+str(c)
        for protein in comp:
            for item in clusters:
                if protein in clusters[item]:
                    if item in chromosomes:
                        if chromosomes[item] != chromosom:
                            chromosomes[item] = "undecided"
                    else:
                        chromosomes[item] = chromosom
    print("Read templates!")


def chromosomize_strain(strain):
    global clusters
    global chromosomes
    s = nx.read_graphml(args.straingraphs+"/"+strain)
    for node in s.nodes():
        chrom = "unknown"
        for c in clusters:
            if node in clusters[c]:
                if c in chromosomes:
                    chrom = chromosomes[c]
        s.node[node]['chromosome'] = chrom
    nx.write_graphml(s, args.outdir + "/" + strain)


def main():
    """
        Main function
    """
    read_clusters(args.clusters)
    with open(args.template) as t:
        for line in t:
            l = re.sub("\n","", line)
            template(args.straingraphs+"/"+l)
    for f in os.listdir(args.straingraphs):
        if f.endswith(".graphML"):
            chromosomize_strain(f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("straingraphs", help="Directory with Straingraph GraphMLs")
    parser.add_argument("outdir", help="Directory with Chromosomized Straingraph GraphMLs")
    parser.add_argument("clusters", help="Cluster file from Panakeia")
    parser.add_argument("template", help="File with names of templates in straingraph directory")
    args = parser.parse_args()
    main()
