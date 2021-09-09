"""
@author Sina Beier

sb59@sanger.ac.uk

Copyright is owned by the Sanger Institute (Genome Research Ltd.) and Sina Beier

This is a script for collecting statistics on the clusters.

"""

import argparse
import networkx as nx
import re
import os
import statistics

global clusters
clusters = dict()


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


def stats():
    global clusters
    with open(args.outfile, 'w') as o:
        num = len(list(clusters.keys()))
        proteins = 0
        sizes = list()
        maximal = 0
        singletons = 0
        for c in clusters:
            proteins = proteins + len(clusters[c])
            sizes.append(len(clusters[c]))
            if len(clusters[c]) >= args.maximum:
                maximal = maximal + 1
            if len(clusters[c]) == 1:
                singletons = singletons + 1
        mean = statistics.mean(sizes)
        median = statistics.median(sizes)

        o.write("Number of clusters: "+str(num)+"\n")
        o.write("Number of proteins: "+str(proteins)+"\n")
        o.write("Mean size of clusters: "+str(mean)+"\n")
        o.write("Median size of clusters: "+str(median)+"\n")
        o.write("Number of full clusters: "+str(maximal)+"\n")
        o.write("Number of singletons: " + str(singletons) + "\n")


def main():
    """
        Main function
    """
    read_clusters(args.clusters)
    stats()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("clusters", help="Cluster file from Panakeia")
    parser.add_argument("outfile", help="Output information")
    parser.add_argument("maximum", type=int, help="Size of a full cluster")
    args = parser.parse_args()
    main()
