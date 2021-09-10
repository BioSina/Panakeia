"""
@author Sina Beier

sb59@sanger.ac.uk

Copyright is owned by the Sanger Institute (Genome Research Ltd.) and Sina Beier

This is a script for highlighting clades  on a pangenome graph based on output from Florents analysis.

"""

import argparse
import networkx as nx
import re
import sys

global clusters
clusters = dict()

global clades
clades = dict()

global mapping
mapping = dict()

global specific
specific = dict()


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


def read_map(map):
    global mapping
    with open(map) as m:
        for line in m:
            line = re.sub('\n', "", line)
            ls = re.split('\t', line)
            my_id = re.split("\.", ls[0])[0]
            mapping[ls[1]] = my_id
    print("Read mapping file!")


def read_specific(myspec):
    global specific
    with open(myspec) as sp:
        my_genes = set()
        for line in sp:
            line = re.sub('\n', "", line)
            if line.startswith("# clade"):
                if len(my_genes) != 0:
                    specific[name] = list(my_genes)
                name = re.split(" ", line)[1]
                name = re.sub("\n", "", name)
                my_genes = set()
            if line.startswith("# no specific gene found"):
                my_genes = set()
                specific[name] = list(my_genes)
                continue
            if not line.startswith("#"):
                if not line.startswith("gene_family_id"):
                    split = re.split("\t", line)
                    pre_gene = split[4]
                    for cl in clusters:
                        if pre_gene in clusters[cl]:
                            my_genes.add(cl)

    print("Read in clade specific genes!")


def read_clades(myclades):
    global clades
    global mapping
    with open(myclades) as cl:
        cl.readline()
        for line in cl:
            line = re.sub('\n', "", line)
            ls = re.split('\t', line)
            name = ls[0]
            name = re.sub(" ", "", name)
            name = re.sub("\n", "", name)
            clade = re.split(",", ls[1])
            my_clade = list()
            for c in clade:
                my_clade.append(mapping[c])
            clades[name] = my_clade
    print("Read in clades!")


def cladeify_pangenome():
    global clades
    global specific
    outfile = args.outdir+"/pangenome_clades.graphML"
    s = nx.read_graphml(args.pangenome)
    for node in s.nodes():
        newnode = re.split("_subcluster", node)[0]
        strains = re.split(",", s.nodes[node]['strains'])
        node_clades = set()
        node_specifics = set()
        for strain in strains:
            for c in clades:
                if strain in clades[c]:
                    node_clades.add(c)
        node_clade = ",".join(node_clades)
        s.nodes[node]['clades'] = node_clade
        for speci in specific:
            if newnode in specific[speci]:
                node_specifics.add(speci)
        if len(node_specifics) > 0:
            node_spec = ",".join(node_specifics)
        else:
            node_spec = "none"

        s.nodes[node]['specifics'] = node_spec

    for u, v in s.edges():
        strains = re.split(",", s[u][v]['strains'])
        edge_clades = set()
        for strain in strains:
            for c in clades:
                if strain in clades[c]:
                    edge_clades.add(c)
        edge_clade = ",".join(edge_clades)
        s[u][v]['clades'] = edge_clade
    nx.write_graphml(s, outfile)
    return s


def mark_pangenome(pangenome, clade):
    outfile = args.outdir+"/"+clade+"_highlighted.graphML"
    outfile2 = args.outdir+"/"+clade+"_cladegraph.graphML"

    s = pangenome.copy()
    for n in s.nodes():
        c = s.nodes[n]['clades']
        clades2 = re.split(",", c)

        p = s.nodes[n]['specifics']
        if p != "none":
            specific2 = re.split(",", p)
        else:
            specific2 = []

        s.nodes[n]['highlight'] = "none"

        if clade in clades2:
            s.nodes[n]['highlight'] = "contained"
        if clade in specific2:
            s.nodes[n]['highlight'] = "specific"

    for u, v in s.edges():
        c = s[u][v]['clades']
        clades3 = re.split(",", c)
        s[u][v]['highlight'] = "none"
        if clade in clades3:
            s[u][v]['highlight'] = "contained"

    nx.write_graphml(s, outfile)

    hardcore = 0
    softcore = 0
    shell = 0
    cloud = 0
    specif = 0

    nodes = list(s.nodes())

    for n in nodes:
        if s.nodes[n]['highlight'] == "none":
            s.remove_node(n)
        else:
            if s.nodes[n]['core'] == "hard core":
                hardcore = hardcore+1
            if s.nodes[n]['core'] == "soft core":
                softcore = softcore+1
            if s.nodes[n]['core'] == "shell":
                shell = shell+1
            if s.nodes[n]['core'] == "cloud":
                cloud = cloud+1
            if s.nodes[n]['highlight'] == "specific":
                specif = specif+1

    edges = list(s.edges())

    for u, v in edges:
        if s[u][v]['highlight'] == "none":
            s.remove_edge(u, v)

    proteins = str(len(s.nodes()))
    connections = str(len(s.edges()))

    nx.write_graphml(s, outfile2)

    s_spec = s.copy()
    nodes2 = list(s_spec.nodes())
    for n2 in nodes2:
        if s_spec.nodes[n2]['highlight'] != "specific":
            s_spec.remove_node(n2)
    ccs = sorted(nx.connected_components(s_spec), key=len, reverse=True)
    spec_long = 0
    if len(ccs) > 0:
        spec_long = len(ccs[0])

    s_cloud = s.copy()
    nodes3 = list(s_cloud.nodes())
    for n3 in nodes3:
        if s_cloud.nodes[n3]['core'] != "cloud":
            s_cloud.remove_node(n3)
    ccs = sorted(nx.connected_components(s_cloud), key=len, reverse=True)
    cloud_long = 0
    if len(ccs) > 0:
        cloud_long = len(ccs[0])

    s_shell = s.copy()
    nodes4 = list(s_shell.nodes())
    for n4 in nodes4:
        if s_shell.nodes[n4]['core'] != "shell":
            s_shell.remove_node(n4)
    ccs = sorted(nx.connected_components(s_shell), key=len, reverse=True)
    shell_long = 0
    if len(ccs) > 0:
        shell_long = len(ccs[0])
    # stats: clade  proteins    edges   hardCore    softCore    shell   cloud  specificProteins maxSynSpecific maxSynCloud maxSynShell
    stats = list()
    stats.append(clade)
    stats.append(proteins)
    stats.append(connections)
    stats.append(str(hardcore))
    stats.append(str(softcore))
    stats.append(str(shell))
    stats.append(str(cloud))
    stats.append(str(specif))
    stats.append(str(spec_long))
    stats.append(str(cloud_long))
    stats.append(str(shell_long))
    stat = "\t".join(stats)

    return stat



def main():
    """
        Main function
    """
    read_map(args.map)
    read_clusters(args.clusters)
    read_clades(args.clades)
    read_specific(args.spec)
    if args.cladepan:
        cladified = nx.read_graphml(args.cladepan)
    else:
        cladified = cladeify_pangenome()

    with open(args.outdir+"/statistics.txt", 'w') as o:
        o.write("clade\tproteins\tedges\thardCore\tsoftCore\tshell\tcloud\tspecificProteins\tmaxSynSpecific\tmaxSynCloud\tmaxSynShell\n")
        for c in clades:
            m = mark_pangenome(cladified, c)
            o.write(m+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pangenome", help="Pangenome GraphML")
    parser.add_argument("outdir", help="Output directory")
    parser.add_argument("clusters", help="Cluster file from Panakeia")
    parser.add_argument("map", help="Mapping file for Vibrio")
    parser.add_argument("clades", help="Clade file")
    parser.add_argument("spec", help="Clade specific genes file")
    parser.add_argument("-cladepan", help="Exisiting cladified pangenome file")

    args = parser.parse_args()
    main()
