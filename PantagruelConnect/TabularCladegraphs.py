"""
@author Sina Beier

sb59@sanger.ac.uk

Copyright is owned by the Sanger Institute (Genome Research Ltd.) and Sina Beier

This is a script for printing tables from cladegraphs.

"""

import argparse
import networkx as nx
import re
import os


def tab_cladegraph(clade, name):
    outprefix = args.indir+"/"+name+"_"
    s = nx.read_graphml(clade)
    s2 = nx.read_graphml(clade)
    nodes = list(s.nodes())
    if len(nodes) > 0:
        for n in nodes:
            if s.node[n]['core'] == "hard core":
                s.remove_node(n)
        ccs = sorted(nx.connected_components(s), key=len, reverse=True)
        with open(outprefix+"softcore+.txt", 'w') as out:
            counter = 1
            for cc in ccs:
                components = ",".join(cc)
                annots = list()
                for c in cc:
                    annots.append(s.node[c]['annotation'])
                annotations = "|".join(annots)
                out.write("component"+str(counter)+"\t"+components+"\t"+annotations+"\n")
                counter = counter+1

        nodes = list(s.nodes())
        for n in nodes:
            if s.node[n]['core'] == "soft core":
                s.remove_node(n)
        ccs = sorted(nx.connected_components(s), key=len, reverse=True)
        with open(outprefix+"shell+.txt", 'w') as out:
            counter = 1
            for cc in ccs:
                components = ",".join(cc)
                annots = list()
                for c in cc:
                    annots.append(s.node[c]['annotation'])
                annotations = "|".join(annots)
                out.write("component"+str(counter)+"\t"+components+"\t"+annotations+"\n")
                counter = counter + 1

        nodes = list(s.nodes())
        for n in nodes:
            if s.node[n]['core'] == "shell":
                s.remove_node(n)
        ccs = sorted(nx.connected_components(s), key=len, reverse=True)
        with open(outprefix+"cloud.txt", 'w') as out:
            counter = 1
            for cc in ccs:
                components = ",".join(cc)
                annots = list()
                for c in cc:
                    annots.append(s.node[c]['annotation'])
                annotations = "|".join(annots)
                out.write("component"+str(counter)+"\t"+components+"\t"+annotations+"\n")
                counter = counter + 1

        nodes = list(s2.nodes())
        for n in nodes:
            if s2.node[n]['highlight'] != "specific":
                s2.remove_node(n)
        ccs = sorted(nx.connected_components(s2), key=len, reverse=True)
        with open(outprefix+"specific.txt", 'w') as out:
            counter = 1
            for cc in ccs:
                components = ",".join(cc)
                annots = list()
                for c in cc:
                    a = re.sub("\n", "", s2.node[c]['annotation'])
                    annots.append(a)
                annotations = "|".join(annots)
                annotations = re.sub("\n", "", annotations)
                out.write("component"+str(counter)+"\t"+components+"\t"+annotations+"\n")
                counter = counter + 1


def main():
    """
        Main function
    """
    for i in os.listdir(args.indir):
        if i.endswith("_cladegraph.graphML"):
            if not i.startswith(r'\_'):
                name = re.split(r'\.', i)[0]
                name = re.split("_", name)[0]
                tab_cladegraph(args.indir+"/"+i, name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("indir", help="Cladegraph Directory")
    args = parser.parse_args()
    main()
