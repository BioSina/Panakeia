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


def stats():
    proteinlist = list()
    connectedlist = list()
    longestlist = list()
    seclongestlist = list()
    with open(args.outtable, 'w') as table:
        table.write("name\tnum_proteins\tnum_connected_components\tlargest_cc\tsecond_largest_cc\n")
        for f in os.listdir(args.straingraphs):
            if f.endswith(".graphML"):
                name = re.sub(".graphML", "", f)
                graph = nx.read_graphml(args.straingraphs+"/"+f)
                proteins = len(graph.nodes())
                temp = nx.MultiGraph(graph)
                edges = list(temp.edges())
                for u, v in edges:
                    if temp[u][v][0]['color'] == "grey":
                        temp.remove_edge(u, v)
                ccs = sorted(nx.connected_components(temp), key=len, reverse=True)
                connection = len(ccs)
                largest_cc = len(ccs[0])
                if len(ccs) > 1:
                    second_largest_cc = len(ccs[1])
                else:
                    second_largest_cc = 0
                info = list()
                info.append(proteins)
                info.append(connection)
                info.append(largest_cc)
                info.append(second_largest_cc)
                table.write(name+"\t"+str(info[0])+"\t"+str(info[1])+"\t"+str(info[2])+"\t"+str(info[3])+"\t"+"\n")
                proteinlist.append(proteins)
                connectedlist.append(connection)
                longestlist.append(largest_cc)
                seclongestlist.append(second_largest_cc)
    with open(args.outfile, 'w') as o:
        o.write("Mean number of proteins: "+str(statistics.mean(proteinlist))+"\n")
        o.write("Median number of proteins: "+str(statistics.median(proteinlist))+"\n")
        o.write("Mean number of connected components: " + str(statistics.mean(connectedlist)) + "\n")
        o.write("Median number of connected components: " + str(statistics.median(connectedlist)) + "\n")
        o.write("Mean length of longest cc: "+str(statistics.mean(longestlist))+"\n")
        o.write("Median length of longest cc: "+str(statistics.median(longestlist))+"\n")
        o.write("Mean length of second longest cc: "+str(statistics.mean(seclongestlist))+"\n")
        o.write("Median length of second longest cc: "+str(statistics.median(seclongestlist))+"\n")


def main():
    """
        Main function
    """
    stats()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("straingraphs", help="Directory of (chromosomized) straingraphs")
    parser.add_argument("outtable", help="Output table")
    parser.add_argument("outfile", help="Output information")
    args = parser.parse_args()
    main()
