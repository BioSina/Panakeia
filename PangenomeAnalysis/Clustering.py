"""
@author Sina Beier

sb59@sanger.ac.uk

Copyright is owned by the Sanger Institute (Genome Research Ltd.) and Sina Beier

This is a script for clustering proteins iteratively. Uses 10 threads

"""

import argparse
import subprocess
import re
import sys


global protein_dict
protein_dict = dict()

global protein_names
protein_names = dict()

global keepdict
keepdict = dict()

global fullClusters
fullClusters = dict()


def filter_clusters(clusters, outfile):
    with open(clusters) as clu:

        cluster = list()
        rep = ""
        for line in clu:
            if line.startswith('>'):
                if len(cluster) >= args.count:
                    if not rep == "":
                        fullClusters[rep] = cluster
                    for c in cluster:
                        protein_dict.pop(c)
                cluster = list()
                rep = ""
            else:
                s = re.split('>', line)
                ident = re.split(r'\.\.\.', s[-1])[0]
                ident = re.sub(" ", "", ident)
                cluster.append(ident)
                if line.endswith('*\n'):
                    rep = ident

        if len(cluster) >= args.count:
            if not rep == "":
                fullClusters[rep] = cluster
            for c in cluster:
                protein_dict.pop(c)

    with open(outfile, 'w') as o:
        for p in protein_dict:
            o.write(">" + p + "\n")
            o.write(protein_dict[p] + "\n")


def last_clusters(clusters):
    with open(clusters) as clu:
        cluster = list()
        rep = ""
        for line in clu:
            if line.startswith('>'):
                if cluster:
                    if not rep == "":
                        fullClusters[rep] = cluster

                    for c in cluster:
                        protein_dict.pop(c)
                cluster = list()
                rep = ""
            else:
                s = re.split('>', line)
                ident = re.split(r'\.\.\.', s[-1])[0]
                ident = re.sub(" ", "", ident)

                cluster.append(ident)
                if line.endswith('*\n'):
                    rep = ident
        if not rep == "":
            fullClusters[rep] = cluster
        for c in cluster:
            protein_dict.pop(c)


def seq2dict(file):
    global protein_names
    outdict = dict()
    with open(file) as f:
        name = None
        seq = ""
        for line in f:
            if line.startswith(">"):
                if name:
                    shortname = re.split(" ", name)[0]
                    protein_names[shortname] = name
                    outdict[shortname] = seq

                    seq = ""
                name = re.sub("\n", "", line)
                name = re.sub(">", "", name)
            else:
                part = re.sub("\n", "", line)
                seq = seq+part
        shortname = re.split(" ", name)[0]
        protein_names[shortname] = name
        outdict[shortname] = seq
    return outdict


def dict2file(indict, outfile):
    with open(outfile, 'w') as outf:
        for x in indict:
            outf.write(">"+x+"\n")
            outf.write(indict[x]+"\n")


def main():
    """
        Main function
    """
    global protein_dict
    global keepdict
    protein_dict = seq2dict(args.infasta)
    #sys.stderr.write(str(list(protein_dict.keys())))
    keepdict = seq2dict(args.infasta)
    # cd-hit -i allproteins.fasta -o cdhit_${ID} -c 0.7 -s 0.6 -n 5 -T 10

    c1 = subprocess.Popen(['cd-hit', '-M', str(0), '-d', str(0), '-i', args.infasta, '-o', 'cdhit90', '-c', str(0.9), '-s', str(0.6), '-n', str(5),
                           '-T', str(10)])
    c1.wait()
    filter_clusters("cdhit90.clstr", "proteins90.fasta")
    c2 = subprocess.Popen(['cd-hit', '-M', str(0), '-d', str(0), '-i', "proteins90.fasta", '-o', 'cdhit80', '-c', str(0.8), '-s', str(0.6), '-n', str(5),
                           '-T', str(10)])
    c2.wait()
    filter_clusters("cdhit80.clstr", "proteins80.fasta")
    c3 = subprocess.Popen(
        ['cd-hit', '-M', str(0), '-d', str(0), '-i', "proteins80.fasta", '-o', 'cdhit75', '-c', str(0.75), '-s', str(0.6), '-n', str(5),
         '-T', str(10)])
    c3.wait()
    filter_clusters("cdhit75.clstr", "proteins75.fasta")
    c4 = subprocess.Popen(
        ['cd-hit', '-M', str(0), '-d', str(0),  '-i', "proteins75.fasta", '-o', 'cdhit70', '-c', str(0.7), '-s', str(0.6), '-n', str(5),
         '-T', str(10)])
    c4.wait()
    last_clusters("cdhit70.clstr")
    co = 1
    clu_out = open("protein_clusters.txt", 'w')
    fa_out = "representative_seq.faa"
    newrecs = dict()
    for ident in fullClusters:
        shortname = "cluster_"+str(co)
        name = shortname+" "+protein_names[ident]
        co += 1
        cjoin = '\t'.join(fullClusters[ident])
        clu_out.write(shortname+": "+cjoin+"\n")
        newrecs[name] = keepdict[ident]
    dict2file(newrecs, fa_out)
    clu_out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infasta", help="Input file of protein sequences")
    parser.add_argument("count", type=int, help="Number of sequence origins (strains)")
    args = parser.parse_args()
    main()


print("done")
