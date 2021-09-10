"""
@author Sina Beier

sb59@sanger.ac.uk

Copyright is owned by the Sanger Institute (Genome Research Ltd.) and Sina Beier

This is a script for generating pangenome networks which are synteny and paralogy aware.

"""


import argparse
import networkx as nx
import re
import os
import sys

global num_genomes
num_genomes = 0

global clusters
clusters = dict()

global strains
strains = set()

global pangenome
pangenome = nx.Graph()

global allstrains
allstrains = dict()

#mapping all loci to clusters
global mapping
mapping = dict()

global annotations
annotations = dict()

global paralogs
paralogs = dict()

global parastrains
parastrains = dict()

global extraedges
extraedges = dict()

global hoods
hoods = dict()

global clusterhoods
clusterhoods = dict()


def read_clusters(cluster, rep):
    global mapping
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
            for loc in loci:
                mapping[loc] = myid
    print("Read clusters!")
    with open(rep) as r:
        for line in r:
            if line.startswith('>'):
                s = re.split(" ", line)
                myid = re.sub('>', "", s[0])
                repseq = s[1]
                product = " ".join(s[2:])
                l = (repseq, product)
                annotations[myid] = l
    print("Annotated clusters!")


def read_gff(gff, strain):
    global allstrains
    strainG = nx.MultiGraph()
    with open(gff) as gff:
        print("Reading in strain "+strain)
        lastregion = ""
        lastnode = ""
        for line in gff:
            if line.startswith("##FASTA"):
                break
            else:
                if not line.startswith("#"):
                    split = re.split('\t', line)

                    s = re.split(';', split[8])
                    mytype = split[2]
                    locus = ""
                    gene = ""
                    annotation = "hypothetical protein"
                    for i in s:
                        if i.startswith("ID"):
                            locus = re.sub("ID=", "", i)
                        if i.startswith("gene"):
                            gene = re.sub("gene=", "", i)
                        if i.startswith("product"):
                            annotation = re.sub("product=", "", i)
                    if locus != "":
                        if mytype == "CDS" and locus in mapping:
                            strainG.add_node(locus, cluster=mapping[locus], paralogs="", strain=strain, stability=0, gene=gene, annotation=annotation, color='black')
                            if not lastnode == "":
                                if split[0] == lastregion:
                                    strainG.add_edge(lastnode, locus, weight=10, color='black', type="local")
                            lastregion = split[0]
                            lastnode = locus
    for n in strainG.nodes():
        for n2 in strainG.nodes:
            if n != n2:
                if strainG.nodes[n]['cluster'] == strainG.nodes[n2]['cluster']:
                    if strainG.has_edge(n, n2):
                        if strainG[n][n2][0]['type'] == "local" and len(strainG[n][n2]) == 1:
                            strainG.add_edge(n, n2, weight=1, color='grey', type="paralogy")
                    else:
                        strainG.add_edge(n, n2, weight=1, color='grey', type="paralogy")
                    p = strainG.nodes[n]['paralogs']
                    if p == "":
                        p = n2
                    else:
                        p_set = set(re.split(",", strainG.nodes[n]['paralogs']))
                        p_set.add(n2)
                        ps = ",".join(p_set)
                        p = ps
                    strainG.nodes[n]['paralogs'] = p
                    p2 = strainG.nodes[n2]['paralogs']
                    if p2 == "":
                        p2 = n
                    else:
                        p2_set = set(re.split(",", strainG.nodes[n2]['paralogs']))
                        p2_set.add(n)
                        p2s = ",".join(p2_set)
                        p2 = p2s
                    strainG.nodes[n2]['paralogs'] = p2
    allstrains[strain] = strainG
    print("Successfully read strain " + strain)


def read_gffs(gffdir):
    global num_genomes
    global mapping
    for i in os.listdir(gffdir):
        if i.endswith(".gff"):
            if not i.startswith(r'\_'):
                name = re.split(r'\.', i)[0]
                read_gff(gffdir + "/" + i, name)
                num_genomes += 1
                strains.add(name)
    mapping = None
    print("Read GFFs!")


def draw_clusters():
    global clusters
    global num_genomes
    global pangenome
    global annotations
    global allstrains
    global paralogs
    global extraedges
    global clusterhoods
    global hoods

    # nodes have weight, strains, annotation, representative and core

    for s in allstrains:
        g = allstrains[s]
        # add the nodes
        for n in g.nodes():
            c = g.nodes[n]['cluster']
            if c in pangenome:
                s_set = set(re.split(",", pangenome.nodes[c]['strains']))
                s_set.add(s)
                st = ",".join(s_set)
                pangenome.nodes[c]['strains'] = st
                if g.nodes[n]['paralogs'] != "":
                    ps = re.split(",", g.nodes[n]['paralogs'])
                else:
                    ps = []
                # weight is here just adding up the copy numbers
                pangenome.nodes[c]['weight'] = pangenome.nodes[c]['weight']+len(ps)+1
                if len(ps)+1 > pangenome.nodes[c]['paralogy']:
                    pangenome.nodes[c]['paralogy'] = len(ps)+1
            else:
                if g.nodes[n]['paralogs'] != "":
                    ps = re.split(",", g.nodes[n]['paralogs'])
                else:
                    ps = []
                pangenome.add_node(c, strains=s, weight=len(ps)+1, paralogy=len(ps)+1, annotation="unknown", representative="unknown", core="shell")

    #annotate all nodes with the representative sequence and annotation and determine their "coreness", also average out the weight
    nl = list(pangenome.nodes())
    for node in nl:
        pangenome.nodes[node]['representative'] = annotations[node][0]
        pangenome.nodes[node]['annotation'] = annotations[node][1]
        pangenome.nodes[node]['weight'] = float(pangenome.nodes[node]['weight'])/float(num_genomes)
        ss = re.split(",", pangenome.nodes[node]['strains'])
        if float(len(ss))/float(num_genomes) >= args.hardcore:
            pangenome.nodes[node]['core'] = "hard core"
        else:
            if float(len(ss))/float(num_genomes) >= args.softcore:
                pangenome.nodes[node]['core'] = "soft core"
            else:
                if float(len(ss)) / float(num_genomes) >= args.shell:
                    pangenome.nodes[node]['core'] = "shell"
                else:
                    pangenome.nodes[node]['core'] = "cloud"

    # add the edges... and handle the paralog nodes
    # and edge has a weight and a list of strains it is found in and how core it is
    for s in allstrains:
        g = allstrains[s]
        # add the nodes
        for u, v in g.edges():
            cl1 = g.nodes[u]['cluster']
            cl2 = g.nodes[v]['cluster']
            for k in g[u][v]:
                if g[u][v][k]['type'] == "local":
                    if pangenome.has_edge(cl1, cl2):
                        pangenome[cl1][cl2]['weight'] = pangenome[cl1][cl2]['weight']+1
                        s_set = set(re.split(",", pangenome[cl1][cl2]['strains']))
                        s_set.add(s)
                        st = ",".join(s_set)
                        pangenome[cl1][cl2]['strains'] = st
                    else:
                        pangenome.add_edge(cl1, cl2, weight=1, strains=s, core="cloud")

    # add the "coreness"
    edges = list(pangenome.edges())

    for u, v in edges:
        # adapt "how core" the edge is
        if float(pangenome[u][v]['weight'])/float(num_genomes) >= args.hardcore:
            pangenome[u][v]['core'] = "hard core"
        else:
            if float(pangenome[u][v]['weight'])/float(num_genomes) >= args.softcore:
                pangenome[u][v]['core'] = "soft core"
            else:
                if float(pangenome[u][v]['weight']) / float(num_genomes) >= args.shell:
                    pangenome[u][v]['core'] = "shell"
                else:
                    pangenome[u][v]['core'] = "cloud"

        if args.weightcutoff:
            if pangenome[u][v]['weight'] < args.weightcutoff:
                pangenome.remove_edge(u, v)

    nodes = pangenome.nodes()

    #This is where the paralogs are handled
    for para in nodes:
        tt = False
        if para == "cluster_1290":
            tt = True
        if pangenome.nodes[para]['paralogy'] > 1:
            neighborhoods = list()
            clusterhood = list()
            # find set of neighborhoods (left&right local neighbor in all straingraphs)
            cls = clusters[para]
            if tt:
                sys.stderr.write("cluster : " + str(para) + " has the following members: \n")
                sys.stderr.write("cls : " + str(cls) + "\n")
            for c in cls:
                if tt:
                    sys.stderr.write("c : " + str(c) + "\n")
                for s in allstrains:

                    strainG = allstrains[s]
                    if c in strainG.nodes():

                        if tt:
                            sys.stderr.write("strain : " + str(s) + "\n")
                        neighbors = strainG.neighbors(c)
                        clus = set()
                        neighs = dict()
                        neighbors = strainG.neighbors(c)
                        for n in neighbors:

                            for j in strainG[c][n]:
                                if strainG[c][n][j]['type'] == "local":
                                    clus.add(strainG.nodes[n]['cluster'])
                                    if tt:
                                        sys.stderr.write("neighbor : " + str(n) + "\n")
                                        sys.stderr.write("neighbor cluster: " + str(strainG.nodes[n]['cluster']) + "\n")
                                    if strainG.nodes[n]['cluster'] in neighs:
                                        s_set = neighs[strainG.nodes[n]['cluster']]
                                    else:
                                        s_set = set()
                                    s_set.add(s)

                                    neighs[strainG.nodes[n]['cluster']] = s_set
                        clus = frozenset(clus)

                        if neighs:
                            neighborhoods.append(neighs)
                            if c in clusterhoods:
                                cluhood = clusterhoods[c]
                            else:
                                cluhood = dict()
                            if s in cluhood:
                                s_s2 = cluhood[s]
                            else:
                                s_s2 = set()
                            s_s2.add(clus)
                            cluhood[s] = s_s2
                            clusterhoods[c] = cluhood

            # find overlaps between neighborhoods
            neighboring = dict()

            strainset = dict()
            # the best match is a) the exact match b) a fully containing match c) the partial match with the highest edge weight

            for a in neighborhoods:
                if tt:
                    sys.stderr.write("s: "+str(list(a.keys()))+"\n")
                # handle exact matches
                a_set = set(a.keys())
                placed = False
                potential = set()
                for the_set in neighboring:

                    if a_set == the_set:
                        if tt:
                            sys.stderr.write("Already in exactly. \n")

                        d = neighboring[frozenset(a_set)]
                        for x in d.keys():
                            d[x] = d[x]+1
                        neighboring[frozenset(a_set)] = d
                        strain_s = strainset[frozenset(a_set)]
                        placed = True
                        for x in a_set:
                            bx = set()
                            for b in neighborhoods:
                                if x in b:

                                    bx = bx|b[x]
                            sstrain = set()
                            if x in strain_s:
                                if type(strain_s[x]) == str:
                                    sstrain = set()
                                    sstrain.add(strain_s[x])


                                else:
                                    sstrain = sstrain | strain_s[x]
                            strain_s[x] = sstrain|bx

                    else:
                        # handle ones that are subsets (fully containing matches)
                        for b in neighboring.keys():
                            if a_set.issubset(b):
                                potential.add(b)
                                if tt:
                                    sys.stderr.write("Is a subset \n")

                        for x in a_set:
                            x_set = set()
                            x_set.add(x)
                            if x_set.issubset(b):
                                potential.add(b)
                                if tt:
                                    sys.stderr.write("Is a partial subset \n")
                # find best potential match
                oldsum = 0
                current = None

                if len(potential) > 0:
                    if tt:
                        sys.stderr.write("Finding best potential match in potential: " + str(potential) + " \n")
                    for pot in potential:
                        for x in a_set:
                            numb = 0
                            if x in pot:
                                numb = numb + neighboring[frozenset(pot)][x]
                                if numb > oldsum:
                                    current = set(pot)
                                    oldsum = numb
                    if current:
                        mydict = neighboring[frozenset(current)]
                        mystrains = strainset[frozenset(current)]

                        newcurrent = current.copy()

                        for x in a_set:
                            if x in current:
                                bx = set()
                                for b in neighborhoods:
                                    if x in b:
                                        bx = bx|b[x]
                                sstrain = set()
                                if x in strain_s:

                                    if type(strain_s[x]) == str:
                                        sstrain = set()
                                        sstrain.add(strain_s[x])

                                    else:
                                        sstrain = sstrain | strain_s[x]

                                mystrains[x] = sstrain | bx

                            else:

                                bx = set()
                                for b in neighborhoods:
                                    if x in b:
                                        bx = bx|b[x]
                                mystrains[x] = bx

                        for x in a_set:
                            if x in current:
                                mydict[x] = mydict[x]+1
                            else:
                                mydict[x] = 1

                                newcurrent.add(x)

                        neighboring.pop(frozenset(current))
                        if tt:
                            sys.stderr.write("Found best placement. \n")
                        neighboring[frozenset(newcurrent)] = mydict
                        placed = True
                        strainset.pop(frozenset(current))
                        strainset[frozenset(newcurrent)] = mystrains

                if not placed:
                    if tt:
                        sys.stderr.write("Not placed.\n")
                    d = dict()
                    for x in a_set:
                        d[x] = 1

                    neighboring[frozenset(a_set)] = d

                    strain_s = dict()
                    for x in a_set:
                        for ax_set in neighborhoods:
                            values = set.union(*ax_set.values())

                            if a[x] == values:
                                strain_s[x] = values
                    if tt:
                        sys.stderr.write("a_set: " + str(a_set) + "\n")
                    strainset[frozenset(a_set)] = strain_s

                # now I have a dictionary which is my neighborhood and includes weights to the existing neighbors
                # I need to create len(neighboring) subclusters and assign one neighborhood to each of them
            count = 1
            paradict = dict()
            straindict = dict()

            for k in neighboring:
                name = para +"_subcluster"+str(count)
                paradict[name] = k
                straindict[name] = strainset[k]
                count = count+1
                if para in hoods:
                    neighborhood = hoods[para]
                else:
                    neighborhood = dict()
                neighborhood[k] = name
                hoods[para] = neighborhood
            paralogs[para] = paradict
            parastrains[para] = straindict
            # For the neighbors which have paralogs I need to somehow handle that problem as well.
            # My current idea is to run through the list of paralog neighborhoods and adapt them.
            # Then run through it again and add all the nodes and edges
    for paral in paralogs.keys():
        testing = False
        #if paral == "cluster_1290":
            #testing = True
        if testing:
            sys.stderr.write("paralogs[paral]: " + str(paralogs[paral]) + " \n")
        for s in paralogs[paral]:
            if testing:
                sys.stderr.write("s: " + str(s) + " \n")
            for neigh in paralogs[paral][s]:
                if testing:
                    sys.stderr.write("neigh: " + str(neigh) + " \n")
                if re.search("_subcluster", neigh) or pangenome.nodes[neigh]['paralogy'] != 1:
                    if testing:
                        sys.stderr.write("Neighbor is a paralog! \n")
                    nneigh = re.split("_subcluster", neigh)[0]
                    for subclust in paralogs[nneigh]:

                        #Find the matching neighborhood and adapt it with the correct subcluster
                        # Find out if this is the neighboring subcluster which has the current paralog as neighbor,
                        # if so, replace the names (of the neighbor in current and the current paralogous subluster s in the neighboring subcluster
                        if paral in paralogs[nneigh][subclust]:
                            if testing:
                                sys.stderr.write("s: " + str(s) + " \n")
                                sys.stderr.write("subclust: " + str(subclust) + " \n")
                                sys.stderr.write("paralogs[neigh][subclust]: " + str(paralogs[neigh][subclust]) + " \n")
                            if testing:
                                sys.stderr.write("s: " + str(s) + " \n")
                            s2 = set(paralogs[paral][s])

                            if nneigh in s2:
                                if testing:
                                    sys.stderr.write("Current neighbor: " + str(nneigh) + " \n")
                                s2.remove(nneigh)
                                s2.add(subclust)
                                paralogs[paral][s] = frozenset(s2)
                                if testing:
                                    sys.stderr.write("Replaced neighbor with: "+str(subclust)+" \n")



                            n2 = set(paralogs[nneigh][subclust])

                            if paral in n2:
                                n2.remove(paral)
                                n2.add(s)
                                paralogs[nneigh][subclust] = frozenset(n2)
                                if testing:
                                    sys.stderr.write("Replaced original paralog with: "+str(s)+" \n")
                    else:
                        if testing:
                            sys.stderr.write("Neighbor is not a paralog! \n")

            # generate a subcluster node for each neighborhood and add the correct neighbors and edges
            # This is difficult if the neighbor is also a paralog!!
            # then remove the "original" node
    extraweights = dict()
    for paral in paralogs.keys():
        testing = False
        #if paral == "cluster_1290":
            #testing = True
        if testing:
            sys.stderr.write("paralogs[paral]: " + str(paralogs[paral]) + " \n")
        edgesets = set()
        for ss in paralogs[paral]:
            if testing:
                sys.stderr.write("ss: " + str(ss) + " \n")
            # There are clusters, not strains in the strainlists??
            # Also, they are empty a lot, but not always for subclusters.
            if re.search("_subcluster", ss):
                strainl = parastrains[paral][ss]
                if testing:
                    sys.stderr.write("subcluster strainl: " + str(strainl) + "\n")

            else:
                strainl = strainset[ss]
                if testing:
                    sys.stderr.write("cluster strainl: " + str(strainl) + "\n")

            if type(strainl.values()) == str:
                strainlist = set()
                strainlist.add(strainl.values())
            else:

                strainlist = set.union(*strainl.values())

            w = float(len(strainlist)) / float(num_genomes)

            straintxt = ",".join(strainlist)
            if ss not in pangenome.nodes():
                pangenome.add_node(ss, strains=straintxt, weight=w, paralogy=pangenome.nodes[paral]['paralogy'],
                                    annotation=pangenome.nodes[paral]['annotation'],
                                    representative=pangenome.nodes[paral]['representative'],
                                    core=pangenome.nodes[paral]['core'])
                if testing:
                    sys.stderr.write("Subcluster node:" + ss + " added.\n")
            else:
                pangenome.nodes[ss]['strains'] = straintxt
                pangenome.nodes[ss]['weight'] = w
                pangenome.nodes[ss]['paralogy'] = pangenome.nodes[paral]['paralogy']
                pangenome.nodes[ss]['annotation'] = pangenome.nodes[paral]['annotation']
                pangenome.nodes[ss]['representative'] = pangenome.nodes[paral]['representative']
                pangenome.nodes[ss]['core'] = pangenome.nodes[paral]['core']

            for neigh in paralogs[paral][ss]:
                if testing:
                    sys.stderr.write("neigh: " + str(neigh) + " \n")
                # Maybe I can add the neighbors without additional info and then add it later when it's their part?
                # Because any nonexisting neighbor would be another paralog.
                if neigh not in pangenome.nodes():
                    pangenome.add_node(neigh, strains="", weight=0,
                                       paralogy=0,
                                       annotation="undefined",
                                       representative="undefined",
                                       core="unknown")

                # I don't know how to get the edge information yet. I will test the rest first.
                edgeset = set()
                edgeset.add(ss)
                edgeset.add(neigh)

                edgesets.add(frozenset(edgeset))

        #sys.stderr.write("edgesets: " + str(edgesets) + " \n")

        for ex in edgesets:
            testing = False
            #for exx in ex:
                #if re.match("cluster_1290", exx):
                    #testing = True
            if testing:
                sys.stderr.write("ex: " + str(ex) + " \n")
            exset = list(ex)

            if len(exset) < 2:
                exset.append(exset[0])
            if testing:
                sys.stderr.write("exset " + str(exset) + "\n")

            olds = re.split("_subcluster", exset[0])[0]
            oldneigh = re.split("_subcluster", exset[1])[0]
            if testing:
                sys.stderr.write("olds: " + str(olds) + " \n")
                sys.stderr.write("oldneigh: " + str(oldneigh) + " \n")

            w = 0
            for st in allstrains:
                e1 = set()
                e2 = set()
                if testing:
                    sys.stderr.write("Teststrain: " + str(st) + " \n")
                s = allstrains[st]
                for n in s.nodes():
                    if s.nodes[n]['cluster'] == olds:
                        e1.add(n)
                    if s.nodes[n]['cluster'] == oldneigh:
                        e2.add(n)
                for a in e1:
                    for b in e2:
                        if testing:
                            sys.stderr.write("Edge: " + str(a) + ", " +str(b) + " \n")
                            sys.stderr.write("Strain: " + str(st) + " \n")
                        if s.has_edge(a, b):
                            if testing:
                                sys.stderr.write("e1_a: " + str(a) + " \n")
                                sys.stderr.write("e2_b: " + str(b) + " \n")
                            w = w+1
                            if ex in extraedges:
                                strainsets = extraedges[ex]
                            else:
                                strainsets = set()
                            strainsets.add(st)
                            extraedges[ex] = strainsets
            extraweights[ex] = w

        pangenome.remove_node(paral)

    for xx in extraedges:
        xeset = list(xx)
        if len(xeset) < 2:
            xeset.append(xeset[0])
        w = extraweights[xx]
        tests = False
        for tx in xeset:
            if re.match("cluster_960", tx):
                tests = True
        edgestrains = ",".join(extraedges[xx])

        if float(w) / float(num_genomes) >= args.hardcore:
            c = "hard core"
        else:
            if float(w) / float(num_genomes) >= args.softcore:
                c = "soft core"
            else:
                if float(w) / float(num_genomes) >= args.shell:
                    c = "shell"
                else:
                    c = "cloud"

        if args.weightcutoff:
            if w > args.weightcutoff:
                pangenome.add_edge(xeset[0], xeset[1], weight=w, strains=edgestrains, core=c)
                if tests:
                    sys.stderr.write("added edge for: " + str(xeset) + "\n")

        else:
            pangenome.add_edge(xeset[0], xeset[1], weight=w, strains=edgestrains, core=c)
            if tests:
                sys.stderr.write("added edge for: " + str(xeset) + "\n")

    # clean up
    thenodes = list(pangenome.nodes())
    for tn in thenodes:
        if not pangenome.nodes[tn].keys():
            pangenome.remove_node(tn)
            sys.stderr.write("Removed node: "+tn+" \n")
        if re.search("_subcluster", tn):
            degree = nx.degree(pangenome, tn)
            if degree == 0:
                #pangenome.remove_node(tn)
                sys.stderr.write("Would remove node for unfilled subcluster: " + tn + " \n")
    draw_genomes()
    pl = set(pangenome.nodes())
    for p in pl:
        s_set = set(re.split(",", pangenome.nodes[p]['strains']))
        pangenome.nodes[p]['straincount'] = len(s_set)

    print("Generated pangenome!")

    nx.write_graphml(pangenome, args.output + "/" + "pangenome.graphML")
    pangenome = None
    print("Saved pangenome!")


def draw_genomes():
    global pangenome
    global allstrains
    global clusterhoods
    global hoods

    thestrains = list(allstrains.keys())

    for s in thestrains:
        straingraph = allstrains[s]
        for n in straingraph.nodes():
            cl = straingraph.nodes[n]['cluster']

            if cl in pangenome.nodes():

                sts = re.split(",", pangenome.nodes[cl]['strains'])

            else:
                for cx in pangenome.nodes():
                    if cx.startswith(cl):
                        teststrains = re.split(",", pangenome.nodes[cx]['strains'])
                        if s in teststrains:
                            sts = teststrains
            straingraph.nodes[n]['stability'] = len(sts)
            if n in clusterhoods:


                if s in clusterhoods[n]:
                    handled = False
                    for cc in clusterhoods[n][s]:

                        #find the neighborhood in the pangenome paralogs
                        for p in hoods:
                            if straingraph.nodes[n]['cluster'].startswith(p):
                                for h in hoods[p]:
                                    if set(h) == set(cc) or set(cc).issubset(set(h)):
                                        handled = True
                                        straingraph.nodes[n]['cluster'] = hoods[p][h]
                        if not handled:
                            sys.stderr.write("node not handled: " + str(n) + "\n")

                            #sys.stderr.write("cc is from strain: " + str(s) + "\n")
                            #sys.stderr.write("cc is from node: " + str(n) + "\n")

        nx.write_graphml(straingraph, args.output+"/"+"strain_"+s+".graphML")
        allstrains.pop(s)
    allstrains = None


def main():
    """
        Main function
    """
    # Run iterative clustering on all proteins (Clustering.py)
    read_clusters(args.clust, args.rep)
    print("Read clusters!")
    # Read (annotated) gffs
    read_gffs(args.gffs)
    draw_clusters()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", help="Output directory")
    parser.add_argument("gffs", help="Directory with input gffs")
    parser.add_argument("clust", help="File of clustered proteins")
    parser.add_argument("rep", help="File of (annotated) representative sequences")
    parser.add_argument("-hardcore", type=float, default=0.99, help="Percentage of genomes to count as hard core, defaults to 0.99")
    parser.add_argument("-softcore", type=float, default=0.95, help="Percentage of genomes to count as soft core, defaults to 0.95")
    parser.add_argument("-shell", type=float, default=0.15, help="Percentage of genomes to count as shell, defaults to 0.15, everything less is cloud")
    parser.add_argument("-weightcutoff", type=int, help="Integer cutoff for minimal weight of an edge to be drawn")
    args = parser.parse_args()
    main()


print("done")
