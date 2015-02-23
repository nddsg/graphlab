#!/usr/bin/python

import networkx as nx
import sys
import getopt
import glob
import time


def main(argv):
    helptext = "nx.py -f <filename> -d <0, 1>"
    filename = ""
    directed = False
    try:
        opts, args = getopt.getopt(argv, "hf:d:", ["file=", "directed="])
    except getopt.GetoptError:
        print helptext
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print helptext
            sys.exit()
        elif opt in ("-f", "--file"):
            filename = arg
        elif opt in ("-d", "--directed"):
            directed = int(arg) != 0
    if (filename == ""):
        print helptext
        sys.exit(2)

    if (directed):
        g = nx.DiGraph()
        base = nx.DiGraph()
    else:
        g = nx.Graph()
        base = nx.Graph()
    for file in (glob.glob(filename) + glob.glob(filename + ".*")):
        g.add_edges_from(nx.read_edgelist(file, create_using=base).edges())

    print "Original:"
    print nx.info(g)

    if (directed):
        h = [x for x in nx.strongly_connected_component_subgraphs(g)]
    else:
        h = [x for x in nx.connected_component_subgraphs(g)]
    print "Largest Connected Component:"
    print nx.info(h[0])
    start = time.time()
    cstart = time.clock()
    print "Diameter:", nx.diameter(h[0])
    print "TIME: ", (time.time() - start)
    print "CLOCKTIME: ", (time.clock() - cstart)

if __name__ == "__main__":
    main(sys.argv[1:])
