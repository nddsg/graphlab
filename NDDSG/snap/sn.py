#!/usr/bin/python

import snap as sn
import sys
import getopt
import glob
import time


def main(argv):
    helptext = "sn.py -f <filename> -d <0, 1>"
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
        g = sn.TNGraph.New()
    else:
        g = sn.TUNGraph.New()
    verts = {}
    for file in (glob.glob(filename) + glob.glob(filename + ".*")):
        f = open(file, 'r')
        for line in f:
            if line[0] == "#":
                continue
            v1, v2 = line.strip().split("\t")[0:2]
            v1, v2 = int(v1), int(v2)
            if not verts.has_key(v1):
                verts[v1] = g.AddNode(v1)
            if not verts.has_key(v2):
                verts[v2] = g.AddNode(v2)
            g.AddEdge(v1, v2)

    start = time.time()
    cstart = time.clock()
    print "Diameter:", sn.GetBfsFullDiam(g, 1, directed)
    print "TIME: ", (time.time() - start)
    print "CLOCKTIME: ", (time.clock() - cstart)

if __name__ == "__main__":
    main(sys.argv[1:])
