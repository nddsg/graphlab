import numpy
from os.path import isfile

def parseFile(type, prefix, batch, nodes, modifier, directory, directed):
    filename = "%s/%s-%s-%s-%s.snap.out" % (directory, prefix, batch, nodes, modifier)
    if isfile(filename):
        file = open(filename, "r")
        time = ""
        clocktime = ""
        system = ""
        memory = ""
        voluntary = ""
        involuntary = ""
        exit = ""
        diameter = ""
        percent = ""
        for line in file:
            line = line.strip()
            if line != "":
                if line[:12] == "RuntimeError":
                    return ""
                token = line.split()
                if (token[0] == "Diameter:") and (len(token) >= 2):
                    diameter = token[1]
                if line[:27] == "The approximate diameter is":
                    diameter = token[4]
                if token[0] == "TIME:": #wall clock
                    time = token[1]
                if token[0] == "CLOCKTIME:": # user time
                    clocktime = token[1]
                if line[:11] == "System time": # system time
                    system = token[3]
                if line[:14] == "Percent of CPU":
                    percent = token[6]
                if line[:25] == "Maximum resident set size":
                    memory = token[5]
                if line[:26] == "Voluntary context switches":
                    voluntary = token[3]
                if line[:28] == "Involuntary context switches":
                    involuntary = token[3]
                if line[:33] == "Major (requiring I/O) page faults":
                    faults = token[5]
                if token[0] == "Exit":
                    exit = token[2]
        file.close()
        return "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" % (type, directed, batch, nodes, modifier, diameter, time, clocktime, system, percent, memory, voluntary, involuntary, faults, exit)
    return ""

file = open("stats.csv","w")
file.write("program,type,directed,batch,nodes,modifier,diameter,time,clocktime,system,percent,memory,voluntary,involuntary,faults,exit\n")
vertex = [10 ** i for i in range(1,6)]
edgefactor = [2 ** i for i in range(1,6)]
batch = range(5)

for v in vertex:
    for f in edgefactor:
        e = v * f
        for b in batch:
            for d in ["0", "1"]:
                output = parseFile(type="er", prefix="er%s" % ("d" if d == "1" else "u"), batch=b, nodes=v, modifier=f, directory="resultsnetworkx", directed=d)
                if output != "":
                    file.write("nx," + output + "\n")
                output = parseFile(type="er", prefix="er%s" % ("d" if d == "1" else "u"), batch=b, nodes=v, modifier=f, directory="resultsgraphlab", directed=d)
                if output != "":
                    file.write("gl," + output + "\n")
                output = parseFile(type="er", prefix="er%s" % ("d" if d == "1" else "u"), batch=b, nodes=v, modifier=f, directory="resultsgraphtool", directed=d)
                if output != "":
                    file.write("gtp," + output + "\n")
                output = parseFile(type="er", prefix="er%s" % ("d" if d == "1" else "u"), batch=b, nodes=v, modifier=f, directory="resultsgraphlabpseudo", directed=d)
                if output != "":
                    file.write("glp," + output + "\n")
                output = parseFile(type="er", prefix="er%s" % ("d" if d == "1" else "u"), batch=b, nodes=v, modifier=f, directory="resultssnap", directed=d)
                if output != "":
                    file.write("snp," + output + "\n")

vertex = [10 ** i for i in range(1,6)]
alpha = numpy.linspace(2.0, 3.0, 11)
batch = range(5)

for v in vertex:
    for a in alpha:
        for b in batch:
            output = parseFile(type="sf", prefix="sf", batch=b, nodes=v, modifier=a, directory="resultsnetworkx", directed="0")
            if output != "":
                file.write("nx," + output + "\n")
            output = parseFile(type="sf", prefix="sf", batch=b, nodes=v, modifier=a, directory="resultsgraphlab", directed="0")
            if output != "":
                file.write("gl," + output + "\n")
            output = parseFile(type="sf", prefix="sf", batch=b, nodes=v, modifier=a, directory="resultsgraphtool", directed="0")
            if output != "":
                file.write("gtp," + output + "\n")
            output = parseFile(type="sf", prefix="sf", batch=b, nodes=v, modifier=a, directory="resultsgraphlabpseudo", directed="0")
            if output != "":
                file.write("glp," + output + "\n")
            output = parseFile(type="sf", prefix="sf", batch=b, nodes=v, modifier=a, directory="resultssnap", directed="0")
            if output != "":
                file.write("snp," + output + "\n")

file.close()
