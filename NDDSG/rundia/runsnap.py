import numpy
from os import system
from os.path import isfile
from subprocess import PIPE, Popen, STDOUT

def cmdline(command):
    process = Popen(
        args=command,
        stdout=PIPE,
        stderr=STDOUT,
        shell=True
    )
    return process.communicate()[0]
    
def getGraphLabCommand(path, filename, directed):
    return "/usr/bin/time -v python ~/diameter/snap/sn.py -f %s/%s -d %s" % (path, filename, directed)

vertex = [10 ** i for i in range(1,6)]
edgefactor = [2 ** i for i in range(1,6)]
batch = range(5)

for v in vertex:
    for f in edgefactor:
        e = v * f
        for b in batch:
            for d in ["0", "1"]:
                filename = "er%s-%d-%d-%d.snap" % ("d" if d == "1" else "u", b, v, f)
                outfile = "resultssnap/%s.out" % (filename)
                if isfile(outfile):
                    print "Skipping:", filename
                    continue
                output = cmdline(getGraphLabCommand("~/gdata/er", filename, d))
                file = open(outfile, 'w')
                file.write(output)
                file.close()
                print "Completed:", filename

vertex = [10 ** i for i in range(1,6)]
alpha = numpy.linspace(2.0, 3.0, 11)
batch = range(5)

for v in vertex:
    for a in alpha:
        for b in batch:
            filename = "sf-%d-%d-%.1f.snap" % (b, v, a)
            outfile = "resultssnap/%s.out" % (filename)
            if isfile(outfile):
                print "Skipping:", filename
                continue
            output = cmdline(getGraphLabCommand("~/gdata/sf", filename, "0"))
            file = open(outfile, 'w')
            file.write(output)
            file.close()
            print "Completed:", filename
