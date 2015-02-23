import numpy
from os import system
from os.path import isfile
from subprocess import PIPE, Popen, STDOUT
import threading

def cmdline(command):
    process = Popen(
        args=command,
        stdout=PIPE,
        stderr=STDOUT,
        shell=True
    )
    return process.communicate()[0]
    
def getGraphLabCommand(path, filename, directed):
    return "/usr/bin/time -v python ../networkx/nx.py -f %s/%s -d %s" % (path, filename, directed)
    
commands = dict()
threads = []

vertex = [10 ** i for i in range(1,6)]
edgefactor = [2 ** i for i in range(1,6)]
batch = range(5)

for v in vertex:
    for f in edgefactor:
        e = v * f
        for b in batch:
            for d in ["0", "1"]:
                filename = "er%s-%d-%d-%d.snap" % ("d" if d == "1" else "u", b, v, f)
                commands[filename] = getGraphLabCommand("~/gdata/er", filename, d)

vertex = [10 ** i for i in range(1,6)]
alpha = numpy.linspace(2.0, 3.0, 11)
batch = range(5)

for v in vertex:
    for a in alpha:
        for b in batch:
            filename = "sf-%d-%d-%.1f.snap" % (b, v, a)
            commands[filename] = getGraphLabCommand("~/gdata/sf", filename, "0")

accessCommands = threading.Lock()
def getNextCommand():
    accessCommands.acquire()
    keys = commands.keys()
    if (len(keys) == 0):
        accessCommands.release()
        return "", ""
    key = keys[0]
    value = commands[key]
    del commands[key]
    accessCommands.release()
    return key, value
            
printLock = threading.Lock()
def donetworkx():
    while True:
        filename, command = getNextCommand()
        if (filename == ""):
            break
        resultfile = "resultsnetworkx/%s.out" % (filename)
        if (isfile(resultfile)):
            printLock.acquire()
            print "Skipping:", filename
            printLock.release()
            continue
        output = cmdline(command)
        f = open(resultfile, 'w')
        f.write(output)
        f.close()
        printLock.acquire()
        print "Completed:", filename
        printLock.release()

for n in range(50):
    thread = threading.Thread(target=donetworkx)
    thread.start()
    threads.append(thread)

for thread in threads:
    thread.join()
