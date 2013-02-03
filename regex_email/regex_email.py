
import sys
import re
import heapq

if __name__  == "__main__":
    if len(sys.argv) != 2:
        print "Parameters"
        print "[1] file"
        quit()

    mail = re.compile(r"\w+@\w+\.\w+")
    f = open(sys.argv[1], "r")
    for line in f.readlines():
        m = re.match(mail, line)
        if m:
            print "true"
        else:
            print "false"
    f.close()
