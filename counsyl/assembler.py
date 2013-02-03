
import sys
import os
import time
import logging
import overlap
log = logging.getLogger("mlearning")


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description="""The Da Vyncy assembler """)
    parser.add_argument("fn",
                    help="File with the set of examples to assemble. One example " \
                    "per line. The fragments need to be separated by ;")
    parser.add_argument("--log",
                    dest="log",
                    default = False,
                    help="Log file")
    args = parser.parse_args()
    if(args.log):
        logging.basicConfig(filename=args.log, filemode="w")
    else:
        logging.basicConfig(stream=sys.stdout)
    logging.root.setLevel(logging.ERROR)

    f = open(args.fn, "r")
    for line in f:
        fragments = line.rstrip().split(";")
        contigs = overlap.assemble(fragments)
        for contig in contigs:
            print contig
    f.close()

