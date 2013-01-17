
import unittest
import os
import csv
import sys
import logging
import overlap
import debruijn
import networkx as nx
log = logging.getLogger("test_database")


class TestAssembleReads(unittest.TestCase):

    def setUp(self):
        self.fn_easy_case = "easy_case.txt"

    def test_easy_case(self):
        """ Test the assembly of the easy case """
        f = open(self.fn_easy_case, "r")
        solution = "O draconian devil! Oh lame saint!"
        for line in f:
            reads = line.rstrip().split(";")
            contigs = overlap.assemble(reads)
            self.assertEqual(contigs[0],solution)
        f.close()


    @unittest.skip("blah")
    def test_build_graph(self):
        """ Test building a basic deBruijn graph """
        G = debruijn.build_graph(self.fn_easy_case, 4)
        print G.nodes()
        print G.edges()
        for n in G.nodes_iter():
            print n, G.degree(n), G.neighbors(n)
#        edges = debruijn.assemble(G)
#        print "Assembled"
#        for e in edges:
#            print e
        print nx.is_eulerian(G)
#        self.assertEqual(len(G.nodes()),9)


class TestOverlap(unittest.TestCase):

    def test_overlap(self):
        """ Test the overlap function with an easy case
        """
        seq1 = "abdc sdf "
        seq2 = "sdf sabd"
        seq3 = "cc"
        n_overlaps = overlap.characters_overlapping(seq1, seq2)
        self.assertEqual(n_overlaps, 4)
        ov = seq2[0:n_overlaps]
        self.assertEqual(ov, "sdf ")
        n_overlaps = overlap.characters_overlapping(seq2, seq1)
        self.assertEqual(n_overlaps, 3)
        ov = seq1[0:n_overlaps]
        self.assertEqual(ov, "abd")
        n_overlaps = overlap.characters_overlapping(seq2, seq3)
        self.assertEqual(n_overlaps, 0)

if __name__ == "__main__":
    logging.basicConfig(stream=sys.stdout)
    logging.root.setLevel(logging.INFO)
    unittest.main()
