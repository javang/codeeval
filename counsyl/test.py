
import unittest
import os
import csv
import sys
import logging
import KMP
import debruijn
import overlap
# import networkx as nx
log = logging.getLogger("test_database")


class TestOverlap(unittest.TestCase):

    def test_overlap(self):
        """ Test the overlap function with an easy case
        """
        seq1 = "abdc sdf "
        seq2 = "sdf sabd"
        seq3 = "cc"
        n_overlaps = KMP.characters_overlapping(seq1, seq2)
        self.assertEqual(n_overlaps, 4)
        ov = seq2[0:n_overlaps]
        self.assertEqual(ov, "sdf ")
        n_overlaps = KMP.characters_overlapping(seq2, seq1)
        self.assertEqual(n_overlaps, 3)
        ov = seq1[0:n_overlaps]
        self.assertEqual(ov, "abd")
        n_overlaps = KMP.characters_overlapping(seq2, seq3)
        self.assertEqual(n_overlaps, 0)

    def test_find_overlap(self):
        """ Test the find/overlap function with an easy case
        """
        seq1 = "abdc sdf"
        seq2 = "dfke "
        seq3 = "dfs abd"
        seq4 = "dc s"
        seq5 = "pp"
        position, n_overlaps  = KMP.find_or_overlap(seq1, seq2)
        self.assertEqual(n_overlaps, 2)
        self.assertEqual(position, 6)
        position, n_overlaps = KMP.find_or_overlap(seq3, seq1)
        self.assertEqual(n_overlaps, 3)
        self.assertEqual(position, 4)
        position, n_overlaps = KMP.find_or_overlap(seq1, seq4)
        self.assertEqual(n_overlaps, 4)
        self.assertEqual(position, 2)
        position, n_overlaps = KMP.find_or_overlap(seq1, seq5)
        self.assertEqual(n_overlaps, 0)
        self.assertEqual(position, len(seq1))


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


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stdout)
    logging.root.setLevel(logging.INFO)
    unittest.main()
