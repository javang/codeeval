
import unittest
import os
import csv
import sys
import logging
import debruijn
import networkx as nx
log = logging.getLogger("test_database")

class TestAssembleReads(unittest.TestCase):

    def setUp(self):
        self.fn_easy_case = "easy_case.txt"

    @unittest.skip("blah")
    def test_easy_case(self):
        """ Test the assembly of the easy case """
        text = debruijn.assemble(self.fn_easy_case)
        solution = "O draconian devil! Oh lame saint!"
        self.assertEqual(text, solution, msg="The assembly is incorrect")

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

if __name__ == "__main__":
    logging.basicConfig(stream=sys.stdout)
    logging.root.setLevel(logging.DEBUG)
    unittest.main()

