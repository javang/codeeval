import networkx as nx
import itertools
import logging
log = logging.getLogger("debruijn")

class DeBruijnGraph(nx.DiGraph):

    def __init__(self, kmer_size):
        nx.DiGraph.__init__(self)
        self.kmer_size = kmer_size

    def add_read(self, read):
        log.debug("Adding read: %s", read)
        L = len(read)
        if len(read) < self.kmer_size:
            raise ValueError("The size of the read must be smaller than the kmer_size")
        k = self.kmer_size-1
        for i in range(L - k + 1):
            log.debug("adding node %s",read[i:i+k])
            self.add_node(read[i:i+k])

    def build_edges(self):
        """
           Build all the edges of the graph. The nodes of the graph
           are k-1-mers. An edge between k-1-mers A and B is drawn if
           the end A is equal to the beginning of B
        """
        log.info("Building the edges of the graph")
        if self.number_of_nodes() == 0 :
            raise ValueError("Can not build edges in an empty graph")
        k = self.kmer_size-1


        for n1 in self.nodes_iter():
            for n2 in self.nodes_iter():
                if n1[1:k] == n2[0:k-1]:
                    log.debug("Adding \"%s\" ==> \"%s\"",n1, n2)
                    self.add_edge(n1,n2)
        """
        n = self.number_of_nodes()
        nodes = self.nodes()
        for i in range(0,n):
            for j in range(0,n):
                if nodes[i][1:k] == nodes[j][0:k-1]:
                    log.debug("Adding \"%s\" ==> \"%s\"",nodes[i], nodes[j])
                    self.add_edge(nodes[i], nodes[j])
        """

    def add_reads(self, reads):
        for read in reads:
            self.add_read(read)





def build_graph(fn, kmer_size):
    """ Build a deBruijn graph with the contents of the
        file. The file is assumed to have multiple lines,
        and each line contains various reads separated
        by an ; character
    """
    dBG = DeBruijnGraph(kmer_size)

    f = open(fn, "r")
    for line in f:
        reads = line.rstrip().split(";")
        log.debug("Reads : %s", reads)
        dBG.add_reads(reads)
    dBG.build_edges()
    f.close()
    return dBG


def assemble(graph):
    edges = nx.eulerian_circuit(graph)
    return edges
