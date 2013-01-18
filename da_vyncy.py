#!/usr/bin/python

import sys
import os
import time
import logging
import itertools
log = logging.getLogger("assembler")




"""
COUNSYL CHALLENGE 17 Jan 2013

JAVIER VELAZQUEZ-MURIEL 
(415) 283 7085
javi.velazquez@gmail.com

Algorithm for assemblying a set of text fragments.

At each iteration, the two fragments that have the maximum overlap are combined
and removed from the set. Ties are broken using a greedy approach. The algorithm
ends when all the fragments have been assembled or there are no more overlaps.
The output is the set of senteces that could be assembled from the fragments 
(There can be more than one).

The algorithm tries to avoid calculating overlaps between fragments when they
are not required. To do so, sorts the fragments by length and starts exploring
overlaps between them. As longer fragments can potentially have longer overlaps,  
there are more chances of finding a long overlap and avoid further calculations
with the shorter fragments. Additionally, a dictionary of overlaps calculated 
so far is kept. 

To save memory, after combining two fragments the entries in the dictionary
involving those fragments are deleted. 


There is no limit on the size of the fragments


Bonus points question:

A typical example of the case where a greedy (or random) approach for merging
fragments can change the final reconstruction is the presence of repeats.
For example, let's say that we have 3 fragments:

1) "bioinformatics is very very"
2) "very very very very"
3) "very very not very very very very cool"

Fragments 1 and 2 have an overlap of 9 characters "very very". The same is true
for fragments 1 and 3. If we choose to assemble in the order: 1-2-3 we would get:

1-2:   "bioinformatics is very very very very"
1-2-3: "bioinformatics is very very very very not very very very very cool"

On the other hand, if we choose to assemble as 1-3-2:

1-3:   "bioinformatics is very very not very very very very cool"
1-3-2: "bioinformatics is very very not very very very very cool"

"very" appears 8 times in the first case, but only 6 in the second case.

"""

def characters_overlapping(left, right):
    """
        Compute the overlap between two sequences.
        Uses a modified version of the Knuth Morris Pratt algorithm
        @param left  Sequence to the left
        @param right Sequence to the right
        left -------------
                      ||||
                      --------------- right
        @return The number of letters that overlap
        (located at the end of "left" and at the beginning of "right")
    """
    if len(left) > len(right):
        k = len(left) - len(right)
        left = left[k:]
    table = compute_back_track_table(right);
    index1 = 0
    index2 = 0
    while (index1 + index2) < len(left):
        if right[index2] == left[index1 + index2]:
            index2 += 1
        else:
            index1 = index1 + index2 - table[index2]
            if index2 > 0:
                index2 = table[index2]
    return index2 # characters matched


def find_or_overlap(left, right):
    """
        Find the sequence "right" in the sequence "left",
        or if not possible, find the overlap.
        Uses a modified version of the Knuth Morris Pratt algorithm
        @param left Sequence to the left
        @param right Sequence to the right

        left -------------
                      ||||
                      --------------- right
        @return The position where the overlap/match starts
            and the number of overlapping/matched characters
    """
    table = compute_back_track_table(right);
    index1 = 0
    index2 = 0
    while (index1 + index2) < len(left):
        if right[index2] == left[index1 + index2]:
            if index2 == len(right) - 1:
                # "right" found in "left"
                return (index1, len(right))
            index2 += 1
        else:
            index1 = index1 + index2 - table[index2]
            if index2 > 0:
                index2 = table[index2]
    # not found.
    index2
    return (len(left) - index2, index2) # characters matched


def compute_back_track_table(seq):
    """ builds the backtrack table of the Knuth Morris Pratt algorithm
        @param seq A string
        @return the backtrack table
    """
    table = [0 for i in range(0,len(seq))]
    table[0] = -1
    table[1] = 0
    position = 2
    cnd = 0
    while position < len(seq):
        if seq[position -1 ] == seq[cnd]:
            table[position] = cnd + 1
            position += 1
            cnd += 1
        elif cnd > 0:
            cnd = table[cnd]
        else:
            table[position] = 0
            position += 1

    return table


class Match:
    """Class to manage the match/overlap between fragments """
    def __init__(self, leftID, rightID, position, n_chars):
        """ Set the initial values
            @param leftID Id of the fragment used as "left"
            @param rightID ID of the fragment used as "right"
            @param position Position in the left fragment where the match starts.
            @param n_chars Number of overlapping characters
        """
        self.leftID = leftID
        self.rightID = rightID
        self.n_chars = n_chars
        self.position = position

    def get_overlaps(self):
        """ return the number of characters overlapping """
        return self.n_chars

    def show(self):
        """ Print values """
        print self.leftID, self.rightID, self.position, self.n_chars

    def __lt__(self, other):
        """ < operator """
        return self.n_chars < other.n_chars

    def __gt__(self, other):
        """ > operator """
        return self.n_chars > other.n_chars



def assemble(fragments):
    """ Assemble a set of fragments
        @param fragments. A set of fragments to assemble
        @return The assembled sentence/s (There can be more than one if full
        assembly is not possible)
    """
    log.info("Assemblying %s sequences", len(fragments))
    n_fragments = len(fragments)
    # sort fragments by increasing length. If there are overlaps between
    # long fragments, this will avoid calculating overlaps between shorter ones
    lenghts = sorted([(len(s),i) for i,s in enumerate(fragments)])
    fragments = [fragments[i] for l,i in lenghts]
    log.debug("Sorted sequences: %s", fragments)

    fragments_dict = {}
    for ID, f in enumerate(fragments):
        fragments_dict[ID] = f
    del fragments

    maxID = n_fragments - 1 # Maximum ID for a fragment so far
    mat = MatchManager()
    fragments_assembled = 0
    while fragments_assembled < n_fragments:
        compute_overlaps(fragments_dict, mat)
        match = mat.get_max_overlapping_pair()
        n_chars = match.get_overlaps()
        if n_chars == 0:
            break
        l = match.leftID
        r = match.rightID
        p = match.position
        maxID += 1
        log.debug("Merging (%s, %s). Overlap: %s", l,r, n_chars)
        log.debug("  * %s",fragments_dict[l])
        log.debug("  * %s",fragments_dict[r])
        # if there is an overlap, concatenate the sentences. Otherwise, do
        # nothing, as the "left" fragment contains "right"
        if (p + n_chars ) == len(fragments_dict[l]):
            fragments_dict[maxID] = fragments_dict[l] + fragments_dict[r][n_chars:]
        else:
            fragments_dict[maxID] = fragments_dict[l]
        log.debug("==> %s", fragments_dict[maxID])
        # remove the fragments just merged
        del fragments_dict[l]
        del fragments_dict[r]
        mat.remove_IDs((r,l))
        fragments_assembled += 1
        log.debug("fragments_assembled so far %s",fragments_assembled)
    return fragments_dict.values()



def compute_overlaps(fragments_dict, mat):
    """ Compute the overlaps between the fragments and store them in a MatchManager
        @param fragments_dict A dictionary of fragments with their ID as key
        @param mat A MatchManager
    """
    log.debug("Computing overlaps")
    # sort the keys in descending order. A greater key will correspond (roughly) to
    # a longer fragment. This condition is only guaranteed during the first iteration
    keys = sorted(fragments_dict.keys())
    keys.reverse()
    for i, j in itertools.product(keys, keys):
        if i == j:
            continue
        if not mat.needs_calculation(i, j, len(fragments_dict[i]), len(fragments_dict[j])):
            continue
        position, n_chars = find_or_overlap(fragments_dict[i], fragments_dict[j])
        log.debug("Overlap between sequences %s (left) and %s (right): %s",i, j, n_chars)
        m = Match(i,j, position, n_chars)
        mat.store(m)


class MatchManager:
    """
        Class to manage the matches between the fragments
    """
    def __init__(self):
        """ Init values """
        self.current_IDs = set()
        self.calculated = dict()
        # match with the maximum absolute overlap
        self.max_match = Match(0,0,0,0)

    def store(self, match):
        """ Store a match between two fragments
            @param match The match between the fragments
        """
        if match > self.max_match:
            self.max_match = match
        self.add_to_calculated(match)

    def add_to_calculated(self, m):
        """ Modifies the internal dictionary to indicate that
            the match (leftID, rightID) has already been calculated.
            @param m a Match
        """
        if m.leftID not in self.calculated:
            self.calculated[m.leftID] = dict()
        self.calculated[m.leftID][m.rightID] = m
        log.debug("Added leftID %s rightID %s",m.leftID,m.rightID)

    def is_calculated(self, leftID, rightID):
        """ Check if the overlap between two
            fragments has been calculated already.
            @param leftID ID of the "left" fragment
            @param rightID ID of the "right" fragment
            @return True if the overlap has been calculated
        """
        if leftID not in self.calculated:
            return False
        if rightID not in self.calculated[leftID]:
            return False
        return True


    def get_max_overlapping_pair(self):
        """ Return the next match with the maximum overlap.
        """
        return self.max_match


    def remove_IDs(self, IDs):
        """ Remove IDs from the dictionary of matches. 
            They are removed both acting as "left" or "right" fragments
            The new max
            @param IDs ids to remove
        """
        log.debug("Removing IDs %s", IDs)
        # remove matches where the IDs is on the left
        for ID in IDs:
            log.debug("Removing %s", ID)
            del self.calculated[ID] # remove "left"
        # remove all matches where the ID appears on the right
        # and find the match with the maximum overlap
        self.max_match = Match(0,0,0,0)
        #import pdb;pdb.set_trace()
        for dictionary in self.calculated.values():
            for ID in IDs:
                if ID in dictionary:
                    del dictionary[ID] # remove "right"
            # calculate max with remaining values
            for m in dictionary.values():            
                self.max_match = max(self.max_match, m)
        m = self.max_match
        log.debug("New best match: leftID %s rightID %s position %s overlap %s",
                m.leftID,m.rightID,m.position,m.n_chars)

    def needs_calculation(self, leftID, rightID, len_left, len_right):
        """ Determines if the overlap between two sequences needs to be done.
            Two conditions are required:

            1) Both fragments need to be longer than the maximum overlap.
            2) The overlap has not been calculated already.

            @param leftID Index of the fragment acting as "left"
            @param rightID Index of the fragment acting as "right"
            @param len_left The length of the fragment acting as "left"
            @param right_left The length of the fragment acting as "right"
            @return True if calculating the overlap is required. False otherwise
        """
        max_overlap = self.max_match.get_overlaps()
        if len_left < max_overlap or len_right < max_overlap:
            return False 
        if self.is_calculated(leftID, rightID):
            return False 
        return True


if __name__ == "__main__":

    import sys
    logging.basicConfig(stream=sys.stdout)
    logging.root.setLevel(logging.ERROR)
    fn = sys.argv[1]

    f = open(fn, "r")
    for line in f:
        fragments = line.rstrip().split(";")
        if len(fragments) == 0:
           continue
        contigs = assemble(fragments)
        for contig in contigs:
            print contig
    f.close()

    exit(0)



