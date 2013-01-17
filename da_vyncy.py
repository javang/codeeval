#!/usr/bin/python

import sys
import os
import time
import logging
import itertools
log = logging.getLogger("assembler")


def characters_overlapping(left, right):
    """
        Compute the overlap between two strings.
        Uses a modified version of the Knuth Morris Pratt algorithm
        @param left Sequence to the left
        @param right Sequence to the right
        left -------------
                      ||||
                      --------------- right
        @return The number of letters that overlap
        (located at the end of left and at the beginning
        of right)
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
    """Class to manage the match/overlap between two text fragments """
    def __init__(self, leftID, rightID, position, n_chars):
        """ Set the initial values
            @param leftID Id of the fragment used as "left"
            @param rightID ID of the fragment used as "right"
            @param position Position in the left fragment where  the Match
                    starts.
            @param n_chars Number of overlapping characters
        """
        self.leftID = leftID
        self.rightID = rightID
        self.n_chars = n_chars
        self.position = position

    def get_overlaps(self):
        """ return the number of characters overlapping.
        """
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
        @param fragmnents. A set of strings to Assemble
        @return The assemble sentence/s (There can be more than one
        sentences if the entire sentence could not be assembled).
    """
    log.info("Assemblying %s sequences", len(fragments))
    # sort fragments by increasing length. If there are overlaps between
    # large sequences it will avoid calculating overlaps between
    # smaller ones
    lenghts = sorted([(len(s),i) for i,s in enumerate(fragments)])
    fragments = [fragments[i] for l,i in lenghts]
    log.debug("Sorted sequences: %s", fragments)
    # put the fragments in a dictionary
    fragments_dict = {}
    for ID, f in enumerate(fragments):
        fragments_dict[ID] = f
    del fragments
    n_fragments = len(fragments_dict)
    maxID = n_fragments - 1 # Maximum ID for a fragment so far
    mat = MatchManager()
    fragments_assembled = 0
    while fragments_assembled < n_fragments:
        compute_overlaps(fragments_dict, mat)
        match = mat.get_max_overlapping_pair()
        n_chars = match.get_overlaps()
        if n_chars == 0:
            break
        leftID = match.leftID
        rightID = match.rightID
        position = match.position
        maxID += 1
        log.debug("Merging (%s, %s). Overlap: %s", leftID, rightID, n_chars)
        log.debug("  * %s",fragments_dict[leftID])
        log.debug("  * %s",fragments_dict[rightID])
        # if there is an overlap, concatenate the sentences. Otherwise, do
        # nothing, as the "left" fragment contains "right"
        if (position + n_chars ) == len(fragments_dict[leftID]):
            fragments_dict[maxID] = fragments_dict[leftID] + \
                                     fragments_dict[rightID][n_chars:]
        else:
            fragments_dict[maxID] = fragments_dict[leftID]
        log.debug("==> %s", fragments_dict[maxID])
        # remove the fragments just merged
        del fragments_dict[leftID]
        del fragments_dict[rightID]
        mat.remove_IDs((leftID, rightID))
        fragments_assembled += 1
        log.debug("fragments_assembled so far %s",fragments_assembled)
    return fragments_dict.values()



def compute_overlaps(fragments_dict, mat):
    """ Compute the overlaps between the fragments and store them in a MatchManager
        @param fragments_dict a dictionary of fragments with their ID as key
        @param mat A MatchManager
    """
    log.debug("Computing overlaps")
    # sort the keys in descending order. Larger key will mean (roughly)
    # larger fragment. This condition is only entirely satisfied during
    # the first iteration
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
            the match (leftID, rightID) has already been calculated
            @param m a Match
        """
        if m.leftID not in self.calculated:
            self.calculated[m.leftID] = dict()
        self.calculated[m.leftID][m.rightID] = m
        log.debug("Added leftID %s rightID %s",m.leftID,m.rightID)
        self.show_calculated_dict()

    def is_calculated(self, leftID, rightID):
        """ Check if the overlap between two
            fragments has been calculated already
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
        """ Remove IDs from the list of active fragments
            @param IDs ids to remove
        """
        log.debug("Removing IDs %s", IDs)
        # remove matches where the IDs is on the left
        for ID in IDs:
            log.debug("Removing %s", ID)
            del self.calculated[ID]
        # remove all matches where the ID appears on the right
        # and find the match with the maximum overlap
        self.max_match = Match(0,0,0,0)
        for ID in IDs:
            for k in self.calculated:
                if ID in self.calculated[k]:
                    del self.calculated[k][ID] # remove "left"
                for s in self.calculated[k]: # remove "rights"
                    self.max_match = max(self.calculated[k][s], self.max_match)
        self.show_calculated_dict()

    def show_calculated_dict(self):
        print "========================"
        print self.calculated
        for k in self.calculated:
            log.debug("DICTIONARY FOR %s", k)
            dic = self.calculated[k]
            for s in dic:
                log.debug("    leftID %s rightID %s position %s overlap %s",
                dic[s].leftID,dic[s].rightID,dic[s].position,dic[s].n_chars)
        print "========================"





    def needs_calculation(self, leftID, rightID, len_left, len_right):
        """ Determines if the overlap between two sequences needs to be done.

            1) If one sequence is shorter than the maximum overlap, there is
            no need to calculate it now
            2) If the overlap has been calculated already, do not do it again

            @param leftID Index of the sequence acting as "left"
            @param rightID Index of the sequence acting as "right"
            @param len_left The length of the fragment acting as "left"
            @param right_left The length of the fragment acting as "right"
            @return True if is worth calculating the overlap between the
            sequences. False otherwise
        """
        max_overlap = self.max_match.get_overlaps()
        if len_left < max_overlap or len_right < max_overlap:
            return False # sequences too short to find a greater maximum overlap
        if self.is_calculated(leftID, rightID):
            return False # calculated already
        return True


if __name__ == "__main__":

    import sys
    logging.basicConfig(stream=sys.stdout)
    logging.root.setLevel(logging.DEBUG)
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



"""
A typical example of the case where multiple reconstructions are possible is the presence of repeats. For example, let's say that we have 3 fragments:

1) "bioinformatics is very very"
2) "very very very very"
3) "very very not very very very very cool"

Fragments 1 and 2 have an overlap of 9 characters "very very". The same is true
for fragments 1 and 3. If we choose to assemble in the order: 1-2-3 we would get:

1-2:  "bioinformatics is very very very very"
1-2-3:  :  "bioinformatics is very very very very not very very very very cool"

On the other hand, if we choose to assemble as 1-3-2:
1-3: "bioinformatics is very very not very very very very cool"
1-3-2 "bioinformatics is very very not very very very very cool"

"very" appears 8 times in the first case, but only 6 in the second case

"""