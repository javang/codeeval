
import heapq
import itertools
import logging
log = logging.getLogger("overlap")

def characters_overlapping(left, right):
    """
        Compute the overlap between two strings
        Using the Knuth Morris Pratt algorithm
        @param left Sequence to the left
        @param right Sequence to the right

        left -------------
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




def compute_back_track_table(seq):
    """ builds the backtrack table of the Knuth Morris Pratt algorithm
        @param seq A strings
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


def build_contigs(seqs, overlaps):

    # find the indices that are only in the left part of the overlaps.
    # These are the indices where the contigs start
    lefts = set()
    rights = set()
    for (n_chars, i, j) in overlaps:
        rights.add(j)
        lefts.add(i)
        if i in rights:
            lefts.remove(i)
        if j in lefts:
            lefts.remove(j)
    contigs = dict()
    for i in lefts:
        contigs[i] = seqs[i]

    while len(overlaps) > 0:
        ov = overlaps.pop(0)
        n_chars = ov[0] # number of overlapping characters
        l = ov[1] # index of the left sequence in the overlap
        r = ov[2] # index of the right sequence in the overlap
        if l in contigs:
            contigs[r] = contigs[l] + seqs[r][n_chars:]
            del contigs[l]
        else:
            # nothing to do with the overlap yet. Put it back
            overlaps.append(ov)
    return contigs.values()



def assemble(seqs):
    log.info("Assemblying %s sequences", len(seqs))
    # sort sequences by length. If there are overlaps between
    # large sequences it will avoid calculating overlaps between
    # smaller ones
    lenghts = [(len(s),i) for i,s in enumerate(seqs)]
    lenghts.sort()
    lengths = lenghts.reverse()
    sorted_seqs = [seqs[i] for l,i in lenghts]
    log.debug("Sorted sequences: %s", sorted_seqs)
    fragments_assembled = 0
    overlaps = []
    n_seqs = len(sorted_seqs)
    mat = OverlapMatrix(n_seqs)
    while fragments_assembled < n_seqs - 1:
        compute_overlaps(sorted_seqs, mat)
        (n_chars, left_index, right_index) = mat.get_max_overlapping_pair()
        if n_chars == 0:
            break
        else:
            overlaps.append((n_chars, left_index, right_index))
            fragments_assembled += 1
            log.debug("fragments assembled %s",fragments_assembled)
    log.info("Assembly %s", overlaps)
    contigs = build_contigs(sorted_seqs,overlaps)
    return contigs



def compute_overlaps(seqs, mat):
    log.info("Computing overlaps")
    n_seqs = len(seqs)
    for i, j in itertools.product(xrange(n_seqs), xrange(n_seqs)):
        if i == j:
            continue
        if not mat.needs_calculation(i, j, len(seqs[i]), len(seqs[j])):
            continue
        n_chars = characters_overlapping(seqs[i], seqs[j])
        log.debug("Overlap between sequences %s (left) and %s (right): %s",i, j, n_chars)
        mat.store(i, j, n_chars)


class OverlapMatrix:

    def __init__(self, n_sequences):
        self.overlaps_heap = []
        self.start_dictionaries(n_sequences)

    def start_dictionaries(self, n_seqs):
        """ Build the dictionary used to check if the overlap
        between 2 sequences has been calculated already

            @param n_seqs Number of sequences considered
        """
        self.used_as_left = dict()
        self.used_as_right = dict()
        self.calculated = dict()
        for i in xrange(n_seqs):
            self.used_as_left[i] = False
            self.used_as_right[i] = False
            self.calculated[i] = dict()
            for j in xrange(n_seqs):
                self.calculated[i][j] = False

    def store(self, i, j, n_chars):
        """ Store the number of oerlapping characters between sequences i and j
            @param i index of the sequence acting as "left"
            @param j index of the sequence acting as "right"
            @param n_chars Number of overlapping characters
        """
        heapq.heappush(self.overlaps_heap, (-1 * n_chars, i, j))
        self.calculated[i][j] = True

    def get_max_overlapping_pair(self):
        """ Return the next pair with the maximum overlap
            Remove pairs that were calculated but need to be descarded
            because they describe impossible overlaps (given the
            previous overlaps)
        """
        try:
            (n_chars, i, j) = heapq.heappop(self.overlaps_heap)
            while self.used_as_left[i] == True and \
                  self.used_as_right[j] == True:
                (n_chars, i, j) = heapq.heappop(self.overlaps_heap)
            self.mark_as_used(i,j)
            return (-1 *n_chars, i, j)
        except IndexError:
            # if the heap fails is because there are no more overlaps
            return (0,-1, -1)


    def mark_as_used(self, left_index, right_index):
        """ Mark that two sequences have been used already in one overlap
            @param left_index Index of the sequence acting as left
            @param rigth_index Index of the sequence acting as right
        """
        log.debug("Cleaning dictionary")
        self.used_as_left[left_index] = True
        self.used_as_right[right_index] = True


    def needs_calculation(self, index_left, index_right, len_left, len_right):
        """ Determines if the overlap between two sequences needs to be done.

            1) If one sequence is shorter than the maximum overlap, there is
            no need to calculate it now
            2) If the overlap has been calculated already, do not do it again

            @param index_let Index of the sequence acting as "left"
            @param index_right Index of the sequence acting as "right"
            @return True if is worth calculating the overlap between the
            sequences i and j. False otherwise
        """
        log.debug("Check if calculating overlap (%s, %s) is required",
                                                 index_left, index_right)
        if len(self.overlaps_heap) == 0:
            return True # nothing calculated yet
        max_overlap = -1 * self.overlaps_heap[0][0]
        if len_left < max_overlap or len_right < max_overlap:
            return False # sequences too short to find larger maximum overlap
        if self.used_as_left[index_left] or self.used_as_right[index_right]:
            return False # the sequences have been assembled already as left or right
        if self.calculated[index_left][index_right]:
            return False # calculated already
        return True
