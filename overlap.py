
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


def build_contigs(fragments, overlaps):

    # find the indices that are only in the left part of the overlaps.
    # These are the indices where the contigs start
    lefts = set()
    rights = set()
    for (n_chars, l, r) in overlaps:
        rights.add(r)
        lefts.add(l)
        if l in rights:
            lefts.remove(l)
        if r in lefts:
            lefts.remove(r)
    log.debug("The contigs start at sequences %s", lefts)
    # set contigs
    contigs = dict()
    for l in lefts:
        contigs[l] = fragments[l]

    while len(overlaps) > 0:
        ov = overlaps.pop(0)
        n_chars = ov[0] # number of overlapping characters
        l = ov[1] # index of the left sequence in the overlap
        r = ov[2] # index of the right sequence in the overlap
        if l in contigs:
            log.debug("Merging:")
            log.debug("    %s",contigs[l])
            log.debug("    %s",fragments[r])
            contigs[r] = contigs[l] + fragments[r][n_chars:]
            log.debug("obtained:")
            log.debug("    %s", contigs[r])
            del contigs[l]
        else:
            # nothing to do with the overlap yet. Put it back
            overlaps.append(ov)
        overlaps.sort()
        log.debug("Remaining overlaps %s. %s", len(overlaps), overlaps)

    return contigs.values()



def assemble(fragments):
    log.info("Assemblying %s sequences", len(fragments))
    # sort sequences by length. If there are overlaps between
    # large sequences it will avoid calculating overlaps between
    # smaller ones
    lenghts = [(len(s),i) for i,s in enumerate(fragments)]
    lenghts.sort()
    lengths = lenghts.reverse()
    sorted_fragments = [fragments[i] for l,i in lenghts]
    log.debug("Sorted sequences: %s", sorted_fragments)
    fragments_assembled = 0
    overlaps = []
    n_fragments = len(sorted_fragments)
    mat = OverlapMatrix(n_fragments)
    while fragments_assembled < n_fragments - 1:
        compute_overlaps(sorted_fragments, mat)
        (n_chars, left_index, right_index) = mat.get_max_overlapping_pair()
        if n_chars == 0:
            break
        else:
            overlaps.append((n_chars, left_index, right_index))
            fragments_assembled += 1
            log.debug("Last overlap had %s characters. " \
            "fragments_assembled so far %s",n_chars,fragments_assembled)
    log.debug("Assembly %s", overlaps)
    contigs = build_contigs(sorted_fragments,overlaps)
    return contigs



def compute_overlaps(fragments, mat):
    log.debug("Computing overlaps")
    n_fragments = len(fragments)
    for i, j in itertools.product(xrange(n_fragments), xrange(n_fragments)):
        if i == j:
            continue
        if not mat.needs_calculation(i, j, len(fragments[i]), len(fragments[j])):
            continue
        n_chars = characters_overlapping(fragments[i], fragments[j])
        #log.debug("Overlap between sequences %s (left) and %s (right): %s",i, j, n_chars)
        mat.store(i, j, n_chars)




class OverlapMatrix:
    def __init__(self, n_sequences):
        self.overlaps_heap = []
        self.start_dictionaries(n_sequences)

    def start_dictionaries(self, n_fragments):
        """ Build the dictionary used to check if the overlap
        between 2 sequences has been calculated already

            @param n_fragments Number of sequences considered
        """
        self.used_as_left = dict()
        self.used_as_right = dict()
        self.calculated = dict()
        for l in xrange(n_fragments):
            self.used_as_left[l] = False
            self.used_as_right[l] = False
            self.calculated[l] = dict()
            for r in xrange(n_fragments):
                self.calculated[l][r] = False

    def store(self, left_index, right_index, n_chars):
        """ Store the number of oerlapping characters  between two fragments
            @param left_index index of the sequence acting as "left"
            @param right_index index of the sequence acting as "right"
            @param n_chars Number of overlapping characters
        """
        heapq.heappush(self.overlaps_heap, (-1 * n_chars, left_index, right_index))
        self.calculated[left_index][right_index] = True

    def get_max_overlapping_pair(self):
        """ Return the next pair with the maximum overlap
            Remove pairs that were calculated but need to be descarded
            because they describe impossible overlaps (given the
            previous overlaps)
        """
        try:
            (n_chars, l, r) = heapq.heappop(self.overlaps_heap)
            while self.used_as_left[l] == True and \
                  self.used_as_right[r] == True:
                (n_chars, l, r) = heapq.heappop(self.overlaps_heap)
            self.mark_as_used(l,r)
            return (-1 *n_chars, l, r)
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
#        log.debug("Check if calculating overlap (%s, %s) is required",
#                                                 index_left, index_right)
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
