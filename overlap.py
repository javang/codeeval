
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


def assemble(fragments):
    log.info("Assemblying %s sequences", len(fragments))
    # sort fragments by increasing length. If there are overlaps between
    # large sequences it will avoid calculating overlaps between
    # smaller ones
    lenghts = sorted([(len(s),i) for i,s in enumerate(fragments)])
    fragments = [fragments[i] for l,i in lenghts]
    log.debug("Sorted sequences: %s", fragments)

    fragments_dict = {}
    for ID, f in enumerate(fragments):
        fragments_dict[ID] = f
    del fragments # not needed anymore
    n_fragments = len(fragments_dict)
    max_ID = n_fragments -1 # Maximum ID for a fragment
    mat = OverlapMatrix()
    fragments_assembled = 0
    while fragments_assembled < n_fragments - 1:
        compute_overlaps(fragments_dict, mat)
        (n_chars, left_ID, right_ID) = mat.get_max_overlapping_pair()
        if n_chars == 0:
            break
        else:
            # build new fragment with the overlapping pair
            max_ID += 1
            log.debug("Merging (%s, %s): \n  * %s\n  * %s",left_ID, right_ID,
                    fragments_dict[left_ID], fragments_dict[right_ID])
            fragments_dict[max_ID] = fragments_dict[left_ID] + fragments_dict[right_ID][n_chars:]
            log.debug("==> %s", fragments_dict[max_ID] )
            # remove the overlapping pair
            del fragments_dict[left_ID]
            del fragments_dict[right_ID]
            mat.remove_IDs((left_ID, right_ID))
            fragments_assembled += 1
            log.debug("Last overlap had %s characters. " \
            "fragments_assembled so far %s",n_chars,fragments_assembled)

    return fragments_dict.values()


def compute_overlaps(fragments_dict, mat):
    log.debug("Computing overlaps")
    # sort the keys in descending order. Larger key will mean (roughly)
    # larger fragment. This condition is only completely satisfied during
    # the first iteration
    keys = sorted(fragments_dict.keys())
    keys.reverse()
    for i, j in itertools.product(keys, keys):
        if i == j:
            continue
        if not mat.needs_calculation(i, j, len(fragments_dict[i]), len(fragments_dict[j])):
            continue
        n_chars = characters_overlapping(fragments_dict[i], fragments_dict[j])
        log.debug("Overlap between sequences %s (left) and %s (right): %s",i, j, n_chars)
        mat.store(i, j, n_chars)


class OverlapMatrix:
    def __init__(self):
        self.overlaps_heap = []
        self.current_IDs = set()
        self.calculated = dict()

    def store(self, left_ID, right_ID, n_chars):
        """ Store the number of oerlapping characters  between two fragments
            @param left_ID index of the sequence acting as "left"
            @param right_ID index of the sequence acting as "right"
            @param n_chars Number of overlapping characters
        """
        heapq.heappush(self.overlaps_heap, (-1 * n_chars, left_ID, right_ID))
        self.add_to_calculated(left_ID, right_ID)
        self.current_IDs.add(left_ID)
        self.current_IDs.add(right_ID)

    def add_to_calculated(self, i, j):
        if i not in self.calculated:
            self.calculated[i] = dict()
        self.calculated[i][j] = True

    def is_calculated(self, left_ID, right_ID):
        if left_ID not in self.calculated:
            return False
        if right_ID not in self.calculated[left_ID]:
            return False
        return self.calculated[left_ID][right_ID]


    def get_max_overlapping_pair(self):
        """ Return the next pair with the maximum overlap
            Remove pairs involving fragments that are not
            active anymore
        """
        try:
            (n_chars, l, r) = heapq.heappop(self.overlaps_heap)
            while l not in self.current_IDs or \
                 r not in self.current_IDs:
                (n_chars, l, r) = heapq.heappop(self.overlaps_heap)
            return (-1 *n_chars, l, r)
        except IndexError:
            # if the heap fails is because there are no more overlaps
            return (0,-1, -1)


    def remove_IDs(self, IDs):
        """ Remove IDs from the list of active fragments
            @param IDs ids to remove
        """
        log.debug("Removing IDs %s", IDs)
        map(self.current_IDs.remove, IDs)


    def needs_calculation(self, left_ID, right_ID, len_left, len_right):
        """ Determines if the overlap between two sequences needs to be done.

            1) If one sequence is shorter than the maximum overlap, there is
            no need to calculate it now
            2) If the overlap has been calculated already, do not do it again

            @param left_ID Index of the sequence acting as "left"
            @param right_ID Index of the sequence acting as "right"
            @return True if is worth calculating the overlap between the
            sequences. False otherwise
        """
#        log.debug("Check if calculating overlap (%s, %s) is required",
#                                                 left_ID, right_ID)
        if len(self.overlaps_heap) == 0:
            return True # nothing calculated yet
        max_overlap = -1 * self.overlaps_heap[0][0]
        if len_left < max_overlap or len_right < max_overlap:
            return False # sequences too short to find larger maximum overlap
        if self.is_calculated(left_ID, right_ID):
            return False # calculated already
        return True
