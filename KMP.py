
import logging
log = logging.getLogger("overlap")

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

