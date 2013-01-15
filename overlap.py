


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
                index2 = table(index2)
    return index2 # characters matched




def compute_back_track_table(seq):
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




def assemble(seqs):

    # sort sequences by length
    lenghts = [(len(s),i) i,s for enumerate(seqs)]
    lenghts.sort()
    lengths = lenghts.reverse()
    sorted_seqs = [seqs(i) for l,i in lenghts]

    # keep calculating the overlaps until there are no more
    # overlaps or all fragments are assembled
    max_overlap = -1
    fragments_assembled = 0
    assembly = []
    n_seqs = len(seqs)
    mat = OverlapMatrix(n_seqs)
    while max_overlap != 0 or fragments_assembled < len(n_seqs) -1:
        compute_overlaps(seqs, mat)
        pair, max_overlap = mat.get_max_overlapping_pair()
        if pair:
            assembly.append(pair)
            fragments_assembled += 1


def compute_overlaps(seqs, mat)
    for i in range(n_seqs):
        if not mat.needs_calculation(i, len(seqs[i])):
            continue
        for j in range(i+1, n_seqs):
            if not mat.needs_calculation(j, len(seqs[j])):
                continue
            if not mat.is_calculated(i,j):
                mat.store(i,j, characters_overlapping(seqs[i], seqs[j]))
            if not mat.is_calculated(j, i):
                mat.store(j,i, characters_overlapping(seqs[j], seqs[i]))




class OverlapMatrix:

    def __init__(self, n_sequences):
        self.max_overlap = -1
        self.max_pair = (-1,-1)
        self.overlaps_dict = dict()
        for i in range(n_seqs):
            self.overlaps_dict[i] = dict()

    def store(self, i, j, n_overlap):
        """
            i index of the left Sequence
            j index of the right sequence
        """
        if i not in self.overlaps_dict:
            self.overlaps_dict[i] = dict()
        self.overlaps_dict[i][j] = n_overlap

        if n_overlap > self.max_overlap:
            self.max_overlap = n_overlap
            self.max_pair = (i,j)

    def is_calculated(self, i, j):
        if i not in self.overlaps_dict:
            return False
        if j not in self.overlaps_dict[i]:
            return False
        return True


    def get_max_overlapping_pair(self):
        if self.max_overlap <= 0:
            return False, self.max_overlap


    def needs_calculation(self, i, length_seq):
        if not i in self.overlaps_dict:
            return False
        if length_seq < self.max_overlap:
            return False
        else:
            return True
