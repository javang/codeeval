


def max_overlap(seq1, seq2):
    """
        Compute the overlap between two strings
        Using the Knuth Morris Pratt algorithm
        @param    
    """
        
    if len(seq1) > len(seq2):
        seq1 = seq1[len(seq2):-1]
    table = compute_back_track_table(seq2);
    index1 = 0
    index2 = 0
    while (index1 + index2) < len(seq1):    
        if seq2[index2] == seq1[index1 + index2]:
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


