
import KMP
import itertools
import logging
log = logging.getLogger("overlap")


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
        position, n_chars = KMP.find_or_overlap(fragments_dict[i], fragments_dict[j])
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


