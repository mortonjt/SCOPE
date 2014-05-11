# NOTE: 'N' represents the ambiguous base -- one whose content has not been 
# determined.  For purposes of compelmentation, the complement of 'N' is 'N'.
BASES = {'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n'}
COMPLEMENT = {x:y for x,y in zip(['A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n'],
                                  ['T', 'G', 'C', 'A', 'N', 't', 'g', 'c', 'a', 'N'])}

class DNAError(Exception):
    def __init__(self, value = "Bad DNA string"):
        self.value = value

    def __str__(self):
        return str(self.value)

class DNAiterator:
    def __init__(self, dnaSeqObj):
        self.it = iter(dnaSeqObj.seq)

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.it)


class dnaSeq:
    """Mutable DNA Sequence class."""

    # 1) Implement the constructor.  Should take one argument (other than self) and throw a ValueError if
    #    the value is not a string, or a DNAError see above) if the value is a string but contains bad 
    #    bases characters.  (That is, characters that are not in the BASES set.)  Otherwise, store the string
    #    is an attribute in whatever for you see fit.
    #    COMMENT: Since this is intended to be a mutable class, yo might want to consider holding
    #    the information in a mutable data strcture -- as opposed to one that you have to copy
    #    every time you want to change the object.
    def __init__(self, s):
        if type(s) != str:
            raise ValueError("Bad value: " + str(s))
        if not set(s) <= BASES:
            raise DNAError
        
        self.seq = list(s);    # I'm storing as an array, as a string isn't mutable.  Will make like easier in future functions.
        
    #---------------------------------
    # 2) Implement the __len__ function so that len(o) will return the length of sequence object o.
    def __len__(self):
        return len(self.seq)
    
    #---------------------------------
    # 3) Implement the __str__ function so that str(o) will return the the sequence as a string object.
    def __str__(self):
        return "".join(self.seq)
    
    """The following was not required, but is useful when working in the Python shell."""
    def __repr__(self):
        return "dnaSeq: " + str(self)
    
    #---------------------------------
    # 5) Implement the __getitem__ function so that we can us the index the operator on a dnaSeq -- 
    #    both single indexing and splicing.  An index should return a single-character string, 
    #    while a slice should return a dnaSeq object.
    #    NOTE: The Python specifications for __getitem__ require that:
    #          1) An illegal index type result in a TypeError being raised.
    #          2) An invalid index value results in a KeyError being raised.
    #    You may assume this for any other object, and implement this for dnaSeq.
    #    There are two legal index types: int and slice.  (The later being the class of the
    #    object __getitem__ gets for a call such as o[5:7].
    def __getitem__(self, i):
        if type(i) == int:
            return self.seq[i]
        if type(i) == slice:
            return dnaSeq("".join(self.seq[i]))
        raise TypeError("Bad index")

    #---------------------------------
    # 6) Implement the __setitem__ function so we can use index-assignment (e.g. o[5] = 'A' or O[5:7] = 'AC').
    #    Should throw the same sorts of errors as __getitem__.  Additionally, it should throw a ValueError
    #    if the value is not either a dnaSeq object or a string.
    #    (If value is a string, throw a DNAErro if it contains any non-base characters.)
    def __setitem__(self, i, value):
        if type(value) == str:
            if set(value) <= BASES:
                self.seq[i] = value
            else:
                raise DNAError
        elif type(value) == dnaSeq:
            self.seq[i] = dnaSeq.seq
        else:
            raise ValueError
 
    #---------------------------------
    # 7) Impelement the __add__ function to enable + to concatenate sequences.
    #    Throw a ValueError if object other is not of class dnaSeq
    #    Make sure it produces a deep copy.
    def __add__(self, other):
        if type(other) != dnaSeq:
            raise ValueError
        return dnaSeq(str(self) + str(other))

    #---------------------------------
    # 8) Implement the necessary method for the "in" operator to to allow operations such as:
    #     s in o (where s is a string or dnaSeq, o is a dnaSeq, and this returns true if either
    #     o contains s as a subsequence.
    #     Raise a ValueError if the query sequence is neither a string nor a dnaSeq.
    def __contains__(self, query):
        if type(query) in {dnaSeq, str}:
            return str(query) in str(self)
        raise ValueError


    #---------------------------------
    # 9) Make the class iterable -- allowing it to be used by for statements and comprehenstions.
    #     e.g. Once this is implemented, this should work: [print(i) for i in o], where o is a dnaSeq
    #     object.  
    #     See project specifications for instructions on how to do this. 
    def __iter__(self):
        return iter(self.seq)         # Just use the list iterator
        #return DNAiterator(self)     # Alternative: return our own iterator class

    def __iadd__(self, o):
        
        if type(o) in {str,dnaSeq}:
            return dnaSeq(str(self) + str(o))
        else:
            raise ValueError
    
    #---------------------------------------------------------
    # Object methods: Mutators
    
    #---------------------------------
    # 11) In-place complement the sequence.  (No return value.)
    def complement(self):
        self.seq = [COMPLEMENT[x] for x in self.seq]
                    
    #---------------------------------
    # 12) method reverse(self): Reverse the order of the sequence in-place.
    def reverse(self):
        self.seq = self.seq[::-1]

    #---------------------------------
    # 13) method reverse_complement(self): Reverse-complements the sequence in-place.  
    def reverse_complement(self):
        self.seq = [COMPLEMENT[x] for x in self.seq[::-1]]
        

    def count(self, c):
        return str(self).count(c)

def readFA(fasta_file):
    """Read sequences in from an .fa file, returning a list of sequence, each
    with an added info attribute containing the corresponding info line."""
    
    S = []

    with open(fasta_file) as fp:
        info = fp.readline().rstrip()
        seq = ""
        while True:
            line = fp.readline()
            if not line or line[0] == ">":
                s = dnaSeq(seq)
                s.info = info
                S.append(s)
                if not line:
                    break
                info = line.rstrip()
                seq = ""
            else:
                seq += line.rstrip()
    return S

def iterateFA(fasta_file):
    """Read sequences in from an .fa file, returning a list of sequence, each
    with an added info attribute containing the corresponding info line."""
    
    S = []

    with open(fasta_file) as fp:
        info = fp.readline().rstrip()
        seq = ""
        while True:
            line = fp.readline()
            if not line or line[0] == ">":
                s = dnaSeq(seq)
                s.info = info
                s.id = s.info[1:s.info.find(" ")]
                yield s
                if not line:
                    break
                info = line.rstrip()
                seq = ""
            else:
                seq += line.rstrip()

def writeFA(S, output_file, col_width = 60):
    """Write a list of fast sequence (with info attributes) to a .fa file"""
    with open(output_file, "w") as fp:
        for s in S:
            fp.write(s.info + "\n" + "\n".join([str(s[i:i+col_width]) for i in range(0, len(s), col_width)]) + "\n")    
