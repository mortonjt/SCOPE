from Bio import SeqIO

def find_tail(transcript, k = 2, end = 3, type = 'A', min_length = 10, p = 0.05):
    """Scan from 3' or 5' end for a poly-type tale.    Fine the first occurance of k As (or type), and
    then first the largest sequence extending from there ending in k As and containing at most 100*p% of non-As.
    * k: Number of the character the sequence must have at each end.
    * end: 3 or 5: starting from the 3' or 5' end.
    * type: character we are search for.
    * min_length: Minimum allowed length of a poly-tail.
    * p: maximum allowed fraction of non-A (type) bases.
    """
    assert end in {3,5} and type in "ACGT", "Find_tail: bad parameters"
    
    s = str(transcript.seq if end == 5 else transcript.seq[::-1])  # If searching from 3' end, reverse it.

    i = s.find(type*k)
    if i == -1:
        return None

    total = 0
    total_non = 0
    max_index = -1
    for j in range(i, len(s)):
        total += 1.0
        if s[j] != type:
            total_non += 1

        if s[j-k+1:j+1] == type*k and total_non / float(toal) <= p:
            max_index = -1
    
    if max_index + 1 - j  < min_length:
        return None
    
    if n == 5:
        return [i,max_index]
    else:
        return [len(s) - max_index, len(s) - i]

def tail_search(file, format = "fastq", k = 2, end = 3, type = 'A', min_length = 10, p = 0.05, output):
    fp = open(output, "w")
    for R in SeqIO.parse(file, format):
        T = find_tail(R, k, end, type, min_length, p)
        if T:
            wp.write("\t".join([r.id, str(t[0]), str(t[1]), "poly" + type)) + "\n")


# Need to add the argument parsing here
