#from Bio import SeqIO
#from Bio.SeqRecord import SeqRecord
import argparse

def find_tail(transcript, x = 2, end = 3, base_type = 'A', min_length = 10, p = 0.05):
    """Scan from 3' or 5' end for a poly-type tale.    Find the first occurance of x As (or type), and
    then first the largest sequence extending from there ending in x As and containing at most 100*p% of non-As.
    * x: Number of the character the sequence must have at each end.
    * end: 3 or 5: starting from the 3' or 5' end.
    * base_type: character we are search for.
    * min_length: Minimum allowed length of a poly-tail.
    * p: maximum allowed fraction of non-A (base_type) bases.
    """
    assert end in {3,5} and base_type in "ACGT", "Find_tail: bad parameters"
    
    s = str(transcript.seq if end == 5 else transcript.seq[::-1])  # If searching from 3' end, reverse it.


    i = s.find(base_type*x)
    if i == -1:
        return None

    total_non = 0
    max_index = -1
    for j in range(i, len(s)):
        if s[j] != base_type:
            total_non += 1
        if (s[j-x+1:j+1] == base_type*x) and (total_non / float(j-i+1) <= p):
            max_index = j
    
    if max_index < 1 or max_index + 1 - i  < min_length:
        return None
    
    if end == 5:
        return [i,max_index]
    else:
        return [len(s) - max_index, len(s) - i]

def find_tail2(transcript, x = 2, end = 3, base_type = 'A', min_length = 10, p = 0.05):
    """Scan from 3' or 5' end for a poly-type tale.  Find the longest sequence that starts and ends with x As has no more that 100*p% non-As."""
    assert end in {3,5} and base_type in "ACGT", "Find_tail: bad parameters"
    
    s = str(transcript.seq if end == 5 else transcript.seq[::-1])  # If searching from 3' end, reverse it.

    count = [0]*(len(s)+1)        # count[i] will be the number of non-A characters in s[0:i]
    start_positions = set()       # coordinates that are the start of a maximal sequence of >=  As.  
    end_positions = set()         # coordinates that are the end of a maximal sequence of >= x As
    continuous = 0

    for i in range(len(s)):
        if s[i] == base_type:
            continuous += 1
            if continuous == x:
                start_positions.add(i-x+1)
            count[i+1] = count[i]
        else:
            if continuous >= x:
                end_positions.add(i)
            continuous = 0
            count[i+1] = count[i] + 1


    best = (0, -1,-1)
    for i in sorted(list(start_positions)):
        for j in sorted(list(end_positions)):
            if i < j and (count[j] - count[i]) / float(j - i) <= p:
                best = max(best, (j-i, i, j))

    print("HERE1")
    if best[0] < min_length:
        return None

    print("HERE2")
    if end == 5:
        print("HERE3")
        return [best[1], best[2]]
    else:
        print("HERE4")
        return [len(s) - best[2], len(s) - best[1]]


def tail_search(input, output, file_format = "fastq", x = 2, end = 3, base_type = 'A', min_length = 10, p = 0.05, search_alg = 'find_tail'):
    fp = open(output, "w")
    for R in SeqIO.parse(input, file_format):
        T = search_alg(R, x, end, base_type, min_length, p)
        if T:
            fp.write("\t".join([R.id, str(T[0]), str(T[1]), "poly" + base_type]) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Simple scan of sequences for tails")

    alg_parser = parser.add_argument_group("Choice of algorithm")
    alg_group = alg_parser.add_mutually_exclusive_group()
    alg_group.add_argument("-b", "--basic_alg", dest = 'alg', action = 'store_const', const = 'find_tail', help = "Basic linear-scan algorithm", default = 'find_tail')
    alg_group.add_argument("-s", "--second_alg", dest = 'alg', action = 'store_const', const = 'find_tail2', help = "More advanced (basic) scan algoritn")
    
    parameters_group = parser.add_argument_group("Algorithmic parameters")
    parameters_group.add_argument('-x', dest = 'x', type = int, help = "Number of bases at each end of tail", default = 2);
    parameters_group.add_argument('-e', dest  = 'end', type = str, help = "End to search from (3 to search from right end, 5 from left)", default = 3)
    parameters_group.add_argument('-t', '--base_type', dest = 'base_type', type = str, help = "Character making up the poly-tail", default = 'A')
    parameters_group.add_argument('-m', '--min_length', dest = 'min', type = int, help = "Minimum tail length", default = 10)
    parameters_group.add_argument('-p', dest = 'p', type = float, help = 'Maximum fraction of impurities', default = 0.05)

    parser.add_argument('-f', '--format', dest = 'format', type = str, help = 'input file format', default = 'fastq')

    parser.add_argument('input', help = 'input file')
    parser.add_argument('output', help = 'output file')
    

    args = parser.parse_args()
    
    eval(args.alg)(args.input, args.output, args.format, args.x, args.end, args.base_type, args.min, args.p)
