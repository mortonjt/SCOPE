#!/usr/bin/python
import sys
from dnaSeq import *
from re import *

D = []

def count_trimPoly(TRIM_file):
    total = 0
    num_found = 0
    mlen = 1000000000

    S = {}   # Set of ids of identified sequences
    for line in open(TRIM_file):
        total += 1
        arr = split("\s+", line.strip())
        info = arr[0]
        start, stop, plen = [int(x) for x in arr[2:]]
        if start > 1 or stop < plen:
            S[info] = [start-1, stop, plen]

    return S

extensions = ["poly_A", "polyA_", "polyAA", "polyAT", "poly_T", "polyT_", "polyTA", "polyTT", "reject"]
def count_scope(scope_file_base):
    S = {}
    for ext in extensions:
        with open(scope_file_base + ".d/" + scope_file_base + ".d." + ext) as fp:
            for line in fp:
                if line[0] == ">":
                    r = search("length=(\d+)\(type=\w,start=(\d+),end=(\d+)\)", line)                        
                    if r:
                        length,start,end = int(r.group(1)), int(r.group(2)), int(r.group(3))
                        info = split("\s+", line)[0][1:]
                        S[info] = [start, end, length]


    return S

def scopeAnalysis(R, s):
    """Add tail length and ptail length to R"""
    print(R)
    R.extend([R[0], R[1]-R[0], R[2]-R[1]])
    print(R)

def trimAnalysis(R, s):
    flip = R[0] == 0
    if flip:
        s = "".join(s[::-1])
        R[0], R[1] = R[2]-R[1], R[2]-R[0]

    right_end = R[0]
    left_end = R[0]-10
    window = s[left_end:right_end]
    c = max([(window.count(x),x) for x in "ACGT"])[1]    # Find most common character in the window
    if not (c in "AT") or window.count(c) < 0.5*le(window):
        R.extend([None]*3)
    
    else:

        # Shift window to the left until we find one with less than 50% characters
        while left_end > 0 and window.count(c) >= 0.5*len(window):
            right_end -= 1
            left_end -= 1
            window = s[left_end:right_end]

        if window.count(c) >= 0.5*len(window):
            right_end += 1
            left_end += 1

        while s[left_end] != 'c':
            left_end += 1

        if right_end - left_end < 10:
            R.extend([None]*3)
        else:
            R.extend(left_end, R[0] - left_end, R[2] - R[0])

    if flip:
        R[0], R[1] = R[2]-R[1], R[2]-R[0]
        R[3], R[5] = R[5], R[3]

        
def analysis(file_base, TRIM, SCOPE):
    count = 0
    for s in iterateFA(file_base + ".fasta"):
        if SCOPE.get(s.id):
            scopeAnalysis(SCOPE[s.id], s)
        if TRIM.get(s.id):
            trimAnalysis(TRIM[s.id], s)
        if count == 9:
            break


if __name__ == "__main__":
    file_base = sys.argv[1]
    
    TRIM = count_trimPoly(file_base + ".trimpoly.out2")
    SCOPE = count_scope(file_base)

    analysis(file_base, TRIM, SCOPE)

    TRIM_SET = set(TRIM.keys())
    SCOPE_SET = set(SCOPE.keys())

    for i in range(10):
        id = TRIM_SET.pop()
        print("TRIM  %s %s" % (id, TRIM[id]))
        id = SCOPE_SET.pop()
        print("SCOPE %s %s" % (id, SCOPE[id]))
    
    #X = TRIM_SET - SCOPE_SET
    
    print("SCOPE ONLY: %d, TRIM ONLY: %d, BOTH: %d" % (len(SCOPE_SET - TRIM_SET), len(TRIM_SET - SCOPE_SET), len(SCOPE_SET & TRIM_SET)))
