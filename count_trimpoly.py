import sys
from re import *

if __name__ == "__main__":
    SCOPE_base, TRIM_file = sys.argv[1:3]
    min_len = int(sys.argv[3])
    

    total = 0
    num_found = 0
    mlen = 1000000000
    for line in open(TRIM_file):
        total += 1
        arr = [int(x) for x in split("\s+", line.rstrip())[2:5]]
        
        flen = arr[0] if arr[0] > 1 else (arr[2] - arr[1])
        if flen >= min_len:
            num_found += 1
            mlen = min(flen, mlen)

    print("TRIM: %d/%d  %5.3f   %d" % (num_found, total, float(num_found)/total, mlen))
