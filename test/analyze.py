"""
Analyzes results from Tail-Seq and compares it to SCOPE
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cbook as cbook
from matplotlib._png import read_png
from matplotlib.offsetbox import OffsetImage 
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg

import numpy
import pylab

from collections import defaultdict
from collections import Counter
from pandas import *

def parseCSV(infile):
    lengths = []
    with open(infile,'r') as handle:
        for ln in handle:
            ln = ln.rstrip()
            toks = ln.split(',')
            if toks[0]=="accession":continue
            pAlength = int(toks[5])
            lengths.append(pAlength)
    return lengths
def parseTab(infile):
    lengths = []
    with open(infile,'r') as handle:
        for ln in handle:
            ln = ln.rstrip()
            toks = ln.split('\t')
            st,end = int(toks[1]),int(toks[2])
            pAlength = end-st
            lengths.append(pAlength)
    return lengths    
def histogram(lengths):
    counter = Counter()
    for L in lengths:
        counter[L]+=1
    return counter.keys(),counter.values()
    
if __name__=="__main__":
    workspace="/media/HD/Documents/Jamie/MiamiBio/SCOPE/workspace"
    scopefile   = "%s/tailseq_out.fa.tab"%(workspace)
    tailseqfile = "%s/GSM1242325_HeLa-polya.csv"%(workspace)
    scopeLens   = parseTab(scopefile)
    tailseqLens = parseCSV(tailseqfile)
    scopeAxis,scopeDensity = histogram(scopeLens)
    tailseqAxis,tailseqDensity = histogram(tailseqLens)
    plt.plot(scopeAxis,scopeDensity,'-r',tailseqAxis,tailseqDensity,'-b')
    plt.show()
    
