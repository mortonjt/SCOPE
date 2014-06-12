"""
This parses the output from SCOPE++
"""
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord

class SCOPEparser(object):
    def __init__(self,fasta,tab):
        self.fasta = fasta
        self.tab = tab
        self.parseTab()
        pass
    def parseFASTA(self):
        pass
    def parseFASTQ(self):
        pass
    def parseTab(self):
        self.titles,self.starts,self.ends,self.ptypes = [],[],[],[]
        with open(self.tab,'r') as handle:
            for ln in handle:
                toks = ln.rsplit()
                title,st,end,ptype = toks[0],int(toks[1]),int(toks[2]),toks[3]
                self.titles.append(title)
                self.starts.append(st)
                self.ends.append(end)
                self.ptypes.append(ptype)

        pass


def result_generator(scopeParserObject):
    """Take a SCOPEparser object and return generator for iteration over the results"""
    for title, start, end, ptype in zip(scopeParserObject.titles, scopeParserObject.starts, scopeParserObject.ends, scopeParserObject.ptypes):
        yield title,start,end,ptype
