"""
Makes sure that all of the output coming from SCOPE++ is expected
"""
import unittest
import subprocess
import os,site
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "parse"))
from scope_parser import *

if __name__=="__main__":
    class TestHomopolyerDetection(unittest.TestCase):
        def setUp(self):
            self.input_fasta = "test.fa"
            seqs = [">seq1",
                    "GGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    ">seq2",
                    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    ">seq3",
                    "GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    ">seq4",
                    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"]
            self.csts  = [20,40,10,30]
            self.cends = [75,75,75,75]
            handle = open(self.input_fasta,'w')
            handle.write('\n'.join(seqs))
            handle.close()
            self.output_fasta = "out.fa"
            self.output_tab = "out.fa.tab"
            self.rootdir = ".."
        def tearDown(self):
            #os.remove(self.input_fasta)
            #os.remove(self.output_fasta)
            #os.remove(self.output_tab)
            pass
        def test1(self):
            cmd = "%s/src/SCOPE++ -i %s -o %s -z"%(self.rootdir,self.input_fasta,self.output_fasta)
            print cmd
            proc = subprocess.Popen(cmd,shell=True)
            proc.wait()
            parser = SCOPEparser(self.output_fasta,self.output_tab)
            self.assertEquals(parser.starts,self.csts)
            self.assertEquals(parser.ends,self.cends)

    unittest.main()


