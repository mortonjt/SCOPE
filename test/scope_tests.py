"""
Makes sure that all of the output coming from SCOPE is expected
"""
import unittest
import subprocess
import os,site
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "parse"))
from scope_parser import *

if __name__=="__main__":
    class TestHomopolyerDetectionFasta(unittest.TestCase):
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
            os.remove(self.input_fasta)
            os.remove(self.output_fasta)
            os.remove(self.output_tab)
            pass
        def test1(self):
            cmd = "%s/src/scope -i %s -o %s -z"%(self.rootdir,self.input_fasta,self.output_fasta)
            print cmd
            proc = subprocess.Popen(cmd,shell=True)
            proc.wait()
            parser = SCOPEparser(self.output_fasta,self.output_tab)
            self.assertEquals(parser.starts,self.csts)
            self.assertEquals(parser.ends,self.cends)
    class TestHomopolyerDetectionFastq1(unittest.TestCase):
        def setUp(self):
            self.input_fasta = "test.fq"
            seqs = ["@seq1",
                    "GGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    "+seq1",
                    "GGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    "@seq2",
                    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    "+seq2",
                    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    "@seq3",
                    "GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    "+seq3",
                    "GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    "@seq4",
                    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    "+seq4",
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
            cmd = "%s/src/scope -i %s -o %s -z"%(self.rootdir,self.input_fasta,self.output_fasta)
            print cmd
            proc = subprocess.Popen(cmd,shell=True)
            proc.wait()
            parser = SCOPEparser(self.output_fasta,self.output_tab)
            self.assertEquals(parser.starts,self.csts)
            self.assertEquals(parser.ends,self.cends)
        def test2(self):
            cmd = "%s/src/scope -i %s -o %s -z --trim"%(self.rootdir,self.input_fasta,self.output_fasta)
            print cmd
            proc = subprocess.Popen(cmd,shell=True)
            proc.wait()
            parser = SCOPEparser(self.output_fasta,self.output_tab)
            self.assertEquals(parser.starts,self.csts)
            self.assertEquals(parser.ends,self.cends)
        def test3(self):
            cmd = "%s/src/scope -i %s -o %s -z --numThreads 4"%(self.rootdir,self.input_fasta,self.output_fasta)
            print cmd
            proc = subprocess.Popen(cmd,shell=True)
            proc.wait()
            parser = SCOPEparser(self.output_fasta,self.output_tab)
            self.assertEquals(set(parser.starts),set(self.csts))
            self.assertEquals(set(parser.ends),set(self.cends))
        def test4(self):
            cmd = "%s/src/scope -i %s -o %s -z --trim --numThreads 4"%(self.rootdir,self.input_fasta,self.output_fasta)
            print cmd
            proc = subprocess.Popen(cmd,shell=True)
            proc.wait()
            parser = SCOPEparser(self.output_fasta,self.output_tab)
            self.assertEquals(set(parser.starts),set(self.csts))
            self.assertEquals(set(parser.ends),set(self.cends))
        def test5(self):
            cmd = "%s/src/scope -i %s -o %s -z --mask --numThreads 4"%(self.rootdir,self.input_fasta,self.output_fasta)
            print cmd
            proc = subprocess.Popen(cmd,shell=True)
            proc.wait()
            parser = SCOPEparser(self.output_fasta,self.output_tab)
            self.assertEquals(set(parser.starts),set(self.csts))
            self.assertEquals(set(parser.ends),set(self.cends))

    class TestHomopolyerDetectionFastq2(unittest.TestCase):
        def setUp(self):
            self.input_fasta = "test.fq"
            seqs = ["@seq1",
                    "GGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    "+seq1",
                    "GGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    "@seq2",
                    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    "+seq2",
                    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    "@seq3",
                    "TTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
                    "+seq3",
                    "TTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
                    "@seq4",
                    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    "+seq4",
                    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"]
            self.csts  = [20,40,0,30]
            self.cends = [75,75,20,75]
            self.casts  = [20,40,30]
            self.caends = [75,75,75]
            self.ctsts  = [0]
            self.ctends = [20]
            handle = open(self.input_fasta,'w')
            handle.write('\n'.join(seqs))
            handle.close()
            self.output_fasta = "out.fa"
            self.output_tab = "out.fa.tab"
            self.rootdir = ".."
        def test1(self):
            cmd = "%s/src/scope -i %s -o %s -z --mask --numThreads 4"%(self.rootdir,self.input_fasta,self.output_fasta)
            print cmd
            proc = subprocess.Popen(cmd,shell=True)
            proc.wait()
            parser = SCOPEparser(self.output_fasta,self.output_tab)
            self.assertEquals(set(parser.starts),set(self.csts))
            self.assertEquals(set(parser.ends),set(self.cends))

        def test2(self):
            cmd = "%s/src/scope -i %s -o %s -z --homopolymer_type A --numThreads 4"%(self.rootdir,self.input_fasta,self.output_fasta)
            print cmd
            proc = subprocess.Popen(cmd,shell=True)
            proc.wait()
            parser = SCOPEparser(self.output_fasta,self.output_tab)
            self.assertEquals(set(parser.starts),set(self.casts))
            self.assertEquals(set(parser.ends),set(self.caends))
        def test3(self):
            cmd = "%s/src/scope -i %s -o %s -z --homopolymer_type T --numThreads 4"%(self.rootdir,self.input_fasta,self.output_fasta)
            print cmd
            proc = subprocess.Popen(cmd,shell=True)
            proc.wait()
            parser = SCOPEparser(self.output_fasta,self.output_tab)
            self.assertEquals(set(parser.starts),set(self.ctsts))
            self.assertEquals(set(parser.ends),set(self.ctends))


    class TestHomopolyerDetectionFasta3(unittest.TestCase):
        def setUp(self):
            self.input_fasta = "test.fa"
            seqs = [">seq1\n",
                    "GGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    ">seq2\n",
                    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n\n",
                    ">seq3\n",
                    "GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n\n",
                    ">seq4\n",
                    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n\n",
                    ">testing\n",
                    "atgcggtaggcccttgatgcggtaggcccttg",
                    ">testing2\n",
                    "aatgcggtaggcccttatgcggtaggcccttggaaaaaaaaaaaaaaaaa\n\n",
            ]
            self.csts  = [20,40,10,30,33]
            self.cends = [75,75,75,75,50]
            handle = open(self.input_fasta,'w')
            handle.write('\n'.join(seqs))
            handle.close()
            self.output_fasta = "out.fa"
            self.output_tab = "out.fa.tab"
            self.rootdir = ".."
        def tearDown(self):
            os.remove(self.input_fasta)
            os.remove(self.output_fasta)
            os.remove(self.output_tab)
            pass
        def test1(self):
            cmd = "%s/src/scope -i %s -o %s -z"%(self.rootdir,self.input_fasta,self.output_fasta)
            print cmd
            proc = subprocess.Popen(cmd,shell=True)
            proc.wait()
            parser = SCOPEparser(self.output_fasta,self.output_tab)
            self.assertEquals(parser.starts,self.csts)
            self.assertEquals(parser.ends,self.cends)

    unittest.main()


