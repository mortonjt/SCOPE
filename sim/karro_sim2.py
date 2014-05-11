#!/usr/bin/python2.7

from Bio import SeqIO
from optparse import OptionParser
from optparse import OptionGroup
from cleaningTools import *
from random import *
import subprocess
import shlex
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIStandalone
import os
import stat
import scipy
import sys
from re import *
from numpy import arange

bases = "ACGT";

def pickFromDistribution(CD):
    r = random()
    return min([i for i,v in enumerate(CD) if r < v])

def randomBase(CD):
    return bases[pickFromDistribution(CD)]

def randomSequence(len, CD):
    return "".join([randomBase(CD) for i in range(len)])

class ESTFrag:
    """Helper class for working with the fragment sequences.
       A fragment consists of three parts:
       * The main fragement.  Length and base distribution specified in parameters.
       * The polyA tail.  Tail is generated and the subject to errors.
       * The "tail" -- any sequence that might come after the polyA tail.
       The fragment is randomly generated in the constructor.  It is mainted as three different sequences, with the
       sequence function returning the concetenation of the three.
       """
    def __init__(self, label, head, polyA, tail):
        self.label = label
        self.head_seq = head
        self.polyA_seq = polyA
        self.tail_seq = tail

    def __str__(self):
        return self.sequence()


    def sequence(self):
        """Returns the entire seqeunces."""
        return self.head_seq + self.polyA_seq + self.tail_seq

    def fastaString(self, print_len = 60):
        """Returns the sequence as a fasta entry"""
        s = self.sequence()

        return ">id%s|frag:%d|polyA:%d|extra:%d" % (self.label, len(self.head_seq), len(self.head_seq) + len(self.polyA_seq), \
                                                    len(s)) + "\n" + "\n".join([s[i:i+60] for i in range(0, len(s), 60)]) + "\n"
    


def stats(V):
    m = (1.0*sum(V)) / len(V)
    v = sum( [(x-m)*(x-m) for x in V] ) / (len(V) - 1.0)
    return (m,pow(v,0.5))


info_re = re.compile("\((\d+)\s*,\s*(\d+)\)$")
def simulation(outputFile, n, DIR, input_file, error_rate = None, block_size = 3, polyType = "A"): 
    """Create a file of simulated sequences"""

    nonPolyType = [x for x in ['A','C','G','T'] if x != polyType]

    # First: Read in actual sequences to a fragment list O

    counter = 0
    O = []
    for r in SeqIO.parse(input_file, "fasta"):
        s = str(r.seq)
        r = info_re.search(r.name)
        if r:
            x = int(r.group(1)) - 1
            y = int(r.group(2))
            O.append(ESTFrag(label = counter, head = s[:x], polyA = s[x:y], tail = s[(y+1):]))
        else:
            O.append(ESTFrag(label = counter, head = s, polyA = "", tail = ""))
        counter = counter + 1
    
    # Second: Computate relevant statistics
 
    
    if not error_rate:
        numA = sum([o.polyA_seq[block_size:-block_size].count(polyType) for o in O])
        total = sum([len(o.polyA_seq[block_size:-block_size]) for o in O])
        error_rate = 1.0 - (1.0*numA)/total;

    # Third: generate and print sequences
    fp = open(DIR + "/" + outputFile, "w")
    
    for index in range(n):        
        # Generate head by sampling by blocks
        ref_seq = choice(O)

        head_len = len(ref_seq.head_seq)
        if head_len < block_size:
            head = ref_seq.head_seq
        else:
            head = ""
            O_head = [o for o in O if len(o.head_seq) > block_size]
            while len(head) < head_len:
                o = choice(O_head)
                position = choice(range(len(o.head_seq) - block_size))
                head = head + o.head_seq[position:(position + block_size)]
            head = head[:head_len][:-block_size] + ref_seq.head_seq[-block_size:]

        # Generate polyA tail    # Improvement: Make error_rate position dependent?
        polyA_len = len(ref_seq.polyA_seq)
        if polyA_len < 2*block_size:
            polyA = ref_seq.polyA_seq
        else:
            polyA = ref_seq.polyA_seq[:block_size] + "".join([choice(nonPolyType) if random() < error_rate else polyType for i in range(polyA_len - 2*block_size)]) + ref_seq.polyA_seq[-block_size:]


        # Generate tail by sampling by blocks
        tail_len = len(ref_seq.tail_seq)
        if tail_len < block_size:
            tail = ref_seq.tail_seq
        else:
            tail = ""
            block_size = min(block_size, tail_len)
            O_tail = [o for o in O if len(o.tail_seq) > block_size]
            while len(tail) < tail_len:
                o = choice(O_tail)
                position = choice(range(len(o.tail_seq)))
                tail = tail + o.tail_seq[position:(position + block_size)]
            tail = ref_seq.tail_seq[:block_size] + tail[tail_len:][block_size:]
        
        fp.write(ESTFrag(str(index), head, polyA, tail).fastaString())
        
    fp.close()

def real_data(outputFile, DIR, input_file):
    counter = 0
    O = []
    for r in SeqIO.parse(input_file, "fasta"):
        s = str(r.seq)
        r = info_re.search(r.name)
        if r:
            x = int(r.group(1)) - 1
            y = int(r.group(2))
            O.append(ESTFrag(label = counter, head = s[:x], polyA = s[x:y], tail = s[(y+1):]))
        else:
            O.append(ESTFrag(label = counter, head = s, polyA = "", tail = ""))
        counter = counter + 1


    fp = open(DIR + "/" + outputFile, "w")
    for o in O:
        fp.write(o.fastaString())

    fp.close()

column_names = ["name", "run_time", "sensitivity", "specificity", "%correct_trim", "avg_trim", "median_trim", "avg_left", "med_left", "avg_right", "median_right", "SoS_trim", "pct_over", "avg_over", "med_over", "pct_under", "avg_under", "median_under"];
def analysis(records, analysis_function):
    """Read in the results of one of the tools and calculate certain statistics.
       fn is a function for reading in the results (e.g. cleangingTools.parseScopaInfo)
       """
    trim = []
    left_trim = []
    right_trim = []

    tp, fp, tn, fn = [0]*4   # true positive, false positive, ...
    for i,seq_record in enumerate(records):
        id, present, actual_start, actual_end, found, predicted_start, predicted_end = analysis_function(seq_record)[-1]

        if present:
            if found:
                tp = tp + 1
                left_trim.append(int(actual_start) - int(predicted_start))
                right_trim.append(int(predicted_end) - int(actual_end))
                trim.append(left_trim[-1] + right_trim[-1])
            else:
                fn = fn + 1
        else:
            if found:
                fp = fp + 1
            else:
                tn = tn + 1

    sensitivity = float(tp) / (tp + fn) if tp + fn > 0 else -1
    specificity = float(tn) / (tn + fp) if tn + fp > 0 else -1
    if len(trim) > 0:
        pct_correct = len(filter(lambda x: x==0, trim)) / float(len(trim))
        avg_trim = scipy.mean(trim)
        median_trim = scipy.median(trim)
        SoS_trim = scipy.mean(map(lambda x : x*x, trim))
        avg_left = scipy.mean(left_trim)
        median_left = scipy.median(left_trim)
        avg_right = scipy.mean(right_trim)
        median_right = scipy.median(right_trim)

        overArr = filter(lambda x : x > 0, trim)
        pct_over = len(overArr) / float(len(trim))
        avg_over = scipy.mean(overArr) if len(overArr) > 0 else -99999
        median_over = scipy.median(overArr) if len(overArr) > 0 else -99999

        underArr = filter(lambda x : x < 0, trim)
        pct_under = len(underArr) / float(len(trim))
        avg_under = scipy.mean(underArr) if len(underArr) > 0 else -99999
        median_under = scipy.median(underArr) if len(underArr) > 0 else -99999
    else:
        return [sensitivity, specificity] + [9999]*10

    return [sensitivity, specificity, pct_correct, avg_trim, median_trim, avg_left, median_left, avg_right, median_right, SoS_trim, pct_over, avg_over, median_over, pct_under, avg_under, median_under]

def printResults(out, name, results):
    arr = results
    out.write('{:<14}'.format(name))
    out.write(('{:<14}'*len(results)).format(*[float("%0.3f" % v) for v in arr]) + "\n")

def basicSim(options, args):
    # Generate files
    if (options.real):
        real_data(*[getattr(options,v) for v in ["sim_file", "outputDir", "real"]])
    else:
        simulation(*[getattr(options,v) for v in ["sim_file", "n", "outputDir", "basis", "polyType"]]);

    # Launch simulations
    SCOPA = scopaJob(inputFile = options.sim_file, e = "", t = options.polyType, w = 8, d = None, m = 10, b = 6, x = 1, n = 10000, k=0, outputFile = "SCOPA.out", DIR=options.outputDir) if options.SCOPA else None
    TRIM = seqTrimJob(options.sim_file, "TRIM.out", DIR = options.outputDir) if options.TRIM else None
    CLEAN = seqCleanJob(options.sim_file, "CLEAN", DIR = options.outputDir) if options.CLEAN else None

    # Get / print results in column format
    fp = sys.stdout if options.output == None else open(options.output, "w")
    if not options.Header_off:
        fp.write(('{:<14}'*len(column_names)).format(*column_names) + "\n")

    if options.SCOPA:
        time, results = scopaCollect(SCOPA)
        stats = analysis(results, parseScopaInfo)
        printResults(fp, "SCOPA", [time] + stats);
    if options.CLEAN:
        time, results = seqCleanCollect(CLEAN)
        stats = analysis(results, parseCleanInfo)
        printResults(fp, "CLEAN", [time] + stats)
    if options.TRIM:
        time, results = seqTrimCollect(TRIM)
        stats = analysis(results, parseTrimInfo)
        printResults(fp, "TRIM", [time] + stats)    


def crossValidation(options, args):
    simulation(*[getattr(options,v) for v in ["sim_file", "n", "outputDir", "basis", "polyType"]]);
    fp = sys.stdout if options.output == None else open(options.output, "w")
    
    if not options.Header_off:
        fp.write(('{:<7}'*4).format('w', 'k', 'b', 'x') + ('{:<14}'*len(column_names)).format(*column_names) + "\n")
    
    currentJobs = [seqCleanJob(options.sim_file, "CLEAN", DIR = options.outputDir)]

    #for w in range(5,20):
    #    for k in range(-100,100,10):
    #        for e in range(0,5):                
    id = 0
    for w1 in range(5,20):
        print "w1: ", w1        
        for k1 in range(0,15):
            for b1 in range(0,w1,2):
                for x1 in range(0,10):
                    job = scopaJob(options.sim_file, e = "", t=options.polyType, d=None, w=w1, k=k1, b=b1, x=x1, outputFile="SCOPA.%d.out" % (id), DIR=options.outputDir)
                    id, job.w, job.k, job.b, job.x = id + 1, w1, k1, b1, x1
                    currentJobs.append(job)

    time, results = seqCleanCollect(currentJobs[0])
    stats = analysis(results, parseCleanInfo)
    fp.write(('{:<7}'*4).format("-1", "-1", "-1", "-1"))
    printResults(fp, "CLEAN", [time] + stats);

    for job in currentJobs[1:]:
        time, results = scopaCollect(job)
        stats = analysis(results, parseScopaInfo)
        fp.write(('{:<7}'*4).format(job.w, job.k, job.b, job.x))
        printResults(fp, "SCOPA", [time] + stats);

def error_variation(options, args):

    fp = sys.stdout if options.output == None else open(options.output, "w")
    
    if not options.Header_off:
        fp.write(('{:<7}'*2).format('e','x') + ('{:<14}'*len(column_names)).format(*column_names) + "\n")
    
    currentJobs = []

    id = 0
    w,k,b = 8,0,1
    for err in arange(*options.error_range):
        sim_file = options.sim_file + "." + str(id)
        simulation(outputFile = sim_file, n = options.n, DIR = options.outputDir, input_file = options.basis, error_rate = err, polyType = options.polyType)

        Sjobs = []
        for x in range(0, 7, 1):
            Sjob = scopaJob(inputFile=sim_file, e = "", t=options.polyType, d=None, w=w, k=k, b=b, x=x, outputFile="SCOPA.%d.out" % (id), DIR=options.outputDir)
            Sjob.err = err
            Sjob.x = x
            Sjobs.append(Sjob)

        Cjob = seqCleanJob(inputFile=sim_file, outputFile = "CLEAN", DIR = options.outputDir, id=id)
        Cjob.err = err
        Cjob.x = -1

        Cjob.wait_on_job()

        currentJobs.extend(Sjobs + [Cjob])
        id = id + 1

    for job in currentJobs:
        if job.name == "SCOPA":
            time, results = scopaCollect(job)
            stats = analysis(results, parseScopaInfo)
            fp.write(('{:<7}'*2).format(job.err, job.x))
            printResults(fp, "SCOPA", [time] + stats);
        elif job.name == "CLEAN":
            time, results = seqCleanCollect(job)
            stats = analysis(results, parseCleanInfo)
            fp.write(('{:<7}'*2).format(job.err, job.x))
            printResults(fp, "CLEAN", [time] + stats);
        else:
            print "BAD JOB NAME: %s" % job.name
            exit(1)

        


        

def tuneK(options,args):
    simulation(*[getattr(options,v) for v in ["sim_file", "n", "frag_len", "base_distribution", "polyA_prob", "polyA_stats", "error_rate", "insert_stats", "delete_prob", "tail_stats", "outputDir"]]);
#    sys.stdout.write(('{:<14}'*len(column_names)).format(*column_names) + "\n")
    CLEAN = seqCleanJob(options.sim_file, "CLEAN", DIR = options.outputDir)
    time, results = seqCleanCollect(CLEAN)
    stats = analysis(results, parseCleanInfo)
    printResults(sys.stdout, "CLEAN", [time] + stats)
    for k in range(-100,100,10):
        SCOPA = scopaRun(options.sim_file,options.polyType,10,10,k,0,"SCOPA.out", DIR=options.outputDir)
        time, results = scopaCollect(SCOPA)
        stats = analysis(results, parseScopaInfo)
        printResults(sys.stdout, "SCOPA", [time] + stats);


if __name__ == "__main__":
    # Parse command-line arguments
    parser = OptionParser()
    simGroup = OptionGroup(parser, "Simulation Switches") 
    parser.add_option_group(simGroup)
    simGroup.add_option("-r", "--real", action="store", type="string", dest="real", help="Use real data from specified file", default=None) 
    simGroup.add_option("-e", "--error_rate", action="store", type="float", dest="error_rate", help="base call error rate (default: %default)", default=0.01)
    simGroup.add_option("-n", "--num_sequences", action="store", type="int", dest="n", help="Number of simulated sequences (default: %default)", default=20000)
    simGroup.add_option("-s", "--sim_file", action="store", type="string", dest="sim_file", help="Output file (default: %default)", default = "sim.fa")
    simGroup.add_option("-b", "--basis", action="store", type="string", dest="basis", help="Basis file for simulation (default: %default)", default = "human_annotated.fa")

    toolGroup = OptionGroup(parser, "Tool Switches")
    parser.add_option_group(toolGroup)
    toolGroup.add_option("-t", "--polyType", action="store", type = "string", dest="polyType", help="homopoloymer type (default: $default)", default = "A")
    toolGroup.add_option("-S", "--SCOPA_off", action="store_false", dest = "SCOPA", help="Supress SCOPA testing (default: %default)", default = True)
    toolGroup.add_option("-T", "--TRIM_off", action="store_true", dest = "TRIM", help="Supress SeqTrim testing (default: %default)", default = False)
    toolGroup.add_option("-C", "--CLEAN_off", action="store_false", dest = "CLEAN", help="Supress SqClean testing (default: %default)", default = True)
    toolGroup.add_option("--DIRECTORY", action="store", type="string", dest="outputDir", metavar="DIR", help="Directory for results (default: %default)", default = "results_files")
    toolGroup.add_option("-c", "--cross_validate", action="store_const", const = "crossValidation", dest="run_type", help="Preform cross validation", default = "basicSim")
    toolGroup.add_option("--tune", action="store_const", const = "tuneK", dest="run_type", help="Perform tuning", default = "basicSim")
    toolGroup.add_option("-E", "--error_variation", action="store_const", const = "error_variation", dest="run_type", help="Vary simulation error rate", default = "basicSim")
    toolGroup.add_option("-R", "--error_range", action="store", type="float", dest="error_range", nargs=3, metavar = "MIN MAX STEP", help = "Vary simulation error", default = [0, 0.101, 0.01])
    toolGroup.add_option("--seed", action="store", dest="seed", help="Set rng seed (default: None)", default = None)

    outputGroup = OptionGroup(parser, "Output Format Switches")
    parser.add_option_group(outputGroup)
    outputGroup.add_option("-o", "--output", action="store", type="string", dest="output", help="Output file (default: stdout)", default = None)
    outputGroup.add_option("-H", "--Header_off", action="store_true", dest = "Header_off", help="Turn reader row off", default = False)

    (options, args) = parser.parse_args()

    if options.seed:
        seed(options.seed)

    eval(options.run_type)(options,args)



