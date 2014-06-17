#!/usr/bin/python2.7

from redhawk import *
from Bio import SeqIO
from optparse import OptionParser
from optparse import OptionGroup
from cleaningTools import *
from random import *
import subprocess
import shlex
from Bio.Blast import NCBIXML
#from Bio.Blast import NCBIStandalone
import os
import stat
import scipy
import sys
from re import *
from numpy import arange
import glob
import shutil

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
    def __init__(self, label, seq_arr, seq_type):
        self.label = label
        self.seq_arr = seq_arr
        self.seq_type = seq_type

    def __str__(self):
        return self.sequence()

    def sequence(self):
        """Returns the entire seqeunces."""
        return "".join(self.seq_arr)

    def fastaString(self, print_len = 60):
        """Returns the sequence as a fasta entry"""
        s = self.sequence().upper()

        return ">id%s" % (self.label) + "".join(["|%s:%d" % (b, len(a)) for a,b in zip(self.seq_arr, self.seq_type)]) + "\n" + "\n".join([s[i:i+60] for i in range(0, len(s), 60)]) + "\n"

    


def stats(V):
    m = (1.0*sum(V)) / len(V)
    v = sum( [(x-m)*(x-m) for x in V] ) / (len(V) - 1.0)
    return (m,pow(v,0.5))


info_re = re.compile("\((\d+)\s*,\s*(\d+)\)$")
def max_letter(s):
    return max([(s.count(c), c) for c in ['A', 'C', 'G', 'T']])[1]


def id_parser(id):
    seq_info = []
    current = 0
    for p in re.findall("(frag|poly\w):(\d+)", id):
        if p[0] == "frag":
            current += int(p[1])
        else:
            t = re.search("poly(\w+)", p[0])
            seq_info.append((t.group(1), current, current + int(p[1])))
            current += int(p[1])

    return seq_info, current

def read_validation_data(input_file, options):
    """Read in and parse data from a specified validation file"""
    O = []
    for mult in range(options.multiplier):
        for counter,r in enumerate(SeqIO.parse(input_file, "fasta")):
            if options.sim_size > 0 and counter >= options.sim_size:
                break
            s = str(r.seq)

            # Induce simulated base-call error in tails
            if options.error_rate >= 0:
                base_set = list(set("ACGT") - set(options.polyType))

                for p in re.findall("\((?:homopolymer=)?(\w)(?:\,|\:)(?:start=)?(\d+),(?:end=)?(\d+)(?:,.*)?\)", r.description):
                    x, y = int(p[1])-1, int(p[2])
                    s = s[:x] + "".join([options.polyType if random() >= options.error_rate else choice(base_set) for c in s[x:y]]) + s[y:]
                #s = "".join([c if random() >= options.error_rate else choice(list(set("ACGT") - set(c))) for c in s])

            #if options.post_tail_extend:
            #    extenion = "".join([choice("ACGT") for i in range(options.post_tail_extend)])
                

            # Break into fragments
            current = 0;
            seq_arr = []
            seq_type = []
            for p in re.findall("\((?:homopolymer=)?(\w)(?:\,|\:)(?:start=)?(\d+),(?:end=)?(\d+)(?:,.*)?\)", r.description):
                type = p[0]
                x = int(p[1]) - 1
                y = int(p[2])
                if (current < x):
                    seq_arr.append(s[current:x])
                    seq_type.append('frag');
                if x < y:
                    seq_arr.append(s[x:y])
                    seq_type.append("poly" + type)
                current = y;

            if current < len(s):
                seq_arr.append(s[current:])
                seq_type.append('frag')

            # Add extend simulated post-fragment length
            if options.post_tail_extend != None:
                if options.polyType == 'A':
                    if seq_type[-1] == 'frag':
                        seq_arr = seq_arr[:-1]
                        seq_type = seq_type[:-1]
                    if options.post_tail_extend > 0:
                        seq_arr.append("".join([choice("ACGT") for i in range(options.post_tail_extend)]))
                        seq_type.append('frag')
                else:
                    if seq_type[0] == 'frag':
                        seq_arr = seq_arr[1:]
                        seq_type = seq_type[1:]
                    if options.post_tail_extend > 0:
                        seq_arr = ["".join([choice("ACGT") for i in range(options.post_tail_extend)])] + seq_arr
                        seq_time = ['frag'] + seq_type

            # Change length of polymer
            if options.poly_len:
                if options.polyType == 'T':
                    i = 0 if seq_type[0] == 'polyT' else (1 if seq_type[1] == 'polyT' else -1)
                elif options.polyType == 'A':
                    n = len(seq_type)
                    i = n-1 if seq_type[n-1] == 'polyA' else (n-2 if seq_type[n-2] == 'polyA' else -1)

                if i >= 0:
                    seq_arr[i] = "".join([choice(seq_arr[i]) for j in range(options.poly_len)])
                
                        

            O.append(ESTFrag(label=str(counter) + ".a", seq_arr = seq_arr, seq_type = seq_type))    
            #print(O[-1].label, O[-1].seq_type, O[-1].seq_arr)

            if options.specificity_test:
                last = 0
                tmp = ""
                for p in re.findall("\((?:homopolymer=)?(\w)(?:\,|\:)(?:start=)?(\d+),(?:end=)?(\d+)(?:,.*)?\)", r.description):
                    tmp += r.seq[last:int(p[1])]
                    last = int(p[2])
                tmp += r.seq[last:-1]
                s = "".join([choice(tmp) for i in range(len(r.seq))])
                O.append(ESTFrag(label=str(counter) + ".b", seq_arr = [s], seq_type = ['frag']))

        #shuffle(O)

    return O

def real_data(outputFile, DIR, inputFile, options):
    """Convert a validation file to the needed format"""
    with open(DIR + "/" + outputFile, "w") as fp:
        for o in read_validation_data(inputFile, options):
            fp.write(o.fastaString())

def timing_data(outputFile, DIR, inputFile, options, n):
    """Create an input file from a validation sequence where we choose n random sequences (with replacement)
    for purposes of checking runtime"""
    with open(DIR + "/" + outputFile, "w") as fp:
        O = read_validation_data(inputFile, options)
        seqs = []
        IDmap = {}
        for i in range(n):
            s = choice(O).fastaString()
            id = search("id(\d+)", s).group(1)
            index = IDmap[id]+1 if IDmap.has_key(id) else 1
            IDmap[id] = index
            s = sub("id(%s)" % (id), "id%s.%d" % (id, index), s)
            fp.write(s)



#column_names = ["name", "run_time", "sens", "spec", "split", "join", "pct_crt_trim", "avg_trim", "med_trim", "avg_left", "med_left", "avg_right", "med_right", "SoS_trim", "pct_over", "avg_over", "med_over", "pct_under", "avg_under", "med_under"]
column_names = ["name", "run_time", "sens", "spec", "pct_crt", "avg_trim", "SoS_trim", "pct_crt_left", "avg_left", "SoS_left", "pct_crt_right", "avg_right", "SoS_right", "total_polyA", "num_split"]
def hasOverlap(i1, i2):
    return not (i1[2] < i2[1] or i2[2] < i1[1])


def mean(A, error_val = -1):
    return scipy.mean(A) if len(A) > 0 else error_val

def analysis(records, analysis_function, polyType, options):
    """Read in the results of one of the tools and calculate certain statistics.
       fn is a function for reading in the results (e.g. cleangingTools.parseScopaInfo)
       """
    trim = []
    left_trim = []
    right_trim = []
    split = 0
    join = 0
    tp, fp, tn, fn = [0]*4   # true positive, false positive, ...

    for i,seq_entry in enumerate(records):
        id, real_list, predicted_list, seq_len = analysis_function(seq_entry, id_parser)
        r, p = 0, 0

        # remove wrong poly type, and polyTypes too far from an edge:
        real_list = [entry for entry in real_list if entry[0] == polyType and
                          ((entry[0] == "T" and entry[1] <= options.poly_dist) or
                          (entry[0] == "A" and seq_len - entry[2] <= options.poly_dist))];
        predicted_list = [entry for entry in predicted_list if entry[0] == polyType and
                          ((entry[0] == "T" and entry[1] <= options.poly_dist) or
                           (entry[0] == "A" and seq_len - entry[2] <= options.poly_dist))];
        

        
        if len(real_list) > 1:
            continue
        if len(real_list)==0:
            if len(predicted_list)==0:
                tn += 1
            else:
                fp += 1
        else:
            if len(predicted_list)==0:
                fn += 1
            elif len(predicted_list) > 1:
                split += 1
            else:
                #print(real_list, predicted_list)
                tp += 1
                left_trim.append(real_list[0][1] - predicted_list[0][1])
                right_trim.append(predicted_list[0][2] - real_list[0][2])
                trim.append(left_trim[-1] + right_trim[-1])


    sensitivity = float(tp) / (tp + fn) if tp + fn > 0 else -1
    specificity = float(tn) / (tn + fp) if tn + fp > 0 else -1

    total_pos = tp + fn
    total_neg = fp + tn
    if len(trim) > 0:
        pct_correct = len(list(filter(lambda x: x==0, trim))) / float(len(trim))
        avg_trim = mean(trim)
        median_trim = scipy.median(trim)
        SoS_trim = mean([t*t for t in trim if t != 0]) / len([t for t in trim]) # if t != 0])
        pct_correct_left = len(list(filter(lambda x: x==0, left_trim))) / float(len(left_trim))
        avg_left = mean(left_trim)

        SoS_left = mean([t*t for t in left_trim if t != 0]) / len([t for t in left_trim])
        median_left = scipy.median(left_trim)
        pct_correct_right = len(list(filter(lambda x: x==0, right_trim))) / float(len(right_trim)) if float(len(right_trim)) else -1
        avg_right = mean(right_trim)
        SoS_right = mean([t*t for t in right_trim if t != 0]) / len([t for t in right_trim]) # if t != 0])
        median_right = scipy.median(right_trim)

        overArr = list(filter(lambda x : x > 0, trim))
        pct_over = len(overArr) / float(len(trim))
        avg_over = mean(overArr) if len(overArr) > 0 else -99999
        median_over = scipy.median(overArr) if len(overArr) > 0 else -99999

        underArr = list(filter(lambda x : x < 0, trim))
        pct_under = len(underArr) / float(len(trim))
        avg_under = mean(underArr) if len(underArr) > 0 else -99999
        median_under = scipy.median(underArr) if len(underArr) > 0 else -99999

    else:
        return [sensitivity, specificity] + [9999]*(len(column_names)-4)

    #return [sensitivity, specificity, float(split)/tp, float(join)/tp, pct_correct, avg_trim, median_trim, avg_left, median_left, avg_right, median_right, SoS_trim, pct_over, avg_over, median_over, pct_under, avg_under, median_under]
    return [sensitivity, specificity, pct_correct, avg_trim, SoS_trim, pct_correct_left, avg_left, SoS_left, pct_correct_right, avg_right, SoS_right, total_pos, float(split)/(tp+fn)]

def printResults(out, name, results):
    arr = results
    out.write('{:<14}'.format(name))
    out.write(('{:<14}'*len(results)).format(*[float("%0.3f" % v) for v in arr]) + "\n")


def launch_jobs(options):
    if not options.last_run:        
        if (options.Liang_data):
            shutil.copyfile(options.Liang_data, options.outputDir + "/" + options.sim_file)
        elif (options.real):
            real_data(*[getattr(options,v) for v in ["sim_file", "outputDir", "real"]], options=options)
        else:
            simulation(*[getattr(options,v) for v in ["sim_file", "n", "outputDir", "basis", "polyType"]]);
    
        m = None
        w,b,x = options.params;

        # Launch simulations
        fp = open("log.txt", "a")
        SCOPA   = scopaJob(inputFile = options.sim_file, e = None, t = options.polyType, w =w, d = None, m = m, b = b, x = x, n = None, k=None, outputFile = "SCOPA.out.%s" % options.id, DIR=options.outputDir, 
                           minIdentity=options.minIdentity, terminate=options.terminate, no_retrain=True, front_gap=options.front_gap, poly=options.poly, numTrain=options.numTrain) if options.SCOPA else None
        SCOPABW = scopaJob(inputFile = options.sim_file, e = "", t = options.polyType, w =w, d = None, m = m, b = b, x = x, n = None, k=None, outputFile = "SCOPA.out.%s" % options.id, DIR=options.outputDir, 
                         minIdentity=options.minIdentity, terminate=options.terminate, no_retrain=False, front_gap=options.front_gap, poly=options.poly, numTrain=options.numTrain) if options.SCOPABW else None
        #SCOPABW = scopaJob(inputFile = options.sim_file, e = "", homopolymer_type = options.polyType, w =w, d = None, m = m, b = b, x = x, n = None, outputFile = "SCOPA.out.%s" % options.id, DIR=options.outputDir, 
        #                 terminate=options.terminate, minIdentity = options.minIdentity, no_retrain=False, left_gap=options.front_gap if options.polyType=='T' else None, right_gap = options.front_gap if options.polyType=='A' else None, poly=options.poly, version = "17", t = None, front_gap = None) if options.SCOPABW else None
        TRIM = seqTrimJob(options.sim_file, "TRIM.out", DIR = options.outputDir, terminate=options.terminate) if options.TRIM else None
        CLEAN = seqCleanJob(options.sim_file, "CLEAN", DIR = options.outputDir, terminate=options.terminate) if options.CLEAN else None
        POLY = polyJob(options.sim_file, "POLY", DIR = options.outputDir, terminate=options.terminate) if options.POLY else None
        TRIMEST = trimestJob(options.sim_file, "TRIMEST", DIR = options.outputDir, terminate=options.terminate) if options.TRIMEST else None
        BASICTOOL1 = basicJob(options.sim_file, outputFile = "BASIC.out.%s" % options.id, DIR = options.outputDir, terminate=options.terminate, p=0.0) if options.BASICTOOL else None
        BASICTOOL3 = basicJob(options.sim_file, outputFile = "BASIC.out.%s" % options.id, DIR = options.outputDir, terminate=options.terminate, p=0.05) if options.BASICTOOL else None
        BASICTOOL1 = basicJob(options.sim_file, outputFile = "BASIC.out.%s" % options.id, DIR = options.outputDir, terminate=options.terminate, p=0.1) if options.BASICTOOL else None

        #for o in [SCOPA, SCOPABW, TRIM, CLEAN, POLY, TRIMEST]:
        #    if o:
        #        fp.write("%s: %s\n" % (o.jobid, o.executable_name))
        fp.close()

        
        # if options.SCOPA:
        #     storePBS([SCOPA], open("last.scopa", "w"))
        # if options.SCOPABW:
        #     storePBS([SCOPABW], open("last.scopabw", "w"))
        # if options.TRIM:
        #     storePBS([TRIM], open("last.trim", "w")) 
        # if options.CLEAN:
        #     storePBS([CLEAN], open("last.clean", "w")) 
        # if options.POLY:
        #     storePBS([POLY], open("last.poly", "w"))
        # if options.TRIMEST:
        #     storePBS([TRIMEST], open("last.trimest", "w"))
        # if options.BASICTOOL:
        #     storePBS([BASICTOOL], open("last.basictool", "w"))

    else:
        SCOPA = loadPBS(open("last.scopa"))[0] if options.SCOPA else None
        SCOPABW = loadPBS(open("last.scopabw"))[0] if options.SCOPABW else None
        TRIM = loadPBS(open("last.trim"))[0] if options.TRIM else None
        CLEAN = loadPBS(open("last.clean"))[0] if options.CLEAN else None
        POLY = loadPBS(open("last.poly"))[0] if options.POLY else None
        TRIMEST = loadPBS(open("last.trimest"))[0] if options.TRIMEST else None
        BASICTOOL = loadPBS(open("last.basictool"))[0] if options.BASICTOOL else None

    return SCOPA, SCOPABW, TRIM, CLEAN, POLY, TRIMEST, BASICTOOL1, BASICTOOL2, BASICTOOL3


def basicSim(options, args):
    # Generate files

    SCOPA, SCOPABW, TRIM, CLEAN, POLY, TRIMEST = launch_jobs(options)


    # Get / print results in column format
    fp = sys.stdout if options.output == None else open(options.output, "w")
    if not options.Header_off:
        fp.write(('{:<14}'*len(column_names)).format(*column_names) + "\n")

    if options.SCOPA:
        time, results = scopaCollect(SCOPA, not options.keep_files)
        stats = analysis(results, parseScopaInfo, options.polyType, options)
        printResults(fp, "SCOPA", [time] + stats);
    if options.SCOPABW:
        time, results = scopaCollect(SCOPABW, not options.keep_files)
        stats = analysis(results, parseScopaInfo, options.polyType, options)
        printResults(fp, "SCOPABW", [time] + stats);
    if options.CLEAN:
        time, results = seqCleanCollect(CLEAN, not options.keep_files)
        stats = analysis(results, parseCleanInfo, options.polyType, options)
        printResults(fp, "CLEAN", [time] + stats)
    if options.TRIM:
        time, results = seqTrimCollect(TRIM, not options.keep_files)
        stats = analysis(results, parseTrimInfo, options.polyType, options)
        printResults(fp, "TRIM", [time] + stats)  
    if options.POLY:
        time, results = polyCollect(POLY, not options.keep_files)
        stats = analysis(results, parsePolyInfo, options.polyType, options)
        printResults(fp, "POLY", [time] + stats)
    if options.TRIMEST:
        time, results = trimestCollect(TRIMEST, not options.keep_files)
        stats = analysis(results, parseTrimestInfo, options.polyType, options)
        printResults(fp, "TRIMEST", [time] + stats)
    if options.BASICTOOL1:
        time, results = basicCollect(BASICTOOL1, not options.keep_files)
        stats = analysis(results, parseBasicInfo, options.polyType, options)
        printResults(fp, "BASICTOOL1", [time] + stats)
    if options.BASICTOOL2:
        time, results = basicCollect(BASICTOOL2, not options.keep_files)
        stats = analysis(results, parseBasicInfo, options.polyType, options)
        printResults(fp, "BASICTOOL2", [time] + stats)
    if options.BASICTOOL3:
        time, results = basicCollect(BASICTOOL3, not options.keep_files)
        stats = analysis(results, parseBasicInfo, options.polyType, options)
        printResults(fp, "BASICTOOL3", [time] + stats)

def scopeVbasic(options, args):
    
    fp = sys.stdout if options.output == None else open(options.output, "w")
    if not options.Header_off:
        fp.write(('{:<7}'*2).format('p', 'e') + ('{:<14}'*len(column_names)).format(*column_names) + "\n")

    error_start, error_stop, error_step = 0, 0.1001, 0.001
    p_start, p_stop, p_step = 0, 0.10001, 0.01

    for error in arange(error_start, error_stop, error_step):
        options.error_rate = error
        real_data(*[getattr(options,v) for v in ["sim_file", "outputDir", "real"]], options=options)

        if (options.Liang_data):
            shutil.copyfile(options.Liang_data, options.outputDir + "/" + options.sim_file)
        elif (options.real):
            real_data(*[getattr(options,v) for v in ["sim_file", "outputDir", "real"]], options=options)
        else:
            simulation(*[getattr(options,v) for v in ["sim_file", "n", "outputDir", "basis", "polyType"]]);
    
        m = None
        w,b,x = options.params;

        SCOPA   = scopaJob(inputFile = options.sim_file, e = None, t = options.polyType, w =w, d = None, m = m, b = b, x = x, n = None, k=None, outputFile = "SCOPA.out.%s" % options.id, DIR=options.outputDir, 
                           minIdentity=options.minIdentity, terminate=options.terminate, no_retrain=True, front_gap=options.front_gap, poly=options.poly, numTrain=options.numTrain) if options.SCOPA else None

        time, results = scopaCollect(SCOPA, False)
        stats = analysis(results, parseScopaInfo, options.polyType, options)
        fp.write(('{:<7}'*2).format(-1, error))
        printResults(fp, "SCOPA", [time] + stats);

        for p in arange(p_start, p_stop, p_step):
            print(error, p)
            BASICTOOL = basicJob(options.sim_file, outputFile = "BASIC.out.%s" % options.id, DIR = options.outputDir, base_type = options.polyType, min_length = 20, p = p, format = 'fasta', end = '3' if options.polyType == 'A' else '5', terminate=options.terminate) if options.BASICTOOL else None

            time, results = basicCollect(BASICTOOL, not options.keep_files)
            stats = analysis(results, parseBasicInfo, options.polyType, options)
            fp.write(('{:<7}'*2).format(round(p,3), error))
            printResults(fp, "BASICTOOL", [time] + stats)


def error_variation(options, args):
    fp = sys.stdout if options.output == None else open(options.output, "w")
    if not options.Header_off:
        fp.write(('{:<7}'*1).format('e') + ('{:<14}'*len(column_names)).format(*column_names) + "\n")

    start = options.error_range[0]
    stop = options.error_range[1] + options.error_range[2]/2.0
    step = options.error_range[2]

    for error in arange(start, stop, step):
        print("ERROR: %s" % str(error))
        options.error_rate = error
        real_data(*[getattr(options,v) for v in ["sim_file", "outputDir", "real"]], options=options)
        SCOPA, SCOPABW, TRIM, CLEAN, POLY, TRIMEST, BASICTOOL = launch_jobs(options)
        print("LAUNCHED\n")

        if options.SCOPA:
            time, results = scopaCollect(SCOPA, False)
            stats = analysis(results, parseScopaInfo, options.polyType, options)
            fp.write(('{:<7}'*1).format(error))
            printResults(fp, "SCOPA", [time] + stats);
            print("SCOPE done")
        if options.SCOPABW:
            time, results = scopaCollect(SCOPABW, False)
            stats = analysis(results, parseScopaInfo, options.polyType, options)
            fp.write(('{:<7}'*1).format(error))
            printResults(fp, "SCOPABW", [time] + stats);
        if options.CLEAN:
            time, results = seqCleanCollect(CLEAN, False)
            fp.write(('{:<7}'*1).format(error))
            stats = analysis(results, parseCleanInfo, options.polyType, options)
            printResults(fp, "CLEAN", [time] + stats)
        if options.TRIM:
            time, results = seqTrimCollect(TRIM, False)
            stats = analysis(results, parseTrimInfo, options.polyType, options)
            fp.write(('{:<7}'*1).format(error))
            printResults(fp, "TRIM", [time] + stats)    
        if options.POLY:
            time, results = polyCollect(POLY, not options.keep_files)
            stats = analysis(results, parsePolyInfo, options.polyType, options)
            fp.write(('{:<7}'*1).format(error))
            printResults(fp, "POLY", [time] + stats)
        if options.TRIMEST:
            time, results = trimestCollect(TRIMEST, not options.keep_files)
            stats = analysis(results, parseTrimestInfo, options.polyType, options)
            fp.write(('{:<7}'*1).format(error))
            printResults(fp, "TRIMEST", [time] + stats)
        if options.BASICTOOL:
            time, results = basicCollect(BASICTOOL, not options.keep_files)
            stats = analysis(results, parseBasicInfo, options.polyType, options)
            fp.write(('{:<7}'*1).format(error))
            printResults(fp, "BASICTOOL", [time] + stats)
            print("BASIC DONE")

def endLengthSim(options, args):
    """Look at the effect of variation in theo length of the post-fragment tail"""
    options.specificity_test = False
    fp = sys.stdout if options.output == None else open(options.output, "w")
    if not options.Header_off:
        fp.write(('{:<7}'*1).format('ptail') + ('{:<14}'*len(column_names)).format(*column_names) + "\n")

    start = options.post_tail_range[0]
    stop  = options.post_tail_range[1] + 1
    step  = options.post_tail_range[2]

    for tail_len in range(start, stop, step):
        options.post_tail_extend = tail_len
        real_data(*[getattr(options,v) for v in ["sim_file", "outputDir", "real"]], options=options)
        SCOPA, SCOPABW, TRIM, CLEAN, POLY, TRIMEST = launch_jobs(options)

        if options.SCOPA:
            time, results = scopaCollect(SCOPA, False)
            stats = analysis(results, parseScopaInfo, options.polyType, options)
            fp.write(('{:<7}'*1).format(tail_len))
            printResults(fp, "SCOPA", [time] + stats);
        if options.SCOPABW:
            time, results = scopaCollect(SCOPABW, False)
            stats = analysis(results, parseScopaInfo, options.polyType, options)
            fp.write(('{:<7}'*1).format(tail_len))
            printResults(fp, "SCOPABW", [time] + stats);
        if options.CLEAN:
            time, results = seqCleanCollect(CLEAN, False)
            fp.write(('{:<7}'*1).format(tail_len))
            stats = analysis(results, parseCleanInfo, options.polyType, options)
            printResults(fp, "CLEAN", [time] + stats)
        if options.TRIM:
            time, results = seqTrimCollect(TRIM, False)
            stats = analysis(results, parseTrimInfo, options.polyType, options)
            fp.write(('{:<7}'*1).format(tail_len))
            printResults(fp, "TRIM", [time] + stats)    
        if options.POLY:
            time, results = polyCollect(POLY, not options.keep_files)
            stats = analysis(results, parsePolyInfo, options.polyType, options)
            fp.write(('{:<7}'*1).format(tail_len))
            printResults(fp, "POLY", [time] + stats)
        if options.TRIMEST:
            time, results = trimestCollect(TRIMEST, not options.keep_files)
            stats = analysis(results, parseTrimestInfo, options.polyType, options)
            fp.write(('{:<7}'*1).format(tail_len))
            printResults(fp, "TRIMEST", [time] + stats)
        if options.BASICTOOL:
            time, results = basicCollect(TRIMEST, not options.keep_files)
            stats = analysis(results, parseBasicInfo, options.polyType, options)
            printResults(fp, "BASICTOOL", [time] + stats)

def lenTest(options, args):
    """Look at the effect of changing poly(A) tail length and post-fragment tail length"""
    options.specificity_test = False
    fp = sys.stdout if options.output == None else open(options.output, "w")
    if not options.Header_off:
        fp.write(('{:<7}'*2).format('tail', 'p-tail') + ('{:<14}'*len(column_names)).format(*column_names) + "\n")

    tail_start, tail_stop, tail_step = options.tail_range
    post_start, post_stop, post_step = options.post_tail_range

    for i in range(tail_start, tail_stop+1, tail_step):
        options.poly_len = i
        for j in range(post_start, post_stop+1, post_step):
            with open("log.txt", "a") as wp:
                wp.write("ITERATION: %d\t%d\n" % (i,j))
            options.post_tail_extend = j
            
            real_data(*[getattr(options,v) for v in ["sim_file", "outputDir", "real"]], options=options)
            SCOPA, SCOPABW, TRIM, CLEAN, POLY, TRIMEST = launch_jobs(options)


            if options.SCOPA:
                try:
                    time, results = scopaCollect(SCOPA, False)
                    stats = analysis(results, parseScopaInfo, options.polyType, options)
                    fp.write(('{:<7}'*2).format(i,j))
                    printResults(fp, "SCOPA", [time] + stats);
                except Exception as e:
                    open("log.txt", "a").write("ERROR: " + str(e) + "\n")
            if options.SCOPABW:
                try:
                    time, results = scopaCollect(SCOPABW, False)
                    stats = analysis(results, parseScopaInfo, options.polyType, options)
                    fp.write(('{:<7}'*2).format(i,j))
                    printResults(fp, "SCOPABW", [time] + stats);
                except Exception as e:
                    open("log.txt", "a").write("ERROR: " + str(e) + "\n")
            if options.CLEAN:
                try:
                    time, results = seqCleanCollect(CLEAN, False)
                    fp.write(('{:<7}'*2).format(i,j))
                    stats = analysis(results, parseCleanInfo, options.polyType, options)
                    printResults(fp, "CLEAN", [time] + stats)
                except Exception as e:
                    open("log.txt", "a").write("ERROR: " + str(e) + "\n")
            if options.TRIM:
                try:
                    time, results = seqTrimCollect(TRIM, False)
                    stats = analysis(results, parseTrimInfo, options.polyType, options)
                    fp.write(('{:<7}'*2).format(i,j))
                    printResults(fp, "TRIM", [time] + stats)    
                except Exception as e:
                    open("log.txt", "a").write("ERROR: " + str(e) + "\n")
            if options.POLY:
                try:
                    time, results = polyCollect(POLY, not options.keep_files)
                    stats = analysis(results, parsePolyInfo, options.polyType, options)
                    fp.write(('{:<7}'*2).format(i,j))
                    printResults(fp, "POLY", [time] + stats)
                except Exception as e:
                    open("log.txt", "a").write("ERROR: " + str(e) + "\n")
            if options.TRIMEST:
                try:
                    time, results = trimestCollect(TRIMEST, not options.keep_files)
                    stats = analysis(results, parseTrimestInfo, options.polyType, options)
                    fp.write(('{:<7}'*2).format(i,j))
                    printResults(fp, "TRIMEST", [time] + stats)
                except Exception as e:
                    open("log.txt", "a").write("ERROR: " + str(e) + "\n")

def trainingSize(options, args):
    (options.SCOPABW, options.TRIM, options.CLEAN, options.POLY, options.TRIMEST) = [False]*5
    
    output = options.outputDir + "/" + options.sim_file
    input = options.real
    sim_file = options.outputDir + "/" + options.sim_file

    fp = sys.stdout if options.output == None else open(options.output, "w")
    if not options.Header_off:
        fp.write(('{:<10}'*3).format('No BW', 'data_size', 'training') + ('{:<14}'*len(column_names)).format(*column_names) + "\n")

    O = []
    for i in range(int(options.training_size[0])):
        O.extend(read_validation_data(input, options))

    with open(sim_file, "w") as wp:
        IDmap = {}
        for o in O:
            s = o.fastaString()
            id = search("id(\d+)", s).group(1)
            index = IDmap[id]+1 if IDmap.has_key(id) else 1
            IDmap[id] = index
            s = sub("id(%s)" % (id), "id%s.%d" % (id, index), s)
            wp.write(s)
    

    for options.numTrain in range(*options.training_size[1:]):
        for i in range(options.num_trials):
            SCOPA = scopaJob(version = options.version, file_type = options.fileType, inputFile = options.sim_file, e = "", homopolymer_type = options.polyType, d = None, k=0, outputFile = "SCOPA.out.%s" % options.id, DIR=options.outputDir, 
                             terminate=options.terminate, poly=options.poly, numTrain=options.numTrain, no_retrain = False)
            time, results = scopaCollect(SCOPA, not options.keep_files)
            stats = analysis(results, parseScopaInfo, options.polyType, options)
            fp.write(('{:<10}').format('FALSE'))
            fp.write(('{:<10}').format(len(O)))
            fp.write(('{:<10}').format(options.numTrain))
            printResults(fp, "SCOPA", [time] + stats)
            SCOPA = scopaJob(version = options.version, file_type = options.fileType, inputFile = options.sim_file, e = "", homopolymer_type = options.polyType, d = None, k=0, outputFile = "SCOPA.out.%s" % options.id, DIR=options.outputDir, 
                             terminate=options.terminate, poly=options.poly, numTrain=options.numTrain, no_retrain = True)
            time, results = scopaCollect(SCOPA, not options.keep_files)
            stats = analysis(results, parseScopaInfo, options.polyType, options)
            fp.write(('{:<10}').format('TRUE'))
            fp.write(('{:<10}').format(len(O)))
            fp.write(('{:<10}').format(options.numTrain))
            printResults(fp, "SCOPA", [time] + stats)

    
        
        
def timing(options, args):
    # Get / print results in column format
    fp = sys.stdout if options.output == None else open(options.output, "w")
    if not options.Header_off:
        fp.write(('{:<10}').format('size') + ('{:<14}'*len(column_names)).format(*column_names) + "\n")


    for n in range(*options.timing):
        print("n: ", n)
        timing_data(*[getattr(options,v) for v in ["sim_file", "outputDir", "real"]], options=options, n=n)

        m,w,b,x = options.params;

        # Launch simulations
        SCOPA = scopaJob(file_type = options.fileType, inputFile = options.sim_file, e = "", t = options.polyType, w =w, d = None, m = m, b = b, x = x, n = 10000, k=0, outputFile = "SCOPA.out.%d" % options.id, DIR=options.outputDir, terminate=options.terminate, no_retrain=options.no_retrain, front_gap=options.front_gap, poly=options.poly) if options.SCOPA else None
        TRIM = seqTrimJob(options.sim_file, "TRIM.out", DIR = options.outputDir) if options.TRIM else None
        POLY = polyJob(options.sim_file, "POLY", DIR = options.outputDir, terminate=options.terminate) if options.POLY else None
        TRIMEST = trimestJob(options.sim_file, "TRIMEST", DIR = options.outputDir, terminate=options.terminate) if options.TRIMEST else None
        CLEAN = seqCleanJob(options.sim_file, "CLEAN", DIR = options.outputDir, terminate=options.terminate) if options.CLEAN else None

        if options.SCOPA:
            storePBS([SCOPA], open("last.scopa", "w"))
        if options.TRIM:
            storePBS([TRIM], open("last.trim", "w")) 
        if options.CLEAN:
            storePBS([CLEAN], open("last.clean", "w")) 
        if options.POLY:
            storePBS([POLY], open("last.poly", "w"))
        if options.TRIMEST:
            storePBS([TRIMEST], open("last.trimest", "w"))

        if options.SCOPA:
            time, results = scopaCollect(SCOPA, not options.keep_files)
            stats = analysis(results, parseScopaInfo, options.polyType, options)
            fp.write(('{:<10}').format(n))
            printResults(fp, "SCOPA", [time] + stats);
        if options.CLEAN:
            time, results = seqCleanCollect(CLEAN, not options.keep_files)
            stats = analysis(results, parseCleanInfo, options.polyType, options)
            fp.write(('{:<10}').format(n))
            printResults(fp, "CLEAN", [time] + stats)
        if options.TRIM:
            time, results = seqTrimCollect(TRIM, not options.keep_files)
            stats = analysis(results, parseTrimInfo, options.polyType, options)
            fp.write(('{:<10}').format(n))
            printResults(fp, "TRIM", [time] + stats)    
        if options.POLY:
            time, results = polyCollect(POLY, not options.keep_files)
            stats = analysis(results, parsePolyInfo, options.polyType, options)
            fp.write(('{:<10}').format(n))
            printResults(fp, "POLY", [time] + stats)
        if options.TRIMEST:
            time, results = trimestCollect(TRIMEST, not options.keep_files)
            stats = analysis(results, parseTrimestInfo, options.polyType, options)
            fp.write(('{:<10}').format(n))
            printResults(fp, "TRIMEST", [time] + stats)        
                
if __name__ == "__main__":
    with open("log.txt", "w") as fp:
        fp.write(" ".join(sys.argv) + "\n")

    # Parse command-line arguments
    parser = OptionParser()
    simGroup = OptionGroup(parser, "Simulation Switches") 
    parser.add_option_group(simGroup)
    simGroup.add_option("-r", "--real", action="store", type="string", dest="real", help="Use real data from specified file", default="validation/Arabidopsis.Illumina.polyA.500.fa") 
    simGroup.add_option('-R', "--Real", action="store", type="string", dest="Real", nargs = 3, help="Use a file from the validation folder: ORG TECH TYPE", default = None) 
    simGroup.add_option("-e", "--error_rate", action="store", type="float", dest="error_rate", help="base call error rate (default: %default)", default=-1)
    simGroup.add_option("-n", "--num_sequences", action="store", type="int", dest="n", help="Number of simulated sequences (default: %default)", default=20000)
    simGroup.add_option("-s", "--sim_file", action="store", type="string", dest="sim_file", help="Output file (default: %default)", default = "sim.fa")
    simGroup.add_option("-P", "--no_pid", action="store_false", dest="pid", help="Do not use PID extension of file names", default = True)
    simGroup.add_option("-L", "--Liang_data", action = "store", dest = "Liang_data", help="Use Liang data", default = "")
    simGroup.add_option("-m", "--multipler", action = "store", type = "int", dest = "multiplier", help = "Multiply the data set size (e.g. each sequence will generated M simulated sequences)", default = 1)
    simGroup.add_option("--post_tail_extend", action = "store", type = "int", dest = "post_tail_extend", help = "Extend post-tail fragment to a specified length", default = None)
    simGroup.add_option("--poly_len", action = "store", type = "int", dest = "poly_len", help = "Length of the polymer", default = None)
    simGroup.add_option("--num_trials", action = "store", type = "int", dest = "num_trials", help = "Numer of trials in repeated experiments (default: %default)", default = 1)
    simGroup.add_option("--sp_off", "--specificity_test_off", action="store_false", dest = "specificity_test", help = "Turn the specificity test off", default = True)
    #simGroup.add_option("-b", "--basis", action="store", type="string", dest="basis", help="Basis file for simulation (default: %default)", default = "human_annotated.fa")

    toolGroup = OptionGroup(parser, "Tool Switches")
    parser.add_option_group(toolGroup)
    toolGroup.add_option("-t", "--polyType", action="store", type = "string", dest="polyType", help="homopoloymer type (default: $default)", default = "A")
    toolGroup.add_option("-f", "--fileType", action="store", type = "string", dest="fileType", help="file type (sanger, 454, illumina)", default = None)
    toolGroup.add_option("--error_range", action="store", type="float", dest="error_range", nargs=3, metavar = "MIN MAX STEP", help = "Vary simulation error", default = [0, 0.10, 0.01])
    toolGroup.add_option("--tail_range", action="store", type="int", dest="tail_range", nargs=3, metavar = "MIN MAX STEP", help = "Vary tail length", default = [0, 100, 10])
    toolGroup.add_option("--post_tail_range", action="store", type="int", dest="post_tail_range", nargs=3, metavar = "MIN MAX STEP", help = "Vary post-tail fragment length", default = [0, 100, 10])
    toolGroup.add_option("--seed", action="store", dest="seed", help="Set rng seed (default: None)", default = None)
    toolGroup.add_option("--timing", action="store", dest="timing", type = "int", nargs=3, metavar = "start stop skip", help = "Run tools at different sizes inputs", default = None)
    toolGroup.add_option("--training_size", action="store", dest="training_size", type = int, nargs  = 4, metavar = "data_mult start stop skip", help = "Run SCOPE with different training size inputs", default = None)


    selectionGroup = OptionGroup(parser, "Tool Selection")
    parser.add_option_group(selectionGroup)
    selectionGroup.add_option("-S", "--SCOPA_off", action="store_false", dest = "SCOPA", help="Supress SCOPA testing (default: %default)", default = True)
    selectionGroup.add_option("-B", "--SCOPABW_off", action="store_false", dest = "SCOPABW", help="Supress SCOPA (Baum Welch) testing (default: %default)", default = True)
    selectionGroup.add_option("-T", "--TRIM_off", action="store_false", dest = "TRIM", help="Supress SeqTrim testing (default: %default)", default = True)
    selectionGroup.add_option("-C", "--CLEAN_off", action="store_false", dest = "CLEAN", help="Supress SeqClean testing (default: %default)", default = True)
    selectionGroup.add_option("--BT", "--BASICTOOL_off", action = "store_false", dest = "BASICTOOL", help="Suprosee BasicTool testing (default: %default)", default = True)
    selectionGroup.add_option("--TP", "--TRIMPOLY_off", action="store_false", dest = "POLY", help="Supress PolyTrim testing (default: %default)", default = True)
    selectionGroup.add_option("--TE", "--TRIMEST_off", action="store_false", dest = "TRIMEST", help="Supress TrimEST testing (default: %default)", default = True)
    selectionGroup.add_option("--so", "--SCOPE_ONLY", action="store_true", dest = "scope_only", help="Supress every tool but SCOPE", default = False)
    selectionGroup.add_option("--AT", "--all_tools", action ="store", dest = "all_tools", type = int, nargs = 7, help = "Designate which of the six tools should be used (scope++, scope++/bw, SeqTrim, SeqClean, PolyTrim, TrimEST, BasicClean)", default = None)

    experimentGroup = OptionGroup(parser, "Experiment Selection")
    parser.add_option_group(experimentGroup)
    toolGroup.add_option("-E", "--error_variation", action="store_const", const = "error_variation", dest="run_type", help="Vary simulation error rate", default = "basicSim")
    toolGroup.add_option("--PTL", "--post_tail_length", action="store_const", const = "endLengthSim", dest="run_type", help="Effect of variation of the post-tail fragment length", default= "basicSim")
    toolGroup.add_option("--len_test", action="store_const", const="lenTest", dest="run_type", help="Effect of variation of both the post-tail fragment length and the poly(A) length", default= "basicSim")
    toolGroup.add_option("--SB", "--scopeVbasic", action="store_const", const="scopeVbasic", dest = "run_type", help = "Test SCOPE against Basic at various levels of p")

    scopeGroup = OptionGroup(parser, "SCOPE++ parameters")
    parser.add_option_group(scopeGroup)
    scopeGroup.add_option("--params", action="store", dest="params", nargs=3, metavar = "w b x", help = "Set SCOPE params: basic run only", default = [None, None, None])
    scopeGroup.add_option("--no_retrain", action="store_true", dest="no_retrain", help = "Supress SCOPE++ retraining", default = False)
    scopeGroup.add_option("--poly_off", action="store_false", dest="poly", help = "Turn off SCOPE++ poly option", default = True)
    scopeGroup.add_option("-p", "--max_poly_dist", action = "store", type = "int", dest="poly_dist", help="Maximum distance of polyA tail from end (default: %default)", default = 100000)
    scopeGroup.add_option("--minIdentity", action = "store", type = "int", dest="minIdentity", help="minIdentity", default = None)
    scopeGroup.add_option("-v", "--version", action = "store", type = "int", dest = "version", help="SCOPE version to use (default: %default)s", default = 22)
    scopeGroup.add_option("--numTrain", action = "store", type = "int", dest = "numTrain", help="Number of training sequences (default = SCOPE++ default)", default = None)

    outputGroup = OptionGroup(parser, "Output Format Switches")
    parser.add_option_group(outputGroup)
    outputGroup.add_option("--DIRECTORY", action="store", type="string", dest="outputDir", metavar="DIR", help="Directory for results (default: %default)", default = "results_files")
    outputGroup.add_option("-o", "--output", action="store", type="string", dest="output", help="Output file (default: stdout)", default = None)
    outputGroup.add_option("-H", "--Header_off", action="store_true", dest = "Header_off", help="Turn reader row off", default = False)
    outputGroup.add_option("-K", "--keep_files", action="store_true", dest = "keep_files", help="Keep extra files", default = False)
    outputGroup.add_option("--keep_sim", action = "store_true", dest = "keep_sim", help="Retain simulation file", default = False)

    debuggingGroup = OptionGroup(parser, "Debugging Switches")
    outputGroup.add_option("--term", action="store_true", dest="terminate", help="Terminate at first system call (basic sim. only)", default = False)
    outputGroup.add_option("-l", "--last_run", action="store_true", dest="last_run", help="Skip simulation and tool application, using last results; assume pid was surpressed", default = False)
    outputGroup.add_option('--sim_size', action="store", dest="sim_size", type="int", help="Limit number of sequences up in the simulation file", default = -1)
    outputGroup.add_option('--id', action="store", dest="id", type="str", help="ID for files", default = str(os.getpid()))

    parser.add_option_group(debuggingGroup)
    (options, args) = parser.parse_args()

    if options.scope_only:
        options.SCOPABW, options.TRIM, options.CLEAN, options.POLY, options.TRIMEST = [None]*5
        options.SCOPE = True
    elif options.all_tools:
        options.SCOPA, optionsBW, options.TRIM, options.CLEAN, options.POLY, options.TRIMEST, options.BASICTOOL = [bool(x) for x in options.all_tools]

    if options.pid:
        options.sim_file += ".%s" % options.id

    if options.Real:
        options.real = "validation/%s.%s.poly%s.500.fa" % options.Real
        options.polyType = options.Real[2]

    if options.seed:
        seed(options.seed)

    if options.timing:
       options.run_type = "timing" 

    if options.training_size:
        options.run_type = "trainingSize"

    # if options.fileType and not options.fileType in {"sanger", "454", "illumina"}:
    #     sys.stderr.write("Bad file type: " + options.fileType + "\n")
    #     sys.exit(1)

    options.front_gap = options.poly_dist

    eval(options.run_type)(options,args)

    if not options.keep_files:
        for file in glob.glob(options.outputDir + "/*%s*" % (options.id)):
            if file !=  options.outputDir + "/" + options.sim_file or not options.keep_sim:
                os.remove(file)



