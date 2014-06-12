#run!/usr/bin/python2.7
import sys
sys.path.append("../parse/")

from scope_parser import *
from Bio import SeqIO
from redhawk import pbsJobHandler
import re
import os
import subprocess
from itertools import izip
from sys import stderr,exit

#scopa_params = ["z", "s", "f", "r", "d", "e", "t", "w", "m", "n", "k", "b", "x", "file_type", "no_retrain", "front_gap", "poly", "left_gap", "right_gap", "homopolymer_type"]
scopa_params = ["z", "s", "f", "r", "d", "e", "t", "w", "m", "n", "k", "b", "x", "file_type", "no_retrain", "poly", "homopolymer_type", "numTrain", "minIdentity"]

class ScopeError(Exception):
    def __init__(self, value = "Bad SCOPE run"):
        self.value = value

    def __str__(self):
        return str(self.value)

#####################################################################
# Launch a basic_clean.py" job
# Timing disabled
def basicJob(inputFile, outputFile, x = 2, DIR = ".", base_type = 'A', min_length = 10, p = 0.05, format = 'fasta', end = 3, terminate = False):

    #basic_command = "module load python.2.7; python2.7 ../basic/basic_clean.py -x {x} -e {e} -t {t} -m {m} -p {p} -f {format} {inputFile} {outputFile}".format(x=x, e=end, t=base_type, m=min_length, p=p, format=format, inputFile=DIR + "/" + inputFile, outputFile=outputFile)
    basic_command = "python2.7 ../basic/basic_clean.py -x {x} -e {e} -t {t} -m {m} -p {p} -f {format} {inputFile} {outputFile}".format(x=x, e=end, t=base_type, m=min_length, p=p, format=format, inputFile=DIR + "/" + inputFile, outputFile=DIR + "/" + outputFile)
    #timed_command = "/usr/bin/time -f {time_format} -o {cmd}".format(time_format = "%U", cmd = "basic_cmmand")

    if (terminate): # For debugging
        exit(1)
    
    basicObj = subprocess.Popen(basic_command, shell = True, stdin = subprocess.PIPE, stdout = subprocess.PIPE);
    basicObj.output = DIR + "/" + outputFile
    return basicObj
    
def basicCollect(basicObject, clean = True):
    out, err = basicObject.communicate()
    if err: 
        print(err.encode())
        assert False

    return -1, scope_result_generator(SCOPEparser(None, basicObject.output))

def parseBasicInfo(next_seq, id_parser):
    seq_info, seq_len = id_parser(next_seq[0])
    return (next_seq[0], seq_info, [(next_seq[3][-1], int(next_seq[1]), int(next_seq[2]))], seq_len)



#####################################################################
# Launch a SCOPA job and return the pbsJobHandler object
# Input file should be in the DIR directory, and all output files will be written there.
def scopaJob(inputFile, outputFile, z="", s = None, f = None, r = None, d = None, e = None, t = None,
             w = None, m = None, n = None, k = None, b = None, x=None, file_type = None, no_retrain = False, DIR = ".", 
             bfile = "test_scopa.job", front_gap = None, poly = True, terminate = False, left_gap = None, right_gap = None, 
             homopolymer_type = None, numTrain = None, minIdentity = None):


    # Compensate for difference in parameters between version <=14 and version >= 15
    no_retrain = " " if no_retrain else None
    if poly:
        poly = " "
    else:
        poly = None

    #cmd = "/usr/bin/time -f %s -o %s ../src/SCOPE++ -i %s  -o %s" % ("%U", DIR + "/" + outputFile + ".time", DIR + "/" + inputFile, DIR + "/" + outputFile) + " "+ "  ".join(["%s%s %s" % ("-" if len(o)==1 else "--", o, str(eval(o))) for o in scopa_params if not eval(o) == None])


    cmd = "../src/SCOPE++ -i {input_file} -o {output_file}".format(input_file = DIR + "/" + inputFile, output_file = DIR + "/" + outputFile)
    for o in scopa_params:
        if not (eval(o) is None):
            cmd += " %s%s %s" % ("-" if len(o) == 1 else "--", o, str(eval(o)))

    if (terminate):   # For debugging
        print cmd[cmd.find("../SCOPA"):]
        exit(1)

    scopa_job = subprocess.Popen(cmd,shell=True, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    scopa_job.output = DIR + "/" + outputFile + ".tab"
    return scopa_job


# Take the SCPA pbsJobHandeler object, wait for completion, and return the runtime and the results as a list of sequence records
scopaCollect = basicCollect   # Same format, so same function
    
        
# Take the list of sequence records and return a list of results
parseScopaInfo = parseBasicInfo   # Same format, so same function


##########################################
# Launch a SeqTrim job and return the pbsJobHandler object
def seqTrimJob(inputFile, outputFile, DIR = ".", bfile = "test_seqtrim.job", terminate = False):
    #cmd = "/usr/bin/time -o %s -f \"%s\" seqtrim.pl -a qv --seqMinLength 1 -f %s > %s" % (DIR + "/" + outputFile + ".time", "%U", DIR + "/" + inputFile, DIR + "/" + outputFile)
    cmd = "/usr/bin/time -o %s -f \"%s\" ./wrapper.py \"seqtrim.pl -a v --seqMinLength 1 -f %s > %s\"" % (DIR + "/" + outputFile + ".time", "%U", DIR + "/" + inputFile, DIR + "/" + outputFile)
    if (terminate):
        print cmd
        exit(1)
    modules = ["blast", "seqtrim"]
    seqtrim_job = pbsJobHandler(batch_file = bfile, executable = cmd, nodes = 1, ppn = 1, output_location = DIR,
                                RHmodules = modules, walltime = "1:00:00")
    seqtrim_job.name = "TRIM"
    seqtrim_job.output = outputFile
    seqtrim_job.DIR = DIR
    return seqtrim_job.submitjob()

# Take the SeqTrim object, wait for completion, and return the runtime and the results as a list of sequence records
def seqTrimCollect(trimObj, clean = True):
    outputFile = trimObj.DIR + "/" + trimObj.output
    timeFile = outputFile + ".time";

    o,e = trimObj.getResults(clean)
    try:
        time = float(open(timeFile).readline().rstrip("\n"))
    except:
        stderr.write("SeqTrim failed\n")
        exit(1)
    if e.rstrip():
        raise RuntimeError("Error in TRIM execution:\n" + e + "\n");

    return [time, SeqIO.parse(outputFile, "fasta")]

# Take the list of sequence records and return a list of results
def parseTrimInfo(seq_record, id_parser):
    r_id = re.search("id\d+", seq_record.description)
    r_type = re.search("poly(\w)", seq_record.description)
    assert(r_id and r_type)
    seq_info, seq_len = id_parser(seq_record.description)

    #id, actual_start, actual_end, tmp, tmp2, predicted_start = [int(x) for x in re.findall("(\d+)", seq_record.description)][:6]
    #predicted_end = actual_end
    #return [(id, actural_start != actual_end, actual_start, actual_end, predicted_start != predicted_end, predicted_start, predicted_end)]

    r = re.search("Cuts:\s+\[(\d+),(\d+)\]", seq_record.description)
    start = int(r.group(1))
    finish = int(r.group(2))
    return (r_id.group(0), seq_info, [[r_type.group(1)] + ([finish,seq_len] if start == 1 else [0,start-1])], seq_len)


##############################################
# Launch a SeqClean job and return the pbsJobHandler object
def seqCleanJob(inputFile, outputPrefix, nodes = 1, ppn = 1, l = 25, DIR = ".", bfile = "test_seqclean.job", id = 0, terminate = False):
    cmd1 = "cd %s;" % (DIR);
    cmd2 = "cp %s %s.%d;" % (inputFile, inputFile, id)
    cmd3 = "/usr/bin/time -o %s -f \"%s\"" % (outputPrefix + "." + str(id) + ".time", "%U")
    cmd4 = "../wrapper.py \"seqclean %s.%d -o %s.%d.out -r %s.%d.cln -l %d -c %d\"; " % (inputFile, id, outputPrefix, id, outputPrefix, id, l, nodes*ppn);
    cmd5 = "rm %s.%s" % (inputFile, id)
    cmd = " ".join([cmd1, cmd2, cmd3, cmd4, cmd5]);

    
    """
    Jamie Code: subprocess
    """
    ##############
    if (terminate):   # For debugging
       print cmd
       exit(1)
    seqclean_job = subprocess.Popen(cmd,shell=True, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    seqclean_job.wait()
    return [DIR + "/" + inputFile, DIR + "/" +outputFile,DIR + "/" +outputFile+".time"]
    ##############

    """
    Karro Code: redhawk script
    """
    # if terminate:
    #     print cmd
    #     exit(1)
    # modules = ["blast", "seqclean"]
  
    # seqclean_job = pbsJobHandler(batch_file = bfile, executable = cmd, nodes = 1, ppn = 1, output_location = DIR,
    #                              RHmodules = modules, walltime = "01:00:00")
    # seqclean_job.name = "CLEAN"
    # seqclean_job.outputPrefix = outputPrefix + "." + str(id)
    # seqclean_job.DIR = DIR
    # seqclean_job.inputFile = inputFile
    # return seqclean_job.submitjob()



# Take the SeqClean pbsJobHandeler object, wait for completion, and return the runtime and the results as a list of sequence records
def seqCleanCollect(cleanObj, clean = True):
    
    inputFile = scopaObj[0]
    outputFile = scopaObj[1]
    timeFile = scopaObj[2]
    return [0.0, izip(SeqIO.parse(open(inputFile)), "fasta"),open(outputFile)]

    inputFile = cleanObj.DIR + "/" + cleanObj.inputFile
    outputFile = cleanObj.DIR + "/" + cleanObj.outputPrefix + ".cln"
    timeFile = cleanObj.DIR + "/" + cleanObj.outputPrefix + ".time"

    o,e = cleanObj.getResults(clean)
    try:
        time = float(open(timeFile, "r").readline().rstrip("\n"))
    except:
        stderr.write("SeqClean failed\n")
        exit(1)
    #if e.rstrip():
    #    raise RuntimeError("Error in CLEAN execution:\n" + e + "\n");
    return [time, izip(SeqIO.parse(open(inputFile), "fasta"), open(outputFile))]

# Take the list of sequence records and return a list of results
seqCleanRE = re.compile("(id\d+).*\d+(?:\.\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)")
def parseCleanInfo(next_seq, id_parser):
    inputSeq, clnOutput = next_seq
    seq_info, seq_len = id_parser(inputSeq.description)

    r = seqCleanRE.search(clnOutput)
    id, x, y, l = r.group(1), int(r.group(2))-1, int(r.group(3)), int(r.group(4))

    if x == 0 and y ==l:
        return (id, seq_info, [], seq_len)

    elif x > y:
        return (id, seq_info, [[max([(inputSeq.seq.count(c), c) for c in ['A', 'C', 'G', 'T']])[1], 0, l]], seq_len)

    elif x > 0 and y < l:
        return (id, seq_info, [[max([(inputSeq.seq[0:x].count(c), c) for c in ['A', 'C', 'G', 'T']])[1], 0, x],
                               [max([(inputSeq.seq[y:l].count(c), c) for c in ['A', 'C', 'G', 'T']])[1], y, l]], seq_len)

    elif x==0 or y==l:
        s,f = (y,l) if x==0 else (0,x)
        type = max([(inputSeq.seq[s:f].count(c), c) for c in ['A', 'C', 'G', 'T']])[1]
        return (id, seq_info, [[type, s, f]], seq_len)

    else:
        assert(False)
    

##############################################
# Launch a trimpoly job and return the pbsJobHandler object
def polyJob(inputFile, outputPrefix, DIR = ".", bfile = "test_trimpoly.job", id = 0, terminate = False):
    cmd1 = "/usr/bin/time -o %s -f \"%s\"" % (DIR + "/" + outputPrefix + ".time", "%U")
    cmd2 = "./wrapper.py \"cat %s | trimpoly -o %s.out > %s_details.out\"" % (DIR + "/" + inputFile, DIR + "/" + outputPrefix, DIR+ "/" + outputPrefix)
    cmd = " ".join([cmd1, cmd2])

    """
    Jamie Code: subprocess
    """
    if (terminate):
        print cmd
        exit(1)
    seqclean_job = subprocess.Popen(cmd, shell=True, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    seqclean_job.wait()
    return [DIR + "/" + inputFile, DIR + "/" + outputPrefix + "_details.out", DIR + "/" + outputPrefix + ".time"]
    ##################


    """
    Karro COde: redhawk
    """
    if terminate:
        print cmd
        exit(1)
    modules = ["seqclean"]
    
    trimpoly_job = pbsJobHandler(batch_file = bfile, executable = cmd, nodes = 1, ppn = 1, output_location = DIR, RHmodules = modules)
    trimpoly_job.name = "TRIMPOLY"
    trimpoly_job.outputPrefix = outputPrefix
    trimpoly_job.DIR = DIR
    trimpoly_job.inputFile = inputFile
    return trimpoly_job.submitjob()

def polyCollect(trimPolyObj, clean = True):
    inputFile = trimPolyObj[0]
    outputFile = trimPolyObj[1]
    timeFile = trimPolyObj[2]
    return [0.0, izip(SeqIO.parse(open(inputFile), "fasta"), open(outputFile))]

    inputFile = trimPolyObj.DIR + "/" + trimPolyObj.inputFile
    outputFile = trimPolyObj.DIR + "/" + trimPolyObj.outputPrefix + "_details.out"
    timeFile = trimPolyObj.DIR + "/" + trimPolyObj.outputPrefix + ".time"

    o,e = trimPolyObj.getResults(clean)
    try:
        time = float(open(timeFile, "r").readline().rstrip("\n"))
    except:
        stderr.write("TrimPoly failed\n")
        exit(1)
    
    return [time, izip(SeqIO.parse(open(inputFile), "fasta"), open(outputFile))]

def parsePolyInfo(next_seq, id_parser):
    return parseCleanInfo(next_seq, id_parser)

##############################################
# Launch a trimpoly job and return the pbsJobHandler object
def trimestJob(inputFile, outputPrefix, DIR = ".", bfile = "test_trimest.job", id = 0, terminate = False):
    outputPrefix = outputPrefix + "." + str(id)
    cmd1 = "/usr/bin/time -o %s -f \"%s\"" % (DIR + "/" + outputPrefix + ".time", "%U")
    cmd2 = "./wrapper.py \"trimest -noreverse -sequence %s -mismatches 10 -outseq %s\"" % (DIR + "/" + inputFile, DIR + "/" + outputPrefix + ".out")
    cmd = " ".join([cmd1, cmd2])
    if terminate:
        print cmd
        exit(1)
    modules = ["emboss"]


    trimest_job = pbsJobHandler(batch_file = bfile, executable = cmd, nodes = 1, ppn = 1, output_location = DIR, RHmodules = modules)
    trimest_job.name = "TRIMEST"
    trimest_job.outputPrefix = outputPrefix
    trimest_job.DIR = DIR
    trimest_job.inputFile = inputFile
    return trimest_job.submitjob()

def trimestCollect(trimestObj, clean = True):
    inputFile = trimestObj.DIR + "/" + trimestObj.inputFile
    outputFile = trimestObj.DIR + "/" + trimestObj.outputPrefix + ".out"
    timeFile = trimestObj.DIR + "/" + trimestObj.outputPrefix + ".time"

    o,e = trimestObj.getResults(clean)
    try:
        time = float(open(timeFile, "r").readline().rstrip("\n"))
    except:
        stderr.write("trimest failed\n")
        exit(1)

    return [time, izip(SeqIO.parse(open(inputFile), "fasta"), SeqIO.parse(open(outputFile), "fasta"))]

embossReg = re.compile("\[poly-([AT]) tail removed\]")
def parseTrimestInfo(next_seq, id_parser):
    inputSeq, trimSeq = next_seq
    seq_info, seq_len = id_parser(inputSeq.description)

    r = embossReg.search(trimSeq.description)
    if r:
        poly_type = r.group(1)
        start, finish = (len(trimSeq), seq_len) if poly_type == 'A' else (0, seq_len - len(trimSeq))
        return "???", seq_info, [[poly_type, start, finish]], seq_len
    else:
        return "???", seq_info, [], seq_len
        
        
        
