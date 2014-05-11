#!/usr/bin/python2.7

from redhawk import *
import argparse
import subprocess
import re

def run_test(error_rate, version, params, file, type, output_file, args, seed, term = False):   
    r = re.match("validation/(\w+)\.", file)
    cmd = "./karro_sim.py --so -K -E " + " ".join(["-v " + str(version),
                                                   "-m " + str(args.multiplier),
                                                   "--no_retrain" if args.no_retrain else "",
                                                   "--id " + r.group(1) + "." + type,
                                                   "--seed " + str(seed),
                                                   "--error_range %5.3f 1 2" % (error_rate),
                                                   "--params " + " ".join([str(x) for x in params]),
                                                   "-r " + file,
                                                   "-t " + type,
                                                   "-o " + output_file,
                                                   "--term" if term else ""])
                                               

    output, error = subprocess.Popen(cmd.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
        raise ValueError(error)

    if term:
        return output
    #print output
    fp = open(output_file)
    line = fp.readline()
    D = {key:i for i,key in enumerate(re.split("\s+", line.rstrip()))}
    arr = re.split("\s+", fp.readline().rstrip())

    return "Sensitivity: %5.2f\tSpecifity: %5.2f\t Percent Correct: %5.2f" % (float(arr[D["sens"]]), float(arr[D["spec"]]), float(arr[D["pct_crt"]]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Validate a run of SCOPE on semi-simulated data derived from a selection of fixed files")
    parser.add_argument('-e', '--error_rate', action = 'store', type = float, dest = 'error_rate', help = 'rate of base-read error (default = %(default)s)', default = 0.05)
    parser.add_argument('-v', '--version', action = 'store', type = int, dest = 'version', help = 'SCOPE version number (default = %(default)s)', default = 16)
    parser.add_argument("--params", action="store", dest="params", nargs=4, metavar = ('m', 'w', 'b', 'x'), help = "Set SCOPE params (default = %(default)s)", default = [7, 12, 1, 1])
    parser.add_argument("-s", "--seed", action='store', type = int, dest = "seed", help = "Set the RNG seed (default = %(default)s)", default = 1)
    parser.add_argument("--term", action="store_true", help = "Print only the SCOPE commands to be executed, without actually executing them", default = False)
    parser.add_argument("--no_retrain", action="store_true", dest="no_retrain", help = "Supress SCOPE++ retraining", default = False)
    parser.add_argument("-m", "--multiplier", action = "store", type = int, dest = "multiplier", default = 1)
    args = parser.parse_args()

    
    files = [("validation/Chlamy.454.polyA.500.fa", 'A'), ("validation/Chlamy.454.polyT.500.fa", 'T'),
             ("validation/Arabidopsis.Illumina.polyA.500.fa", 'A'), ("validation/Arabidopsis.Illumina.polyT.500.fa", 'T'),
             ("validation/Human.Sanger.polyA.500.fa", 'A'), ("validation/Human.Sanger.polyT.500.fa", 'T')] 
    for f,t in files:
        print("%s: %s" % (f, run_test(error_rate = args.error_rate, version = args.version, params = args.params, file = f, type = t, output_file = "output.tmp", args = args, term = args.term, seed=args.seed)))



