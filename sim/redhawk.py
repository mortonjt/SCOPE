############################################################################
# redhawk.py
# Class for submitting / managing redhawk jobs through a python script
# Written by: John Karro, Jens Muller
# Date: Nov. 1, 2011
#

import re, string, sys
import subprocess
import getopt
import time
import os
import cPickle
from argparse import *
import pwd
import datetime

current_user = pwd.getpwuid(os.getuid())[0]  # When run on redhawk with a nohup, os.getlogin() does not work

redhawkStatsRe = re.compile("\s+C\s+[^\s]+\s*$")
redhawkInQueueRe = re.compile("\s+[R|Q]\s+[^\s]+\s*$")

log_file = "/usr/local/torque/current/var/spool/torque/server_priv/accounting"

class RedhawkError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class pbsJobHandler:
    """A pbsJobHandler corresponds to a job launched (or to be launched) on redhawk.  Once the object is created (and provided with a command-line execution command),
       the user can extract various inforamtion about the job (current status, output, etc...) and cleanup files."""
    def __init__(self, batch_file, executable, use_pid = True, job_name = None, nodes = 1, ppn = 1, walltime = "40:00:00", address = None, join = False, env = None, queue = None, mail = None, output_location = None, chdir = None, RHmodules = None, file_limit = 6, file_delay = 5, timing = False): 
        """Constructor.  Requires a file name for the batch file, and the execution command.  Optional parmeters include:
           * use_pid: will embded a process id into the batch file name if true.  Default = true.
           * job_name: A name for the redhawk name.  Default = the batch file name.
           * nodes: number of nodes required for the job.   Default = 1.
           * ppn: number of processors needed for the job.  Default = 1.
           * walltime: Maximum allowed runtime for the job (hours:minutes:seconds).  Default = 40:00:00.  Max. allowed: 400:00:00.
           * mail = when to send email.  Any combination of:
             b   send mail when job begins
             e   send mail when job ends
             a   send mail when job aborts
           * address: additional email addresses to send notification (comma seperated)
           * join: If true, the stdout and stderr files are combined into one file
           * queue: redhawk queue to run on.  Default: redhawk chooses.
           * output_location: Directory to place output files.  Default: current directory.
           * RHmodules: A list of redhawk modules to be loaded before run (e.g. ['Blast+']).  Default: none.
           * Timing: Will track the time.  If command is a sequence, only works on the first command (up to the first semi-colon)."""
        self.batch_file_name = batch_file
        if use_pid:
            self.batch_file_name = self.batch_file_name + "." + str(os.getpid())
        self.executable_name = executable
        self.jobname = job_name if job_name else batch_file
        self.join = 'n'
        self.file_limit = file_limit
        self.file_delay = file_delay
        self.status = "unstarted"

        f = open(self.batch_file_name, 'w')

        f.write("#!/bin/bash -l\n")
        s="#PBS -N "+ self.jobname +"\n"
        f.write(s)

     
        #some defaults:
        self.nodes = nodes
        self.ppn = ppn
        self.walltime = walltime
        self.modules = RHmodules
        self.output_location = output_location if output_location else "."

        s="#PBS -l nodes="+ str(self.nodes)+":ppn="+str(self.ppn)+"\n"
        f.write(s)
        s="#PBS -l walltime="+self.walltime+"\n"
        f.write(s)
        
        if join:
            f.write("#PBS -j oe\n")
            
        if address:    
            s="#PBS -M "+address+"\n"
            f.write(s)

        if queue:
            s="#PBS -q "+queue+"\n"
            f.write(s)

        if env:    
            f.write("#PBS -V\n")

        if mail:
            s="#PBS -m "+mail+"\n"
            f.write(s)

        if output_location:
            s="#PBS -o "+output_location+"\n"
            f.write(s)
            s="#PBS -e "+output_location+"\n"
            f.write(s)

        if chdir:
            s="cd "+chdir+"\n"
            f.write(s)

        else:
            s="cd $PBS_O_WORKDIR\n";
            f.write(s);

        if timing:
            self.timing = self.output_location + "/" + self.batch_file_name + ".timing"
            self.executable_name = "/usr/bin/time -o " + self.timing + " -f \"%s\" " + self.executable_name
        else:
            self.timing = None

        if self.modules != None:
            self.executable_name = "; ".join(["module load " + x for x in self.modules]) + "; " + self.executable_name

        f.write(self.executable_name)

        f.close()
        self.jobid=0;

### submitjob
### Parameters:
###   file is the job script file
###   if preserve is True, don't delete the job script. Delete otherwise
###   optionalFlag is the flag after qsub
###   retry (default set to retry 600 times), the number of times, the job will be submitted in retry
###   seconds between retry (default is 10 seconds)
### return job id if successful
### return -1 if not
    def submit(self, preserve=True, print_qsub = False, job_limit = 200, delay=10, user=current_user ):
        """Submit job to redhawk.  Optional parameters:
           * preserve: if False, delete the batch file.  Default = true.
           * job_limit: If the user currently has this many jobs on the batch, wait until one finishes.
           """
        if job_limit > 0:
            limit_jobs(limit=job_limit, delay=delay, user=user)

        optionalFlag=''
        retry=600
        sleepTimeBetweenRetry=10
        trial=1;
        cmd = "qsub " + optionalFlag + " " + self.batch_file_name

        while (trial < retry): 
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (output, error) = p.communicate()
            if p.returncode == 0:
                break
            trial = trial + 1

        if trial == retry:
            return -1
            
        if not preserve:
            os.remove(self.batch_file_name)
                
        t=re.split('\.',output)
        self.jobid=t[0]
        self.ofile = self.output_location + "/" + self.jobname + ".o" + str(self.jobid)
        self.efile = self.output_location + "/" + self.jobname + ".e" + str(self.jobid)

        if print_qsub:
            print 'qsub jobid', self.jobid
                
        self.status = "running"
        return self

    def submitjob(self, preserve=False, print_qsub = False, job_limit = 200, delay=10, user=current_user ):
        """Depricated: replaced with submit()"""
        return self.submit(preserve, print_qsub, job_limit, delay, user)

### isJobRunning
### This is primarily useful for waiting a _submitted_ job to finish
###    return False if the job is done, completed for sure
###    return True if the job is in Q, R states [ or that PBS/Torque is not available ]
### Prereq is jobid must be a submitted job
    def isJobRunning ( self, numTrials = 3, delay = 5 ):
        """Query of the object represented by the job is still running."""
        #cmd = "qstat " + str(self.jobid)
        
        #magicString='Unknown Job Id'  ### magicString _might_ need to be changed if Torque version changes
        #(output, error) = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()


        #for i in range(numTrials):
        counter = 1
        while True:
            file_exists = False
            if self.ofile_exists():  #output.find(magicString) >=0 or redhawkStatsRe.search(output):
                self.status = "finished"
                return False

            cmd = "qstat " + str(self.jobid)
            magicString = 'Unknown Job ID'
            (output, error) = subprocess.Popen(cmd.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
            if redhawkInQueueRe.search(output):
                return True

            time.sleep(delay)
            #if counter % 10 == 0:
            #   print ("isJobRunning: %d, %s, %s, %s" % (counter, self.ofile, os.getcwd(), str(os.path.isfile(self.ofile))))
            counter += 1

        raise RedhawkError("RedHawk error: out of queue, no output file.  OFILE: %s" % (self.ofile))

    
    def wait(self, delay=10):
        """Spin until job completes."""
        while self.isJobRunning() == True:
            time.sleep(delay)
        return self.ofile_exists()  

    def wait_on_job(self, delay=10):  
        """Depricated: replace with wait"""
        return self.wait(delay)

    def ofile_name(self):
        """Get the name of the file containing the job stdio output."""
        return self.ofile

    def efile_name(self):
        """Get the name of the file containing the job stderr output."""
        return self.efile

    def timing_name(self):
        """Get the name of the timing file."""
        return self.timing

    def ofile_exists(self):
        """Does the file contiining the job stdout output exist?"""
        return os.path.isfile(self.ofile)

    def efile_exists(self):
        """Does the file contiining the job stderr output exist?"""
        return os.path.isfile(self.efile)


    def ofile_handle(self):
        """Return a handle to the file containing the job stdout output."""
        if not self.status == "finished":
            raise NameError("redhawk: unfinished ofile check")
        tries = 0
        while not self.ofile_exists() and tries < self.file_limit:
            time.sleep(self.file_delay)
            tries = tries+1
        
        if os.path.isfile(self.ofile_name()):
            return open(self.ofile_name(), "r")

        raise NameError("redhawk: unfound ofile")

    def efile_handle(self):
        """Return a handle to the file containing the job stderr output."""
        if not self.status == "finished":
            raise NameError("redhawk: unfinished efile check")

        tries = 0
        while not self.efile_exists() and tries < self.file_limit:
            time.sleep(self.file_delay)
            tries = tries+1
        
        if os.path.isfile(self.efile_name()):
            return open(self.efile_name(), "r")

        raise NameError("redhawk: unfinished efile check")


    def ofile_string(self):
        """Return the entire contents of the stdout file as a single string."""
        fp = self.ofile_handle()
        if (fp):
            return "\n".join([line.rstrip() for line in fp]) + '\n'
        return None

    def efile_string(self):
        """Return the entire contents of the stderr file as a single string."""
        fp = self.efile_handle()
        if (fp):
            return "\n".join([line.rstrip() for line in fp]) + '\n'
        return None

    def get_timing(self, delete_timing = True):
        """Get the time-generated user runtime and delete the timing file.
        (Assumed to be the last line of the timing file.)
        By default, erases the timing file."""
        if not self.timing:
            return None
        else:
            try:
                with open(self.timing) as fp:
                    lines = fp.readlines()
                if delete_timing:
                    os.remove(self.timing)
                return float(lines[-1])
            except:
                sys.stderr.write("Redhawk.runtime: invalid runtime (%s)\n" % (self.timing))
                sys.exit(1)
                return None
                                
    def get_memory(self, num_days  = 31, days_ago = 0):
        """Get memory usage of process.
        Looks up the information by date.  Will start looking days_ago days, and
        check the num_days previous consecutive days for the process."""

        s = str(self.jobid) + ".*" + "mem=(\d+\wb)"
        reg = re.compile(s)
        for d in range(days_ago, num_days+1):
            t = datetime.datetime.now() - datetime.timedelta(days = d)
            accounting_file = log_file + "/" + t.strftime("%Y%m%d")
            try:
                with open(accounting_file) as fp:
                    line = fp.readline()
                    while line:
                        r = reg.search(line)
                        if r:
                            break
                        line = fp.readline()
                    if r:
                        return r.group(1)
            except IOError:
                pass
        return None
                


    def erase_files(self):
        """Erase the stdio and stderr files."""
        self.ofile_handle()
        self.efile_handle()

        os.remove(self.ofile_name())
        os.remove(self.efile_name())

        return None
    
    def get_results(self, cleanup=True):
        """Retrieve strings and cleanup.  Return truple will include runtime
        if it was recorded and not suppressed."""
        self.wait()
        stdout_str = self.ofile_string()
        stderr_str = self.efile_string()
        if cleanup:
            self.erase_files()
                
        return stdout_str, stderr_str

    def getResults(self, cleanup=True):
	"""Legacy"""
	return self.get_results(cleanup)
	

    def getMemory(self):
        """Return the amount of memory used by process"""

    def wait_on_job_limit(self, limit=200, delay=10, user=current_user):
        """Depricated: use the stand-along function job_limit."""
        limit_jobs(self, limit, delay, user)

### Hold until the user has < limit jobs in the circulation
def limit_jobs(limit=200,delay=10,user=current_user):
    """Spin until the user has less < limit jobs in circulation.
        limit = 0 signals no limit."""
    if limit == 0:
        return None
    while 1==1:
        numJobsInQueue = get_number_of_jobs_in_queue(current_user)

        if (numJobsInQueue < limit):
            return None
        time.sleep(delay)


### return the number of jobs in queue, whatever the state is
def get_number_of_jobs_in_queue(user=current_user):
    """Get the number of user jobs currently sitting in the queu or running."""
    cmd = "qstat -u "+user + " 2>/dev/null | grep " + user

    output,error = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

    return len([True for line in output.split("\n") if line and not redhawkStatsRe.search(line)])
                                                        
         

def storePBS(pbsList, fp):
    """Pickle and save a PBS job to a specified file pointer."""
    for p in pbsList:
        cPickle.dump(pbsList, fp)

def loadPBS(fp):
    """Recover a list of pickled jobs from a specified file pointer."""
    return [obj for obj in cPickle.load(fp)]


def relaunch(args = sys.argv, force = False, walltime = "40:00:00", python = "python"):
    """Detects whether program is being run on the head node.  If so, relaunch identical program on a compute node and quit."""
    o, e = subprocess.Popen(["hostname"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if force or re.match("mualhpcp", o):
        o = pbsJobHandler(batch_file = "relaunch" + str(os.getpid()),  executable = python + " " + " ".join(args), walltime=walltime)
        o, e = o.submit().getResults()
        sys.stdout.write("STDOUT: " + o)
        sys.stderr.write("STDERR: " + e)
        return True
    return False
        
############## Allow for a direct launch of a program
if __name__ == "__main__":
    parser = ArgumentParser(description='Launch a redhawk job.')

    parser.add_argument('command', action = "store", type = str, nargs = '+')


    term = parser.add_argument_group("Input/Output switches")
    term.add_argument('--wait', action = "store_true", dest = "wait", help="wait on job completion", default = False)
    term.add_argument('-o', '--output', action = "store", dest = "output", help = "output file (stdout by default)", default = None)
    term.add_argument('-e', '--error', action = "store", dest = "error", help = "error file (stderr by default)", default = None)
    term.add_argument('-r', '--resources', action = "store", dest = "resources", help = "resource file (same as output by default)", default = None)
    term.add_argument('-S', '--suppress_output', action = "store_true", dest = "suppress", help="Quit without output", default = False)

    settings = parser.add_argument_group("Main job-related settings")
    settings.add_argument('-n', '--nodes', action = "store", type = int, dest = "nodes", help = "Number of nodes", default = 1)
    settings.add_argument('-p', '--ppn', action = "store", type = int, dest = "ppn", help = "Number of processors per node", default = 1)
    settings.add_argument('-w', '--walltime', action = "store", type = str, dest = "walltime", help = "Reserved walltime", default = "10:00:00")
    settings.add_argument('-m', '--modules', action = "store", type = str, nargs = "+", dest = "RHmodules", help = "required redhawk modules", default = None)
    settings.add_argument('-O', '--output_location', action = "store", type = str, dest = "output_location", help = "Output location", default = None)
    settings.add_argument('-d', '--dir', action = "store", type = str, dest = "target_directory", help = "target directory", default = None)


    settings2 = parser.add_argument_group("Less important job-related settings")
    settings2.add_argument('-b', '--batch', action = "store", type = str, dest = "batch", help="Batch file name", default = "redhawk")
    settings2.add_argument('-P', '--pid_off', action = "store_false", dest = "pid", help = "Surpress use of pid in file names", default = True)
    settings2.add_argument('-T', '--time_off', action = "store_false", dest = "timing", help = "Suppress runtime reporting", default = True)
    settings2.add_argument('-M', '--memory_off', action = "store_false", dest = "memory", help = "Suppress memory usage reporting", default = True)
    settings2.add_argument('-K', '--keep_files', action = "store_true", dest = "keep", help = "Keep files generated", default = False)




    args = parser.parse_args()
    o = pbsJobHandler(batch_file = args.batch, executable = " ".join(args.command), use_pid = args.pid, nodes = args.nodes, ppn = args.ppn, timing = args.timing,
                      walltime = args.walltime, output_location = args.output_location, chdir = args.target_directory, RHmodules = args.RHmodules).submit(preserve = args.keep)

    if args.suppress:
        sys.exit(0)

    o.wait()
    ofp = sys.stdout if not args.output or args.output == '-' else open(args.output, "w")
    efp = sys.stderr if not args.error or args.error == '-' else open(args.error, "w")
    rfp = ofp if not args.resources else (sys.stdout if args.resources == '-' else open(args.resources, "w"))
    out, err = o.get_results(cleanup = not args.keep)
    ofp.write(out)
    ofp.write(err)

    t = o.get_timing(delete_timing = not args.keep)
    if args.timing:
        rfp.write("Time: " + str(t) + "\n")
    if args.memory:
        rfp.write("Memory: " + o.get_memory() + "\n")
        



################################
### Sample code: Running a single job
### exe = "ls *"     # The command we want to run (to run multiple commands, seperate with semi-colons)
### o = pbsJobHandler(batch_file = "batch.txt", executable = exe, mail = "bea");    # Set up the job
### o.submit()                           # Submit the job to a redhawk queue
### if o.isJobRunning(): print "yes"     # Check to see if job is still in the queue or running
### o.wait()                             # "spin" until job is finished
### output_string = o.ofile_string()     # Get the ouput
### o.erase_files()                      # Erase the output files




