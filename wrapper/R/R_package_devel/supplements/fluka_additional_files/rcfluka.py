#! /usr/bin/env python
#########################################################################
# Copyright (C) 2009 Niels Bassler, bassler@phys.au.dk
#
# This program is free software; you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the 
# Free Software Foundation; either version 3 of the License, or (at your 
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along 
# with this program; if not, see <http://www.gnu.org/licenses/>.
#
# rcfluka.py
#
# Script which eased the submission of a fluka job to a
# condor cluster.
#
# NB: Niels Bassler, bassler@phys.au.dk,
# RH: Rochus Herrmann, rhe@phys.au.dk
# GPL License applies.
#
# This script submits a FLUKA input file to a condor
# cluster. The amount of nodes to use are specified
# with the -M x option.
# If other executables than the default rfluka is needed,
# then these are relinked on each node sperately,
# to avoid problems with heterogenous libraries across
# the condor pool.
# Several unserroutines can be specified as well, and
# these will be submitted to the individual nodes.
# After the calculations have completed, the files are
# all transferred back to the condor client where the
# job were submitted from, and renamed as if the job
# was submitted with the rfluka. The output can
# therfore directly be used by the uswxxx routines
# afterwards.
# Note that the statistics is increased according to the
# amount of nodes required.
#
# Change log:
#   - 24.09.2009, RH changed sequence, first send mail, than cleanup
#   - 26.08.2009, NB Version 1.45, fix -N option
#   - 03.08.2009, NB Version 1.44, single job submissions (-x)
#   - 10.06.2009, RH Version 1.43, option to send email after
#                    job is finished
#   - 27.02.2009, NB Version 1.42, fix executable bit of .sh
#   - 09.02.2009, NB Version 1.41, fix flush and timout
#   - 06.02.2009, NB Version 1.4, ID-hashing and inifinite loop fix.
#   - 25.12.2008, NB Version 1.35, options -n added
#   - 22.11.2008, NB options -i and -o added.
#   - 20.11.2008, NB individualize stdout from condor submissions.
#                 Reinserted sleep().
#   - 18.11.2008, NB wait some more time, before copying files
#   - 13.11.2008, NB increased max number of submits to 9999
#                 and removed sleep()
#   - 12.11.2008, NB Version 1.3, added -f option
#                 .inp suffix is now optional
#   - 05.09.2008, NB Version 1.2, added testing, cleanup and determinism 
#   - 04.09.2008, RH little correction of randomization 
#   - 01.09.2008, RH Version 1.1, for every run a new random
#                 seed is generated 
#   - 31.07.2008, NB First release 1.0
#
# TODO:
#   - Consequent checking existence of files etc.
#   - better  complete clean "-C"
#   - standard universe
#   - test timeouts
#   - reduce long temporary filenames

import os, re, sys
import shutil # for fast file copy
import time, random, glob
from subprocess import Popen
from optparse import OptionParser

version = "1.45"

dir_rfluka = '$FLUPRO/flutil/rfluka'



def condor_ids():  # returns list of running processes
    active_ids = {}
    pipe = os.popen("condor_q -format \"%f\\n\" clusterid",'r')
    output = pipe.readlines()
    pipe.close()
    for lines in output:
        a_id = float(lines)
        active_ids[a_id] = 1
    # this is a list of active keys now:
    return(active_ids.keys())

def cleanup():
    # cleanup, improve this
    print sys.argv[0] + ": ---- Cleanup  -------------------------------"
    for filename in glob.glob('_rcfluka_*') :
        os.remove( filename ) 



nr_jobs = 1

seed = 0.0    # init this value for deterministic runs
# parse options:
args = sys.argv[1:]
if args.__len__() < 1:
    print ""
    print "rcfluka.py Version", version
    print ""
    print "Error: no input file stated."
    print ""
    print "Usage example: run c12_400.inp on 5 nodes with user routines:"
    print "rcfluka.py -M5 c12_400 -s fluscw.f,comscw.f -l ldpm3qmd"
    print ""
    print "Type rcfluka.py -h for complete option list."
    print ""
    sys.exit(0)

parser = OptionParser()
parser.add_option("-f", "--file", dest="FILE",
                  help="additional optional files needed by user routines", metavar="FILE")
parser.add_option("-s", "--source", dest="SOURCE",
                  help="Filename(s) to be compiled and linked in.", metavar="FILE")
parser.add_option("-e", "--executable", dest="EXE",
                  help="Set name of executeable generated (optional).", metavar="string")
parser.add_option("-M", "--MACHINES", dest="max",
                  default=1,help="Number of jobs to submit.", metavar="int")
parser.add_option("-N", "--NUM", dest="num",
                  help="Start at job number N.", metavar="int")
parser.add_option("-n", "--nice", action="store_true", dest="nice",
                  default=False, help="Submit as nice job.")
parser.add_option("-d", "--deterministic", action="store_true", dest="deter",
                  default=False, help="Deterministic run, seeds are not randomzed.")
parser.add_option("-c", "--clean", action="store_true", dest="clean",
                  default=False, help="Remove temporary files. Use carefully.")
parser.add_option("-C", "--bigclean", action="store_true", dest="bigclean",
                  default=False, help="Remove almost all files, including data. Dangerous!")
parser.add_option("-t", "--test", action="store_true", dest="test",
                  default=False, help="Test, do not submit to condor.")
parser.add_option("-l", "--linker", dest="LINKER",
                  help="Linker to use, e.g. ldpm3qmd.", metavar="FILE")
parser.add_option("-o", "--optioncondor", dest="OPTC",
                  help="Additional optional string to be included in condor submit file.", metavar="string")
parser.add_option("-i", "--includecondor", dest="INCC",
                  help="Include contents of FILE in condor submit file.", metavar="FILE")
parser.add_option("-m", "--mail", dest="mail", 
		  help="Sends report to given e-mail address when job is finished", metavar="string")
parser.add_option("-x", "--single", action="store_true", dest="SINGLE",
                  default=False, help="Submit only once and exit after submission. This overrides -MAX and --mail option.")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()


if options.bigclean == True:
    print "Removing all generated files."
    for filename in glob.glob('_rcfluka_*') :
        os.remove( filename )
    for filename in glob.glob('*_fort*') :
        os.remove( filename )
    for filename in glob.glob('*_node*') :
        os.remove( filename )
    for filename in glob.glob('*ranc*') :
        os.remove( filename )
    for filename in glob.glob('*ranc*') :
        os.remove( filename )
    for filename in glob.glob('*log') :
        os.remove( filename )
    for filename in glob.glob('*out') :
        os.remove( filename )
    for filename in glob.glob('*err') :
        os.remove( filename )
    for filename in glob.glob('*_usrbin*') :
        os.remove( filename )
    sys.exit()
    

# hash table for condor_submissions:
sub_list = {}

nr_jobs = int(options.max)

start_job = 0
if options.num != None:
    start_job = int(options.num)

full_input_filename = args[0]
input_filename = os.path.splitext(os.path.basename(full_input_filename))[0]

#print "Inputfilename = "
#print input_filename
#print type(nr_jobs), nr_jobs
print "Preparing ", nr_jobs, " condor submissions..."

if os.path.isfile(input_filename+".inp") is False:
    print sys.argv[0] +" Error: ---- Could not find file "+ input_filename + ".inp Exiting."
    sys.exit()

# increase if you want. This is just for protection of typing errors etc.
if nr_jobs > 9999:
    print sys.argv[0] +" Error: ---- No more than 9999 submissions supported. Exiting."
    sys.exit()

# generate new file names for each node.

if options.INCC != None:
    if os.path.isfile(options.INCC) is False:
        print sys.argv[0] +" Error: ---- Could not find file "+ options.INCC + "Exiting."
    	sys.exit()
    # read files
    inc_condor = open(options.INCC, 'r').readlines()

if options.SINGLE == True:
    nr_jobs = 1

for ii in range(nr_jobs) :
    i = ii + start_job
    temp_filename = input_filename + "_node%04d" % i
    temp_filename_inp = temp_filename + ".inp"
#    print temp_filename
    # copy original input file to new node
    # TODO: check if input file exists
    if options.SINGLE == False:
        temp_filename = input_filename + "_node%04d" % i
	temp_filename_inp = temp_filename + ".inp"
	shutil.copyfile(input_filename+".inp", temp_filename_inp)
    else:
        temp_filename = input_filename
	temp_filename_inp = temp_filename + ".inp"
	

    # random seeds for input-files
    if options.deter == False:
        seed = float(random.randint(1,9E+5))
    else:
        seed += 1
    flkinp = open( temp_filename_inp , 'r')
    content = flkinp.readlines()
    flkinp.close()
    newflkinp = open( temp_filename_inp , 'w+')
    j = 0
    while j < len(content):
        if re.match('RANDOMIZ',content[j]) is not None:
            newflkinp.write('RANDOMIZ         1.0' + str( seed ).rjust(10)+'\n')
        else:
            newflkinp.write(content[j])
        j = j + 1
    newflkinp.close()


    shell_filename = "_rcfluka_"+temp_filename+".sh"

    ##################################### build shell script
    file = open(shell_filename, "w")
    file.write( '#!/bin/bash -i\n')

    # if we need to compile
    if options.SOURCE != None:
        if options.LINKER == None:
            print sys.argv[0] +" Error: Please specify linker to use with the -l option. Exiting."
            sys.exit()
    
    if options.LINKER != None:        
        file.write( '$FLUPRO/flutil/'+options.LINKER)
        if options.SOURCE != None:
                for fff in options.SOURCE.split(","):
                    file.write( " " + fff)
        if options.EXE == None:
            options.EXE = "./_rcfluka_generetated_flukaexec"
        file.write(" -o " + options.EXE + '\n') 
        # after compilation, execute FLUKA with these options:
        file.write( '$FLUPRO/flutil/rfluka -N0 -M1 -e ' + options.EXE +" " + temp_filename + "\n")
    else:
        file.write( '$FLUPRO/flutil/rfluka -N0 -M1 ' + temp_filename + "\n")
    file.close()    
    # something is bad with the ownerships of the files when transferring. this is a fix:
    os.chmod(shell_filename,0744)

    condor_universe = "vanilla"
    condor_filename = "_rcfluka_condor_submit"+temp_filename+".cdr"
    condor_stdout = "_rcfluka_condor_stdout_"+temp_filename+".txt"
    condor_log = "_rcfluka_condor.log"
    

    #######################  generate condor script
    #    if nr_jobs == 0:
    file = open(condor_filename, "w")
    condor_string =  "##################################################\n"
    condor_string += "# rcfluka.py autogenerated condor job description.\n"    
    condor_string += "Executable      = "  + shell_filename + "\n"
    condor_string += "Universe        = " + condor_universe + "\n"
    condor_string += "Output          = " + condor_stdout + "\n"
    # just on abacus 
    condor_string += "Requirements = Arch == \"X86_64\"          \n"
    #
    if options.OPTC != None:
        condor_string += options.OPTC + "\n"
    # nice job?
    if options.nice == True:
	condor_string += "nice_user = True\n"

    # add input condor file bla bla here
    if options.INCC != None:
        for lines in inc_condor:
            condor_string += lines

    # rh
    # condor just notifies if a job experiences an error, not at every
    # successfully completed job
    condor_string += "Notification = Error \n"
    # hr

   
    condor_string += "Log             = " + condor_log + "\n"
    condor_string += "should_transfer_files  = YES\n"
    condor_string += "when_to_transfer_output = ON_EXIT\n"
    condor_string += "transfer_input_files = " + temp_filename + ".inp"
    if options.SOURCE != None:
        for fff in options.SOURCE.split(","):
            condor_string += ", " + fff
    if options.FILE != None:
        for fff in options.FILE.split(","):
            condor_string += ", " + fff
    condor_string +="\n"
    condor_string += "Queue\n"
##     else:
##         condor_string += "Executable = " + shell_filename + "\n"
##         condor_string += "transfer_input_files = " temp_filename + ".inp"
##         for fff in options.SOURCE-split(","):
##             condor_string += " " + fff
##             condor_string +="\n"
##             condor_string += "Queue\n"
    file.write(condor_string)
    file.close()
     
    # execute next line with popen
    exec_string = "condor_submit "+condor_filename+ "\n"
#    print exec_string

    if options.test==False :
        # submission process is unparallelized
        sss = os.popen(exec_string)
        # grab stdout

        sss_lines = sss.readlines()
        sss.close() # is this correct?
        sss_id = float(sss_lines[-1].split()[5]) # needs float becasause of .        
        print "Submitted ", condor_filename, "with Condor id: ", sss_id
        sys.stdout.flush() # flush output for nohup
        sub_list[sss_id] = 1 # 1 for running
        

    #    p.append(Popen(exec_string, shell=True))
    # aparantly there is some problem, if submitted to fast?
    #    time.sleep(1)

    # wait for PIDs to time out and remove old .inp files

if options.SINGLE == True:
    #cleanup no good idea, since it will remove open files which is still used by condor.
    #if options.clean == True:
    #    time.sleep(5)
    #    cleanup()
    print " ---- Single condor job is now running. Exit. ----- "
    sys.exit()

print "-------------- Now waiting for completed jobs ------------------"

if options.test == False:
    finito = False
    timeout = 0
    # spaghetthi time!
    while finito is False:
        timeout += 1
        time.sleep(15) # check every 15 seconds.
        # update sub_list with processes started by this invokation        
        # mark all as done
        finito = True
        for b in sub_list.keys():
            if sub_list[b] == 1:
                sub_list[b] = 2
            else:
                sub_list[b] = 0
        # get current list of active processes:    
        id_list = condor_ids()
        # and remark list
        for b in id_list:
            if (b in sub_list):
                sub_list[b] = 1
                finito = False
        # now all keys with 2 recently completed
        for b in sub_list.keys():
            if sub_list[b] == 2:
                print "rcfluka submission ",b," terminated."
                sys.stdout.flush()
                sub_list[b] = 0

        if timeout > 200000:
            print "rcfluka TIMEOUT Error!"
            print "The submitted job took more than 34 days to complete."
            finito = True

        # or stop if no condor jobs are running
        if id_list == []:
            finito = True

# then start looping over all data stuff


for ii in range(nr_jobs):
    i = ii + start_job
    # old version, get pids
    # os.waitpid(p[i].pid,0)
    temp_filename     = input_filename + "_node%04d" % i
    temp_filename_out = input_filename + "_node%04d001" % i + ".out"
    temp_filename_inp = input_filename + "_node%04d" % i + ".inp"
    print temp_filename_out
    print input_filename
    if options.test == False:

        # wait for file transfer to complete
        timeout = 0 # timout 
        while timeout < 600:
            if os.path.isfile(temp_filename_out) is True:
                timeout = 600
            timeout += 1
            time.sleep(1)
        os.remove(temp_filename_inp)    

    
print sys.argv[0] + ": ---- Done calculations -------------------------------"
print sys.argv[0] + ": ---- Copy files        -------------------------------"
# sleep a lot of seconds, for condor to complete file transfer
time.sleep(30)

# after execution, rename all files to normal terminology
current_dirlist = os.listdir(os.getcwd())

#foo_node01001_fort.XX -> foo001_fort.XX
#foo_node02001_fort.XX -> foo002_fort.XX
#foo_node03001_fort.XX -> foo003_fort.XX
#etc

for ii in range(nr_jobs) :
    i = ii + start_job
    pattern_nodeid = "_node%04d" %i + "001"
    pattern_nodeid_new = "0%04d" %i
    
    #attern = "^"+input_filename + "_node%04d" % i + "001"
    pattern =     input_filename + "_node%04d" % i + "001"
    
    allowed_name = re.compile(pattern).match
    print pattern
    for fname in current_dirlist:
        if allowed_name(fname):
            new_filename = fname.replace(pattern_nodeid, pattern_nodeid_new)
            print "Renameing", fname, "->", new_filename
            new_filname = os.rename(fname, new_filename)


# send mail
if options.mail != None:
    print sys.argv[0] + ": ---- Send Report to: "+options.mail+" --------------------"
    print "mail -s \'"+input_filename +": job_finished\' -c \'\' "+ options.mail+ " <_rcfluka_condor.log"
    os.system("mail -s \'"+input_filename +": job_finished\' -c \'\' "+ options.mail+ " <_rcfluka_condor.log") 

if options.clean == True:
    cleanup()


print sys.argv[0] + ": ---- Finished  -------------------------------"

