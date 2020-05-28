#!/usr/bin/python


import os
import time
import glob
import sys
import datetime
import subprocess
from cproc import *
from PIL import Image
import numpy
from operator import mul
from math import floor, ceil


#these are the default folder settings
rc_store = "/n/contefs1/backup/ysigal/"
rc_home = "/n/home12/ysigal/"
rc_exp = "GLD3_bistratified_RGC/elastic_thresh"
local_exp = rc_store + rc_exp + '/'

blocksize = 5000;
overlapsize = 2000;
gauss  =3
    
    
impath = local_exp + 'storm_merged/*'
files = glob.glob(impath)
file1 = files[0]
im = Image.open(file1)
height = im.size[0]
width = im.size[1]
block = floor(blocksize/15.8);
overlap = floor(overlapsize/15.8);
numxblocks = ceil((width-2*overlap)/block);
numyblocks = ceil((height-2*overlap)/block);

pathout =  local_exp + 'parseyfp/'
if not os.path.exists(pathout):
    os.mkdir(pathout)
counter = int(numxblocks*numyblocks)
#name of slurm partition(s) to submit jobs (serial_requeue , general , conte)
partition = 'conte,serial_requeue'
#partition = 'conte,general,serial_requeue'
   
#define and make folder to record job processing
ccproc_folder = make_proc(rc_exp + '/', rc_store)

finished = 0
num_run=0
while finished ==0:
    submitted = 0
    num_run = num_run + 1
    joblist = []
    nl = " \n "
    partition = partition + ' '
    rerun = 0
    for c in range(counter):
        var = int(c+1)
        currstat = pathout + 'B_' + str(var) + '.mat'
        if not os.path.isfile(currstat):
            matsub = ("""matlab -nosplash -nodisplay -r "arg1='""" +
                      local_exp + """'; blocksize='""" + str(blocksize) +
                      """'; overlapsize='""" + str(overlapsize) +
                      """'; slicenum='""" + str(var) +
                      """'; gauss='""" + str(gauss) +
                      """'; parseyfp_surfdist" """)
       
            job =  "Yseg_AC_" + os.path.basename(rc_exp) + "_" + str(var) + '_' + str(num_run)
            improcbatch = ccproc_folder + '/b_' + job + '.txt'
            out_fp = open(improcbatch, "w")
            out_fp.write("#!/bin/sh" + nl)
            joblist.append(job)
            jobout = (ccproc_folder + "/" + job + ".txt ")
            joberr = (ccproc_folder + "/err" + job + ".txt ")
            out_fp.write("srun -n 1 -c 2 " + matsub + nl)
            out_fp.close()
            mem_amt = str(6000 + num_run*2000)
            time_amt = str(30 + num_run*60)
            currsub = 0
            while currsub == 0:
                try:
                    sub1 = subprocess.check_output("sbatch -J " + job +
                                           " -o " + jobout +
                                           " -e " + joberr +
                                           " -n 1 -c 2 -t " + time_amt +
                                           " --mem=" + mem_amt  +
                                           " -p " + partition +
                                           improcbatch ,shell=True)
                    submitted = 1
                    currsub = 1
                    sublist = sub1.split( )
                    print "ysurfdist submitted as " + str(sublist[-1:])[2:-2]
                except:
                    print "error submitting job..."
            time.sleep(0.2)
                
                            
    if not submitted:
        finished = 1
        rerun = 1
    while rerun == 0 :
        count = 0
        time.sleep(60)
        joball = 0
        jobr = 0
        jobpd = 0
        for job in joblist:
            try:
                ql = subprocess.check_output('''squeue -u ysigal -h -n "''' + job + '''"''',shell=True)
                ql1 = subprocess.check_output('''squeue -u ysigal -h -t R -n "''' + job + '''"''',shell=True)
                ql2 = subprocess.check_output('''squeue -u ysigal -h -t PD -n "''' + job + '''"''',shell=True)
            except:
                q1 = 1
                q11 = 0
                q12 = 0
                print "error calling squeue..."
            if ql:
                joball = joball + 1
            if ql1:
                jobr = jobr + 1
            if ql2:
                jobpd = jobpd + 1


        if not joball:
            print "no current jobs are running or in queue"
            #return
            rerun = 1
        else:
            print ("at " + time.strftime("%I_%M_%p_%S",time.localtime())
                   + " there are " + str(jobr) + " jobs running and " +
                   str(jobpd) + " jobs pending")
            
        if num_run == 4:
            print 'has run 3 times without completion, ending anyway'
            rerun = 1
            finished = 1
# test for each coverslip completion
# start running thresholding and normalization
# rigid stack align
# elastic stack align
print "finshed running alignment"    
#sys.stdout = oldstdout
