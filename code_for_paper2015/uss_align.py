#!/usr/bin/python

from Tkinter import *
import os
import time
import glob
import sys
import datetime
import subprocess
from cproc import *
from rel_conv_ints import conv_ints


top = Tk()
top.title("Universal Serial Section Alignment")

#these are the default folder settings
ddrive = "/n/contefs1/backup2/colenso"
homedrive = "/n/home12/ysigal"
dexp = "06_17_15_sample2/"


# make an initial input screen to set new experiment path
window1 = Frame(top,height = 50,width=100)
window1.pack()
L0 = Label(window1, text="Please Enter Drive and experiment name").grid(row=1,column=2)
drivename = StringVar()
homedrivename = StringVar()
cs_start = StringVar()
cs_end = StringVar()
sm_dim = StringVar()
lg_dim = StringVar()
redovar = StringVar()
tilevar = StringVar()
expname = StringVar()
scope = StringVar()
wga488 = StringVar()

L1 = Label(window1, text="Storage Drive")
L1.grid(row=2,column=2)
E1 = Entry(window1, bd =5, textvariable=drivename,width=50)
E1.grid(row=2,column=3)
L1 = Label(window1, text="Home Drive")
L1.grid(row=3,column=2)
E1 = Entry(window1, bd =5, textvariable=homedrivename,width=50)
E1.grid(row=3,column=3)
L2 = Label(window1, text="Experiment Name")
L2.grid(row=4,column=2)
E2 = Entry(window1, bd =5, textvariable=expname,width=50)
E2.grid(row=4,column=3)

L2 = Label(window1, text="Select Coverslip start")
L2.grid(row=5,column=2)
E2 = Entry(window1, bd =5, textvariable=cs_start,width=10)
E2.grid(row=5,column=3)
L2 = Label(window1, text="Select Coverslip end")
L2.grid(row=6,column=2)
E2 = Entry(window1, bd =5, textvariable=cs_end,width=10)
E2.grid(row=6,column=3)
L2 = Label(window1, text="Small x-y dimension")
L2.grid(row=7,column=2)
E2 = Entry(window1, bd =5, textvariable=sm_dim,width=10)
E2.grid(row=7,column=3)
L2 = Label(window1, text="Large x-y dimension")
L2.grid(row=8,column=2)
E2 = Entry(window1, bd =5, textvariable=lg_dim,width=10)
E2.grid(row=8,column=3)
L2 = Label(window1, text="Redo?")
L2.grid(row=9,column=2)
E2 = Entry(window1, bd =5, textvariable=redovar,width=10)
E2.grid(row=9,column=3)
L2 = Label(window1, text="Align all tiles? (0/1)")
L2.grid(row=10,column=2)
E2 = Entry(window1, bd =5, textvariable=tilevar,width=10)
E2.grid(row=10,column=3)
L2 = Label(window1, text="Which Scope?")
L2.grid(row=13,column=1)
E2 = Entry(window1, bd =5, textvariable=scope,width=10)
E2.grid(row=13,column=2)
L2 = Label(window1, text="Is WGA in 488?")
L2.grid(row=13,column=3)
E2 = Entry(window1, bd =5, textvariable=wga488,width=10)
E2.grid(row=13,column=4)

L3 = Label(window1, text="Coverslips")
L3.grid(row=14,column=2, columnspan=2)
log = Text(window1, state='disabled',height=20,width=80, wrap='none')
log.grid(row=15,column=2, columnspan=2)

drivename.set(ddrive)
homedrivename.set(homedrive)
expname.set(dexp)
cs_start.set("1")
cs_end.set("1")
sm_dim.set("1")
lg_dim.set("4")
redovar.set("001")
tilevar.set("0")
scope.set("STORM2")
wga488.set("0")

def writeToLog(msg):
    numlines = log.index('end - 1 line').split('.')[0]
    log['state'] = 'normal'
    if int(numlines)>=1:
        log.delete(1.0, 20.0)
    if log.index('end-1c')!='1.0':
        log.insert('end', '\n')
    log.insert('end', msg)
    log['state'] = 'disabled'
def endprog() :
    top.destroy()
def findcoverslips() :
    dir_name = drivename.get() + "/" + expname.get() + "/"
    if os.path.exists(dir_name):
        global cslist
        cslist = glob.glob(dir_name + 'coverslips/*')
        if not cs_end.get() == "-1":
            #cslist = cslist[' '.join(sub_csname.get())]
            cslist = cslist[(int(cs_start.get())-1):int(cs_end.get())]
        writeToLog('\n'.join(cslist))
        
        if not cslist:
            writeToLog("no coverslips were selected")
    else:
        writeToLog("This experiment doesn't exist... try again")
def start_exp() :
    top.destroy()
    #name of experiment
    rc_exp = expname.get() + "/"
    #drive location of data
    rc_store = drivename.get() + "/"
    #drive location of code
    rc_home = homedrivename.get() + "/"
    scopename = (scope.get())
    wgaswitch = wga488.get()
    #name of slurm partition(s) to submit jobs (serial_requeue , general , conte)
    partition = 'conte,general,serial_requeue'
    #location of mufit launch program
    mufit = rc_home + 'code_addpath_rc/code_for_odyssey/python/storm-analysis/3d_daostorm/'
    #mufit = rc_home + 'code_addpath_rc/code_for_odyssey/python/storm_analysis_old/trunk/3d_daostorm/'
    #location of mufit xml files
    if scopename == 'STORM2':
        xmlvar1 = rc_home + 'code_addpath_rc/code_for_odyssey/python/storm2_xml/'
        imscale = 128
    elif scopename == 'STORM2U':
        xmlvar1 = rc_home + 'code_addpath_rc/code_for_odyssey/python/storm2_xml_ultra/'
        imscale = 128
    else:
        xmlvar1 = rc_home + 'code_addpath_rc/code_for_odyssey/python/'
        imscale = 64
    #location of mufit launch program
    pyimproc = rc_home + 'code_addpath_rc/code_for_odyssey/python/storm-analysis/image_proc_all.py '
    pyimproc = rc_home + 'code_addpath_rc/code_for_odyssey/python/storm_analysis_old/trunk/image_proc_all.py '

    #define the folder full of coverslips or other experiments containing an "acquisition" folder with storm dax files
    coverslips = cslist
    #coverslips = [rc_store + rc_exp + 'coverslips/sample1B_glyR_geph_cs08']

    #define dimensions of section (this can be found in recording from exp)
    global xdim
    global ydim
    xdim = sm_dim.get()
    ydim = lg_dim.get()
    #define and make folder to record job processing
    cproc_folder = make_cproc(rc_exp, rc_store,coverslips)

    #start recon_start.py for each coverslip
    finished = 0
    redo = str(redovar.get())
    num_run=0
    while finished ==0:
        submitted = 0
        num_run = num_run + 1
        rerun = 0
        joblist = []
        nl = " \n "
        partition = partition + ' '
        
        for c in coverslips:
            local_exp = c + "/"
            leseg =  local_exp.split('/')
            lefolder =  str(leseg[-2:-1])[1:-1]
            coverslip = int(lefolder[-3:-1])
            acq_folder = local_exp + 'acquisition/'
            analysis_folder = local_exp + 'analysis/'
            ISanalysisfolder = analysis_folder + 'individual_sections/'
            rawstorm = analysis_folder + 'stormpngs/'
            exp_folder = '/'.join(leseg[:-3])
            ccproc_folder = cproc_folder + '/' + str(coverslip)
        
            #make necessary rc folders
            if not os.path.exists(rawstorm):
                os.mkdir (rawstorm)
            if not os.path.exists(ISanalysisfolder):
                os.mkdir(ISanalysisfolder)
                #create individual section folders
                files = glob.glob(acq_folder + "IRconv*0.dax")
                if len(files)==0:
                    files = glob.glob(acq_folder + "Visconv*0.dax")
                files.sort()
                prevsequence = -1
                for file in files:
                    sectionfolder = os.path.basename(file)
                    name = os.path.basename(file)
                    idx = name.split('_')
                    index = (int(idx[1]))
                    strsequence = "%04d" % index
                    if index != prevsequence:
                        os.mkdir (ISanalysisfolder + strsequence)
                        os.mkdir (ISanalysisfolder + strsequence + "/rawimages/")
                        os.mkdir (ISanalysisfolder + strsequence + "/aligned/")
                        os.mkdir (ISanalysisfolder + strsequence + "/rawimages/for_matlab/")

                        prevsequence = index
            if not os.path.exists(local_exp + 'analysis/bins/'):                
                os.mkdir(local_exp + 'analysis/bins/')
            if not os.path.exists(acq_folder + "dist_corr/"):                                
                os.mkdir (acq_folder + "dist_corr/")
            
            print "running coverslip " + str(coverslip)
            rel_conv_ints = conv_ints(local_exp)
            rel_conv_ints = map(str,rel_conv_ints)    
            rel_conv_ints = '_'.join(rel_conv_ints)
            
            print rel_conv_ints
            #submit bead fitting code
            if int(redo[0]) or not(os.path.isfile(acq_folder + "dist_corr/distCorr.txt")):
                pysub_bead = ('python recon_bead.py ' + mufit + ' ' + xmlvar1 + ' ' + cproc_folder  +
                              ' ' + str(c) + ' ' + str(xdim) + ' ' + str(ydim) + ' '
                              + str(redovar.get()))
                ccproc_folder = cproc_folder + '/' + str(coverslip)
                job =  "Y_" + os.path.basename(exp_folder) + str(coverslip)  + "_bead_" + str(num_run)
                improcbatch = ccproc_folder + '/b_' + job + '.txt'
                out_fp = open(improcbatch, "w")
                out_fp.write("#!/bin/sh" + nl)
                out_fp.write("#SBATCH -n 1" + nl)
                out_fp.write("#SBATCH -c 1" + nl)
                out_fp.write("#SBATCH --ntasks-per-core 3" + nl)
                
                joblist.append(job)
                jobout = (ccproc_folder + "/" + job + ".txt ")
                joberr = (ccproc_folder + "/err" + job + ".txt ")
                out_fp.write("srun -n 1 -c 1 " + pysub_bead + nl)
                out_fp.close()
                mem_amt = str(2000 + num_run*2000)
                time_amt = str(num_run*300)
		currsub = 0
            	while currsub == 0:
                	try:
                		sub1 = subprocess.check_output("sbatch -J " + job +
                                               " -o " + jobout +
                                               " -e " + joberr +
                                               " -n 1 -c 1 -t " + time_amt +
                                               " --mem=" + mem_amt  +
                                               " -p " + partition +
                                               improcbatch ,shell=True)
                		submitted = 1
                    		currsub = 1
                		sublist = sub1.split( )
                		print "bead batch is submitted as " + str(sublist[-1:])[2:-2]                
                		time.sleep(0.2)
			except:
                    		print "error submitting job..."
            #submit storm molecule fitting
            files = os.listdir(ISanalysisfolder)
            if not int(tilevar.get()):
                files = files[0:1]
            for slicenum in range(len(files)):
                depend = []
                if (int(redo[1])) and (num_run==1):
                    pysub_tile = ('python recon_tile1.py ' + mufit + ' ' + xmlvar1 + ' ' + ccproc_folder  +
                                  ' ' + str(c) + ' ' + str(xdim) + ' ' + str(ydim) + ' '
                                  + str(redovar.get()) + ' ' + str(slicenum) + ' ' + str(rel_conv_ints))
                    job =  "Y1_" + os.path.basename(exp_folder) + str(coverslip)  + "_" + str(slicenum) + '_' + str(num_run)
                    improcbatch = ccproc_folder + '/b_' + job + '.txt'
                    out_fp = open(improcbatch, "w")
                    out_fp.write("#!/bin/sh" + nl)
                    joblist.append(job)
                    jobout = (ccproc_folder + "/" + job + ".txt ")
                    joberr = (ccproc_folder + "/err" + job + ".txt ")
                    out_fp.write("srun -n 1 -c 2 " + pysub_tile + nl)
                    out_fp.close()
                    mem_amt = str(2000 + num_run*1000)
                    time_amt = str(num_run*60*10)
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
                    		print "storm fitting batch is submitted as " + str(sublist[-1:])[2:-2]
                    		depend = str(sublist[-1:])[2:-2]
                    		time.sleep(0.2)
			except:
                    		print "error submitting job..."
                else:
                    stormdax = glob.glob(acq_folder + '*storm_' + "%03d" % int(slicenum) + '_*.dax')
                    numstormdax = len(stormdax)
                    numalist = len(glob.glob(local_exp + 'analysis/bins/*storm_' + "%03d" % int(slicenum) + '_*alist.bin'))
                    if not numstormdax<=numalist:
                        pysub_tile = ('python recon_tile1.py ' + mufit + ' ' + xmlvar1 + ' ' + ccproc_folder  +
                                      ' ' + str(c) + ' ' + str(xdim) + ' ' + str(ydim) + ' '
                                      + str(redovar.get()) + ' ' + str(slicenum) + ' ' + str(rel_conv_ints))
                        job =  "Y1_" + os.path.basename(exp_folder) + str(coverslip)  + "_" + str(slicenum) + '_' + str(num_run)
                        improcbatch = ccproc_folder + '/b_' + job + '.txt'
                        out_fp = open(improcbatch, "w")
                        out_fp.write("#!/bin/sh" + nl)
                        joblist.append(job)
                        jobout = (ccproc_folder + "/" + job + ".txt ")
                        joberr = (ccproc_folder + "/err" + job + ".txt ")
                        out_fp.write("srun -n 1 -c 2 " + pysub_tile + nl)
                        out_fp.close()
                        mem_amt = str(2000 + num_run*1000)
                        time_amt = str(num_run*60*10)
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
                        		print "storm fitting batch is submitted as " + str(sublist[-1:])[2:-2]
                        		depend = str(sublist[-1:])[2:-2]
                        		time.sleep(0.2)
                		except:
                    			print "error submitting job..."

                if depend:
                    if (int(redo[2])) and (num_run==1):
                        pysub_tile = ('python recon_tile2.py ' + mufit + ' ' + xmlvar1 + ' ' + ccproc_folder  +
                                      ' ' + str(c) + ' ' + str(xdim) + ' ' + str(ydim) + ' '
                                      + str(redovar.get()) + ' ' + str(slicenum) + ' ' + str(rel_conv_ints) +
                                      ' ' + str(wgaswitch))
                        job =  "Y2_" + os.path.basename(exp_folder) + str(coverslip)  + "_" + str(slicenum) + '_' + str(num_run)
                        improcbatch = ccproc_folder + '/b_' + job + '.txt'
                        out_fp = open(improcbatch, "w")
                        out_fp.write("#!/bin/sh" + nl)
                        joblist.append(job)
                        jobout = (ccproc_folder + "/" + job + ".txt ")
                        joberr = (ccproc_folder + "/err" + job + ".txt ")
                        out_fp.write("srun -n 1 -c 2 " + pysub_tile + nl)
                        out_fp.close()
                        mem_amt = str(1500 + num_run*3000)
                        time_amt = str(num_run*40 + 60)
			currsub = 0
            		while currsub == 0:
                		try:
                        		sub1 = subprocess.check_output("sbatch -J " + job +
                                                       " -o " + jobout +
                                                       " --dependency=afterok:" + depend +
                                                       " -e " + joberr +
                                                       " -n 1 -c 2 -t " + time_amt +
                                                       " --mem=" + mem_amt  +
                                                       " -p " + partition +
                                                       improcbatch ,shell=True)
                        		submitted = 1
                    			currsub = 1
                        		sublist = sub1.split( )
                        		print "storm alignment batch is submitted as " + str(sublist[-1:])[2:-2] + " depending on " + depend  
                        		time.sleep(0.2)
				except:
                    			print "error submitting job..."
                    else:
                        testobj = exp_folder + "/unaligned/storm_merged_ds/" + "%02d" % (coverslip) + "%03d" % int(slicenum) + '.tif'
                        if not os.path.exists(testobj):
                            pysub_tile = ('python recon_tile2.py ' + mufit + ' ' + xmlvar1 + ' ' + ccproc_folder  +
                                          ' ' + str(c) + ' ' + str(xdim) + ' ' + str(ydim) + ' '
                                          + str(redovar.get()) + ' ' + str(slicenum) + ' ' + str(rel_conv_ints) +
                                          ' ' + str(wgaswitch))
                            job =  "Y2_" + os.path.basename(exp_folder) + str(coverslip)  + "_" + str(slicenum) + '_' + str(num_run)
                            improcbatch = ccproc_folder + '/b_' + job + '.txt'
                            out_fp = open(improcbatch, "w")
                            out_fp.write("#!/bin/sh" + nl)
                            joblist.append(job)
                            jobout = (ccproc_folder + "/" + job + ".txt ")
                            joberr = (ccproc_folder + "/err" + job + ".txt ")
                            out_fp.write("srun -n 1 -c 2 " + pysub_tile + nl)
                            out_fp.close()
                            mem_amt = str(1500 + num_run*3000)
                            time_amt = str(num_run*40 + 60)
            		    currsub = 0
            		    while currsub == 0:
                			try:
                           			sub1 = subprocess.check_output("sbatch -J " + job +
                                                           " -o " + jobout +
                                                           " --dependency=afterok:" + depend +
                                                           " -e " + joberr +
                                                           " -n 1 -c 2 -t " + time_amt +
                                                           " --mem=" + mem_amt  +
                                                           " -p " + partition +
                                                           improcbatch ,shell=True)
                            			submitted = 1
                    				currsub = 1
                            			sublist = sub1.split( )
                            			print "storm alignment batch is submitted as " + str(sublist[-1:])[2:-2]  + " depending on " + depend                           
                            			time.sleep(0.2)
                			except:
                   				print "error submitting job..."

                else:
                    if (int(redo[2])) and (num_run==1):
                        pysub_tile = ('python recon_tile2.py ' + mufit + ' ' + xmlvar1 + ' ' + ccproc_folder  +
                                      ' ' + str(c) + ' ' + str(xdim) + ' ' + str(ydim) + ' '
                                      + str(redovar.get()) + ' ' + str(slicenum) + ' ' + str(rel_conv_ints) +
                                      ' ' + str(wgaswitch))
                        job =  "Y2_" + os.path.basename(exp_folder) + str(coverslip)  + "_" + str(slicenum) + '_' + str(num_run)
                        improcbatch = ccproc_folder + '/b_' + job + '.txt'
                        out_fp = open(improcbatch, "w")
                        out_fp.write("#!/bin/sh" + nl)
                        joblist.append(job)
                        jobout = (ccproc_folder + "/" + job + ".txt ")
                        joberr = (ccproc_folder + "/err" + job + ".txt ")
                        out_fp.write("srun -n 1 -c 2 " + pysub_tile + nl)
                        out_fp.close()
                        mem_amt = str(1500 + num_run*3000)
                        time_amt = str(num_run*40 + 60)
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
                        		print "storm alignment batch is submitted as " + str(sublist[-1:])[2:-2]                   
                        		time.sleep(0.2)
                		except:
                    			print "error submitting job..."
                    else:
                        testobj = exp_folder + "/unaligned/storm_merged_ds/" + "%02d" % (coverslip) + "%03d" % int(slicenum) + '.tif'
                        if not os.path.exists(testobj):
                            pysub_tile = ('python recon_tile2.py ' + mufit + ' ' + xmlvar1 + ' ' + ccproc_folder  +
                                          ' ' + str(c) + ' ' + str(xdim) + ' ' + str(ydim) + ' '
                                          + str(redovar.get()) + ' ' + str(slicenum) + ' ' + str(rel_conv_ints) +
                                          ' ' + str(wgaswitch))
                            job =  "Y2_" + os.path.basename(exp_folder) + str(coverslip)  + "_" + str(slicenum) + '_' + str(num_run)
                            improcbatch = ccproc_folder + '/b_' + job + '.txt'
                            out_fp = open(improcbatch, "w")
                            out_fp.write("#!/bin/sh" + nl)
                            joblist.append(job)
                            jobout = (ccproc_folder + "/" + job + ".txt ")
                            joberr = (ccproc_folder + "/err" + job + ".txt ")
                            out_fp.write("srun -n 1 -c 2 " + pysub_tile + nl)
                            out_fp.close()
                            mem_amt = str(1500 + num_run*3000)
                            time_amt = str(num_run*40 + 60)
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
                            			print "storm alignment batch is submitted as " + str(sublist[-1:])[2:-2]                            
                            			time.sleep(0.2)
                			except:
                    				print "error submitting job..."

        finished=1                
        if not submitted:
            break
        while rerun == 1 :
            count = 0
            time.sleep(180)
            joball = 0
            jobr = 0
            jobpd = 0
            for job in joblist:
                time.sleep(60)
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
                return
    # test for each coverslip completion
    # start running thresholding and normalization
    # rigid stack align
    # elastic stack align
    print "finshed running alignment"    
    #sys.stdout = oldstdout


B1 = Button(window1, text = "Find Coverslips", command = findcoverslips).grid(row=16,column=2)
B2 = Button(window1, text = "Start Alignment", command = start_exp).grid(row=16,column=3)
B3 = Button(window1, text = "Quit", command = endprog).grid(row=16,column=4)



top.mainloop()


#things to add
#1)input what channels will be acquired for storm and conv images...
#   based on this, I want to make the appropriate xml files that will be read into acq_xml_generator
#)windows to input information for dave_xmls (focal planes, movie lengths,movie names)
#several xmls, reg_bead,ffc,conv,storm (how many passes depends on if tiling)
