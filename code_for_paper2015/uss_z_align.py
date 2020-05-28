import time
import glob
import sys
import shutil
import datetime
import subprocess
from cproc import *
from reorder import *

#before running this, *at least* inspect the wga images
#for sections that are out of order or out of focus

rc_store = "/n/contefs1/backup2/colenso"
rc_home = "/n/home12/ysigal"
rc_exp = "06_17_15_sample2/"
tiles = 1
out_of_order = 0 #if this is 1, make sure to change the starts,stops excludes below
starts = [49,100,150,201,246,299]
stops = [1,50,101,151,202,247]
excludes = [246]


exp_folder = rc_store + "/" + rc_exp + "/"
wgaswitch = 0
memlim = 8*4000-1000
#before manipulating data, first makecopy of original
if not os.path.exists(exp_folder + "unaligned_original/"):
    print 'copying original data'
    shutil.copytree(exp_folder + "unaligned/",exp_folder + "unaligned_original/")

#if sections are out of order or contain bad data, fix here
if out_of_order:
    print 'correcting order'
    reorder_sections(exp_folder,starts,stops,excludes)

#for both conv and storm images, normalize intensity across sections
#and threshold images
if not os.path.exists(exp_folder + "unaligned/for_align/"):
    print 'norm_and_thresh'
    matsub = ("""matlab -nosplash -nodisplay -r "arg1='""" +
          exp_folder + """'; wga_norm_and_thresh" """)
    subprocess.check_call(matsub ,shell=True)

#perform rigid align
if not os.path.exists(exp_folder + "rigid_align"):
    fijisub = ('/n/home12/ysigal/Fiji2.app/ImageJ-linux64 ' +
               '-Xms' + str(memlim) + 'm -Xmx' + str(memlim) + 'm -Xincgc -XX:MaxPermSize=256m ' + 
               '-XX:PermSize=256m -XX:NewRatio=5 -XX:CMSTriggerRatio=50 ' +
               ' -- --no-splash -macro fiji_z_align1.py "' +
               exp_folder + '"')
    subprocess.check_call(fijisub ,shell=True)

#open matlab **want to run crop_datastack.m
#manually determine angle of rotation and cropping
if not os.path.exists(exp_folder + "cropped"):
    matsub2 = ("""matlab-default -nosplash """)
    subprocess.check_call(matsub2 ,shell=True)
    
#perform elastic align
if not os.path.exists(exp_folder + "elastic_align"):
    fijisub = ('/n/home12/ysigal/Fiji2.app/ImageJ-linux64 ' +
               '-Xms' + str(memlim) + 'm -Xmx' + str(memlim) + 'm -Xincgc -XX:MaxPermSize=256m ' + 
               '-XX:PermSize=256m -XX:NewRatio=5 -XX:CMSTriggerRatio=50 ' +
               ' -- --no-splash -macro fiji_z_align2.py "' +
               exp_folder + '"')
    subprocess.check_call(fijisub ,shell=True)
    
#if elastic align ran, but images didn't finish saving out, finish saving images
if os.path.exists(exp_folder + "elastic_align"):
    cropfiles = os.listdir(exp_folder + "cropped/storm_merged_ds/")
    elasticfiles = os.listdir(exp_folder + "elastic_align/storm_merged_ds/")
    if len(cropfiles)>len(elasticfiles):
        fijisub = ('/n/home12/ysigal/Fiji2.app/ImageJ-linux64 ' +
                   '-Xms' + str(memlim) + 'm -Xmx' + str(memlim) + 'm -Xincgc -XX:MaxPermSize=256m ' +
                   '-XX:PermSize=256m -XX:NewRatio=5 -XX:CMSTriggerRatio=50 ' +
                   ' -- --no-splash -macro fiji_z_align3.py "' + exp_folder + '"')
        subprocess.check_call(fijisub ,shell=True)

#parameters to use for alignment
    #SIFT
        #Sigma-1.0
        #Steps-5
        #minSize-64
        #maxSize-400
        #Bins+Size = 8
        #ratio = 0.98
        #maxEpsilon = 30
        #InlierRatio = 0.03
        #NumInlier = 5
    #Elastic
        #Downsample = 0.1
        #Search Size = 350
        #Block Size = 350
        #Mesh = 36
        #geometric constraints 0.1 10 0.9
        #Sigma = 100
        #MaxDistance = 30
        #MaxRelative = 3
    #spring const try 0.5
        #Everything else stays the same
