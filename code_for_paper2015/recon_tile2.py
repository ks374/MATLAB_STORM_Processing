import glob
import sys
import time, datetime, os,random
import shutil
import subprocess
from mufit_analysis import *
from image_proc_all import *
from cproc import *

def recon_tile(mufit,xmlvar1,cproc_folder,c,xdim,ydim,redo,slicenum,rel_conv_ints, wgaswitch):
#submit one job for all bead fitting 1)separate out python fcc and fiji dcorr
    local_exp = c + "/"
    print local_exp
    acq_folder = local_exp + 'acquisition/'
    beadfolder = local_exp + 'analysis/bead_fit/'
    ISanalysisfolder = local_exp + 'analysis/individual_sections/'
    leseg =  local_exp.split('/')
    print leseg
    lefolder =  str(leseg[-2:-1])[1:-1]
    print lefolder
    coverslip = int(lefolder[-3:-1])
    print coverslip
    exp_folder = '/'.join(leseg[:-3])
    print exp_folder
    irpres = len(glob.glob(acq_folder + "IRconv*0.dax"))
    vispres = len(glob.glob(acq_folder + "Visconv*0.dax"))
    #adjust xdim and ydim if only one image per tile
    if len(glob.glob(acq_folder + 'IRconv_' + "%03d" % int(slicenum) + '_*.dax'))==1:
        xdim=str(1)
        ydim=str(1)

    #linearize bead mufit... sequential on single core, not lots of jobs
    stormdax = glob.glob(acq_folder + '*storm_' + "%03d" % int(slicenum) + '_*.dax')
    numalist = len(glob.glob(local_exp + 'analysis/bins/*storm_' + "%03d" % int(slicenum) + '_*alist.bin'))
    numstormdax = len(stormdax)
    while not numstormdax<=numalist:
        numalist = len(glob.glob(local_exp + 'analysis/bins/*storm_' + "%03d" % int(slicenum) + '_*alist.bin'))
        time.sleep(10)
        print "waiting until fitting is done"
        print "currently " + str(numalist) + " storm files done out of " + "%03d" % int(slicenum) + "movies"
    print "mufit done!"
    while not os.path.isfile(beadfolder + 'Cumulative_distribution_for_registration_self.tif'):
        time.sleep(10)
        print "waiting for beadfit"
    print str(coverslip)  + " bead_warp is present"  

    while not os.path.isfile(acq_folder + "dist_corr/distCorr.txt"):
        time.sleep(10)
        print "waiting for distcorr"
    print str(coverslip)  + " ffc distcorr is present"  
        
    #processes storm and convimages to tiff files
    print "py_image_proc_tile"
    if int(redo[2]):
        align_tile(local_exp, str(slicenum),rel_conv_ints)      
        subfiles = glob.glob(ISanalysisfolder + "%04d" % int(slicenum)  + "/aligned/*")
        for s in subfiles:
            if os.path.isfile(s):
                os.unlink(s)
            if os.path.isdir(s):
                shutil.rmtree(s)
    #elif not os.path.exists(local_exp + "analysis/stormpngs/488storm_" + slicenum + "_0_1.png"):        
    #        align_tile(local_exp, slicenum, imscale,rel_conv_ints)            
    #        subfiles = glob.glob(ISanalysisfolder + "%04d" % int(slicenum) + "/aligned/*_3.tif")
    #        for s in subfiles:
    #            os.remove(s)
    elif not os.path.exists(exp_folder + "/unaligned/storm_merged_ds/" +
                            "%02d" % (coverslip) + "%03d" % int(slicenum) + '.tif'):        
        align_tile(local_exp, str(slicenum),rel_conv_ints)            
        subfiles = glob.glob(ISanalysisfolder + "%04d" % int(slicenum) + "/aligned/*")
        for s in subfiles:
            if os.path.isfile(s):
                os.unlink(s)
            if os.path.isdir(s):
                shutil.rmtree(s)
            
    #else:
    #    if not os.path.exists(ISanalysisfolder + "%04d" % int(slicenum) + "/rawimages/for_matlab/750storm_" +
    #                      "%03d" % int(slicenum) + "_00_1.tif"):        
    #        align_tile(local_exp, slicenum, imscale,rel_conv_ints)
    #        subfiles = glob.glob(ISanalysisfolder + "%04d" % int(slicenum) + "/aligned/*_3.tif")
    #        for s in subfiles:
    #            os.remove(s)
 
    try:
        #apply warping transform to images
        print "starting matlab_image_warp"
        if irpres == 0:  
            if not len(glob.glob(ISanalysisfolder + "%04d" % int(slicenum) + "/aligned/*_3.tif"))==(int(xdim)*int(ydim)):
                matsub = ("""matlab -nosplash -nodisplay -r "arg1='""" +
                          local_exp + """'; slicenum='""" + str(slicenum) +
                          """'; image_chrom_align_visonly" """)
                subprocess.check_call(matsub ,shell=True)
        else:
            if not len(glob.glob(ISanalysisfolder + "%04d" % int(slicenum) + "/aligned/*_3.tif"))==(int(xdim)*int(ydim)):
                matsub = ("""matlab -nosplash -nodisplay -r "arg1='""" +
                          local_exp + """'; slicenum='""" + str(slicenum) +
                          """'; image_chrom_align" """)
                subprocess.check_call(matsub ,shell=True)            
              
        #while not len(glob.glob(ISanalysisfolder + "%04d" % int(slicenum) + "/aligned/*_3.tif")) == (int(xdim)*int(ydim)):
        #    time.sleep(10)
        print str(coverslip)  + " matlab_image_warp is finished"  
             
        #apply distortion correction and tiling in fiji
        print "starting fiji_dcorr"
        testobj = exp_folder + "/unaligned/storm_merged_ds/" + "%02d" % (coverslip) + "%03d" % int(slicenum) + '.tif'
        print testobj
        if int(redo[2]):
            if os.path.exists(testobj):
                os.remove(testobj)
        if not os.path.exists(testobj):
            if (int(xdim)*int(ydim))==1:
                print "only one field"
                fijisub = ('xvfb-run -n ' + str(random.randint(40000,60000)) + ' /n/home12/ysigal/Fiji.app/ImageJ-linux64 ' +
                           '-Xms6000m -Xmx6000m -Xincgc -XX:MaxPermSize=256m ' + 
                           '-Djava.util.prefs.syncInterval=2000000 ' + 
                           '-XX:PermSize=256m -XX:NewRatio=5 -XX:CMSTriggerRatio=50 ' +
                           '-XX:+UseCompressedOops -- -macro fiji_distcorr1field.py "' +
                           "%03d" % int(slicenum) + ' ' + local_exp + ' ' + str(xdim) + ' ' +
                           str(ydim) + ' ' + cproc_folder + ' ' + str(wgaswitch) + '"')
                subprocess.check_call(fijisub ,shell=True)
            else:
                fijisub = ('xvfb-run -n ' + str(random.randint(40000,60000)) + ' /n/home12/ysigal/Fiji.app/ImageJ-linux64 ' +
                           '-Xms6000m -Xmx6000m -Xincgc -XX:MaxPermSize=256m ' + 
                           '-Djava.util.prefs.syncInterval=2000000 ' + 
                           '-XX:PermSize=256m -XX:NewRatio=5 -XX:CMSTriggerRatio=50 ' +
                           '-XX:+UseCompressedOops -- -macro fiji_distcorr.py "' +
                           "%03d" % int(slicenum) + ' ' + local_exp + ' ' + str(xdim) + ' ' +
                           str(ydim) + ' ' + cproc_folder + ' ' + str(wgaswitch) + '"')
                subprocess.check_call(fijisub ,shell=True)
              
    #    while not os.path.exists(testobj):
    #        time.sleep(10)
        print str(coverslip)  + " fiji_dcorr_tile number " + slicenum + " is finished! XY alignment is done!"  
    except:
        print "Failed"
# Peform analysis when called from the command line
   
if __name__ == "__main__":
	recon_tile(sys.argv[1],sys.argv[2],sys.argv[3],
                   sys.argv[4],sys.argv[5],sys.argv[6],
                   sys.argv[7],sys.argv[8],sys.argv[9],sys.argv[10])
