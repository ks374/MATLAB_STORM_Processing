import glob
import sys
import time, datetime, os
import subprocess
import random
from mufit_analysis import *
from add_bead_images import sum_beads
from image_proc_all import *

def recon_bead(mufit,xmlvar1,cproc_folder,c,xdim,ydim,redo):
#submit one job for all bead fitting 1)separate out python fcc and fiji dcorr
    local_exp = c + "/"
    acq_folder = local_exp + 'acquisition/'
    beadfolder = local_exp + 'analysis/bead_fit/'
    leseg =  local_exp.split('/')
    lefolder =  str(leseg[-2:-1])[1:-1]
    coverslip = int(lefolder[-3:-1])
    cproc_folder = cproc_folder + '/' + str(coverslip)
    irpres = len(glob.glob(acq_folder + "IRbeads*0.dax"))
    vispres = len(glob.glob(acq_folder + "Visbeads*0.dax"))
    #merge and fit bead movies for registration
    sum_beads(c,redo)
    #linearize bead mufit... sequential on single core, not lots of jobs
    print "running mufit_bead"
    
    if not os.path.exists(local_exp + 'analysis/bins'):                          
        os.mkdir(local_exp + 'analysis/bins')
    stormdax = glob.glob(acq_folder + '*beads*new.dax')

    for d in stormdax:
        dvar = os.path.basename(d[:-4])
        if int(redo[0]):
            if os.path.exists(local_exp + 'analysis/bins/' + dvar + 'mlist.bin'):
                os.remove(local_exp + 'analysis/bins/' + dvar + 'mlist.bin')
            if os.path.exists(local_exp + 'analysis/bins/' + dvar + 'alist.bin'):
                os.remove(local_exp + 'analysis/bins/' + dvar + 'alist.bin')
        if not os.path.exists(local_exp + 'analysis/bins/' + dvar + 'mlist.bin'):
            print ("starting bead movie " + dvar + " at " +
                   time.strftime("%I_%M_%p_%S",time.localtime()))
            if dvar[0]=='V':
                xmlvar =  xmlvar1 + 'Visbead.xml'
            if dvar[0]=='I':
                xmlvar =  xmlvar1 + 'IRbead.xml'                        
            binvar = local_exp + 'analysis/bins/' + dvar + 'mlist.bin'

            mufit_run(mufit, d, binvar, xmlvar)
    
        numstormdax = len(glob.glob(acq_folder + '*beads*new.dax'))
        numalist = len(glob.glob(local_exp + 'analysis/bins/*newmlist.bin'))
        if numstormdax<=numalist:
            print str(coverslip) + " is finished"
            break
        else:
            print str(numalist) + " bead files done out of " + str(numstormdax) + "movies"

    # generate warping transform from beads
    if int(redo[0]):
        
        print "deleting current bead fits..."
        if os.path.exists(beadfolder):
            listf = os.listdir(beadfolder)
            for l in listf:
                os.remove(beadfolder + l)
        print "deleting tiffs and current dcorr..."
        if os.path.exists(acq_folder + 'dist_corr/'):
            listf = os.listdir(acq_folder + 'dist_corr/')
            for l in listf:
                if os.path.isfile(acq_folder + 'dist_corr/' + l):
                    os.remove(acq_folder + 'dist_corr/' + l)
                    print(l)
        proctiffs = glob.glob(acq_folder + '*ff*.tif')
        for p in proctiffs:
            print(p)
            os.remove(p)

    print "starting matlab_gen_bead_warp"
    if not os.path.isfile(beadfolder + 'Cumulative_distribution_for_registration_self.tif'):
        matsub = ("""matlab -nosplash -nodisplay -r "arg1='""" + local_exp + """'; gen_bead_warp" """)
        subprocess.Popen(matsub ,shell=True)
          
    while not os.path.isfile(beadfolder + 'Cumulative_distribution_for_registration_self.tif'):
        time.sleep(10)
    print str(coverslip)  + " matlab_gen_bead_warp is finished"  

    #test warping transform on beads (only if two bead passes)
    if vispres == 0:
        if len(glob.glob(local_exp + 'analysis/bins/IR*' +
                         '755_*' + '_0*' + 'mlist.bin' )) == 2:
            print "starting matlab_test_bead_warp"
            if not os.path.isfile(beadfolder + '750_2_647_total_offsets.tif'):

                matsub = ("""matlab -nosplash -nodisplay -r "arg1='""" + local_exp + """'; test_bead_warp" """)
                subprocess.Popen(matsub ,shell=True)
                  
            while not os.path.isfile(beadfolder + '750_2_647_total_offsets.tif'):
                time.sleep(10)
            print str(coverslip)  + " matlab_gen_bead_warp is finished"
    else:
        if len(glob.glob(local_exp + 'analysis/bins/Vis*' +
                         '540_*' + '_0*' + 'mlist.bin' )) == 2:
            print "starting matlab_test_bead_warp"
            if not os.path.isfile(beadfolder + '488_2_647_total_offsets.tif'):

                matsub = ("""matlab -nosplash -nodisplay -r "arg1='""" + local_exp + """'; test_bead_warp" """)
                subprocess.Popen(matsub ,shell=True)
                  
            while not os.path.isfile(beadfolder + '488_2_647_total_offsets.tif'):
                time.sleep(10)
            print str(coverslip)  + " matlab_gen_bead_warp is finished"

    #process ffc bead images
    print "py_image_proc_beads"
    if vispres == 0:
        if not os.path.exists(acq_folder + "avg750ffc.tif"):
            align_bead(local_exp)
    else:
        if not os.path.exists(acq_folder + "avg561ffc.tif"):
            align_bead(local_exp)
    #process distcorr for cs
    print "starting fiji_dcorr_bead"
    
    if not os.path.isfile(acq_folder + "dist_corr/distCorr.txt"):
        fijisub = ('xvfb-run -n ' + str(random.randint(40000,60000)) + ' /n/home12/ysigal/Fiji.app/ImageJ-linux64 ' +
                    '-Xms3500m -Xmx3500m -Xincgc -XX:MaxPermSize=256m ' +
                    '-Djava.util.prefs.syncInterval=2000000 ' +
                    '-XX:PermSize=256m -XX:NewRatio=5 -XX:CMSTriggerRatio=50 ' +
                    '-XX:+UseCompressedOops -- -macro fiji_distcorr_make.py "' +
                    local_exp + '"')
        subprocess.Popen(fijisub ,shell=True)
          
    while not os.path.isfile(acq_folder + "dist_corr/distCorr.txt"):
        time.sleep(10)
    print str(coverslip)  + " fiji_dcorr_bead is finished"  

# Peform analysis if called from the command line 
if __name__ == "__main__":
	recon_bead(sys.argv[1],sys.argv[2],sys.argv[3],
                   sys.argv[4],sys.argv[5],sys.argv[6],
                   sys.argv[7])
