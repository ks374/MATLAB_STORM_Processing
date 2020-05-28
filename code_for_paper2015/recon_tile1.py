import glob
import sys
import time, datetime, os
import subprocess
from mufit_analysis import *
from image_proc_all import *
from cproc import *

def recon_tile(mufit,xmlvar1,cproc_folder,c,xdim,ydim,redo,slicenum):
#submit one job for all bead fitting 1)separate out python fcc and fiji dcorr
    local_exp = c + "/"
    acq_folder = local_exp + 'acquisition/'
    beadfolder = local_exp + 'analysis/bead_fit/'
    ISanalysisfolder = local_exp + 'analysis/individual_sections/'
    leseg =  local_exp.split('/')
    lefolder =  str(leseg[-2:-1])[1:-1]
    coverslip = int(lefolder[-3:-1])
    exp_folder = '/'.join(leseg[:-3])

    #adjust xdim and ydim if only one image per tile
    if len(glob.glob(acq_folder + 'IRconv_' + slicenum + '_*.dax'))==1:
        xdim=str(1)
        ydim=str(1)

    #linearize bead mufit... sequential on single core, not lots of jobs
    print "running mufit_storm"
    stormdax = glob.glob(acq_folder + '*storm_' + "%03d" % int(slicenum) + '_*.dax')

    for d in stormdax:
        dvar = os.path.basename(d[:-4])
        if int(redo[1]):
            if os.path.exists(local_exp + 'analysis/bins/' + dvar + 'mlist.bin'):
                os.remove(local_exp + 'analysis/bins/' + dvar + 'mlist.bin')
            if os.path.exists(local_exp + 'analysis/bins/' + dvar + 'alist.bin'):
                os.remove(local_exp + 'analysis/bins/' + dvar + 'alist.bin')
        if not os.path.exists(local_exp + 'analysis/bins/' + dvar + 'alist.bin'):
            if os.path.exists(local_exp + 'analysis/bins/' + dvar + 'mlist.bin'):
                if os.path.getsize(local_exp + 'analysis/bins/' + dvar + 'mlist.bin')<1000:
                    os.remove(local_exp + 'analysis/bins/' + dvar + 'mlist.bin')
            
            print ("starting storm movie " + dvar + " at " +
                   time.strftime("%I_%M_%p_%S",time.localtime()))
    
            xmlvar =  xmlvar1 + dvar[:3] + dvar[-2:] + '.xml'
            binvar = local_exp + 'analysis/bins/' + dvar + 'mlist.bin'    
            mufit_run(mufit, d, binvar, xmlvar)

        numstormdax = len(stormdax)
        numalist = len(glob.glob(local_exp + 'analysis/bins/*storm_' + "%03d" % int(slicenum) + '_*alist.bin'))
        if not int(redo[1]):
            if numstormdax<=numalist:
                print str(coverslip) + " is finished"
                break
            else:
                print str(numalist) + " storm files done out of " + str(numstormdax) + "movies"

# Peform analysis when called from the command line
   
if __name__ == "__main__":
	recon_tile(sys.argv[1],sys.argv[2],sys.argv[3],
                   sys.argv[4],sys.argv[5],sys.argv[6],
                   sys.argv[7],sys.argv[8])
