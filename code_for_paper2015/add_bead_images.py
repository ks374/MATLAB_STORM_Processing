##correct bead images
#1-load in bead dax
#2-crop each frame and save each quadrant as a tiff
#3-apply lenscorrection in fiji to each image and resave
#4-load tiff for each frame and recombine
#5-resave as dax
#

import glob
import os
from PIL import Image, ImageFilter
import library.datareader as daxspereader
import library.daxwriter as daxwriter
#import daxspereader
import subprocess
import numpy

#define variables


def sum_beads(c,redo):
    print "running sum_beads"
    
    local_exp = c + "/"
    cvar = os.path.basename(c)
    acq_folder = local_exp + 'acquisition/'
    files = glob.glob(acq_folder + "*beads*.dax")
    #print files
    for file in files:
        #print file
        name = os.path.basename(file)
        name = name[:-4]
        if int(redo[0]):
            if name[-3:]=='new':
                os.remove(file)
            else:
                pass
        
        if name[-3:]=='new':
            pass
            #print "this is single frame dax"
        elif os.path.exists(file[:-4] + 'new.dax'):
            pass
            #print "this was already processed"
        else:
            for i in range (20):

                # load dax file
                dax_file = daxspereader.DaxReader(file)
                image = dax_file.loadAFrame(i).astype(numpy.float)
                if file == files[0]:
                    maximage = numpy.max(image)
                # create the image
                img = numpy.zeros((image.shape[0], image.shape[1]), numpy.uint32)
                if i==0:
                    img2 = numpy.zeros((image.shape[0], image.shape[1]), numpy.uint32)
                #image = 255.0 * image/maximage
                img[:,:] = image.astype(numpy.uint32)
                img2 = img2+img
            img3 = img2/20

            dax_file = daxwriter.DaxWriter(acq_folder + name + "new.dax", 0, 0)
            dax_file.addFrame(img3)
            dax_file.close()
            print "summing File:", os.path.basename(file)
