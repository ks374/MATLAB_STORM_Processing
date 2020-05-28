import glob
import os
import sys
from PIL import Image, ImageFilter, ImageOps
import sa_library.daxspereader as daxspereader
import numpy
import subprocess
from scipy import misc, ndimage
import re
import sa_library.arraytoimage as arraytoimage
import sa_library.i3togrid as i3togrid
import math
from rel_conv_ints import ffc_ints


def align_bead(local_exp):
    print local_exp
    leseg =  local_exp.split('/')
    rc_exp =  '/'.join(leseg[:-3]) + '/'
    storm_image_scale = int(10)
    acq_folder = local_exp + 'acquisition/'
    analysis_folder = local_exp + 'analysis/'
    ISanalysisfolder = analysis_folder + 'individual_sections/'
    rawstorm = analysis_folder + 'stormpngs/'
    rel_ffc_ints = ffc_ints(local_exp)
    #save 488ffc image -modified from Hazen's dax_to_png.py
    files = glob.glob(acq_folder + "Visffc*.dax")
    if len(files)>0:
        for file in files:
            print "File:", os.path.basename(file)

            # load dax file
            dax_file = daxspereader.DaxReader(file)
            image = dax_file.loadAFrame(2).astype(numpy.uint16)
            image = numpy.divide(image,int(rel_ffc_ints[0]))
            pilimage = Image.fromarray(image,'I;16')
            pilimage = pilimage.convert('L')
            pilimage = pilimage.rotate(-90)
            pilimage = pilimage.transpose(Image.FLIP_LEFT_RIGHT)

            # save the result
            name = os.path.basename(file)
            pilimage.save(acq_folder + "488" + name[:-4] + ".tif")

         #save 561ffc image -modified from Hazen's dax_to_png.py
           
        files = glob.glob(acq_folder + "Visffc*.dax")
        for file in files:
            print "File:", os.path.basename(file)

            # load dax file
            dax_file = daxspereader.DaxReader(file)
            image = dax_file.loadAFrame(1).astype(numpy.uint16)
            image = numpy.divide(image,int(rel_ffc_ints[1]))

            pilimage = Image.fromarray(image,'I;16')
            pilimage = pilimage.convert('L')
            pilimage = pilimage.rotate(-90)
            pilimage = pilimage.transpose(Image.FLIP_LEFT_RIGHT)

            # save the result
            name = os.path.basename(file)
            pilimage.save(acq_folder + "561" + name[:-4] + ".tif")

        #save 647Visffc image -modified from Hazen's dax_to_png.py

        files = glob.glob(acq_folder + "Visffc*.dax")
        for file in files:
            print "File:", os.path.basename(file)

            # load dax file
            dax_file = daxspereader.DaxReader(file)
            image = dax_file.loadAFrame(0).astype(numpy.uint16)
            image = numpy.divide(image,int(rel_ffc_ints[2]))
            pilimage = Image.fromarray(image,'I;16')
            pilimage = pilimage.convert('L')
            pilimage = pilimage.rotate(-90)
            pilimage = pilimage.transpose(Image.FLIP_LEFT_RIGHT)

            # save the result
            name = os.path.basename(file)
            pilimage.save(acq_folder + "647" + name[:-4] + ".tif")

        #488ffc
        for i in range(9):
            im = Image.open(acq_folder + '488Visffc_' + str(i) + '.tif')
            crop = im.crop((256,256,512,512))
            imnp = numpy.array(crop)
            imnp = numpy.reshape(imnp,(256,256,1))
            if i == 0:
                imstack = imnp
            else:  
                imstack = numpy.concatenate((imstack, imnp), axis=2)
        avgim = numpy.average(imstack,axis=2)
        pilimage = Image.fromarray(avgim)

        ffc488np = ndimage.gaussian_filter(pilimage,20)
        ffc488np[ffc488np == 0] = 1
        ffc488mean = numpy.mean(ffc488np)
        print ffc488mean
        pilimage = Image.fromarray(ffc488np)
        pilimage = pilimage.convert('L')

        pilimage.save(acq_folder + "avg488ffc.tif")
                        
        
        #561ffc
        for i in range(9):
         im = Image.open(acq_folder + '561Visffc_' + str(i) + '.tif')
         crop = im.crop((256,0,512,256))
         imnp = numpy.array(crop)
         imnp = numpy.reshape(imnp,(256,256,1))
         if i == 0:
             imstack = imnp
         else:  
             imstack = numpy.concatenate((imstack, imnp), axis=2)
        avgim = numpy.average(imstack,axis=2)
        pilimage = Image.fromarray(avgim)
        ffc561np = ndimage.gaussian_filter(pilimage,20)
        ffc561np[ffc561np == 0] = 1
        ffc561mean = numpy.mean(ffc561np)
        print ffc561mean
        pilimage = Image.fromarray(ffc561np)
        pilimage = pilimage.convert('L')
        pilimage.save(acq_folder + "avg561ffc.tif")
            
        
        #647Visffc
        for i in range(9):
         im = Image.open(acq_folder + '647Visffc_' + str(i) + '.tif')
         crop = im.crop((0,0,256,256))
         imnp = numpy.array(crop)
         imnp = numpy.reshape(imnp,(256,256,1))
         if i == 0:
             imstack = imnp
         else:  
             imstack = numpy.concatenate((imstack, imnp), axis=2)
        avgim = numpy.average(imstack,axis=2)
        pilimage = Image.fromarray(avgim)
        ffcVis647np = ndimage.gaussian_filter(pilimage,20)
        ffcVis647np[ffcVis647np == 0] = 1
        ffcVis647mean = numpy.mean(ffcVis647np)
        print ffcVis647mean
        pilimage = Image.fromarray(ffcVis647np)
        pilimage = pilimage.convert('L')
        pilimage.save(acq_folder + "avgVis647ffc.tif")


        #save 647IRffc image -modified from Hazen's dax_to_png.py
    files = glob.glob(acq_folder + "IRffc*.dax")
    if len(files)>0:
        for file in files:
            print "File:", os.path.basename(file)

            # load dax file
            dax_file = daxspereader.DaxReader(file)
            image = dax_file.loadAFrame(1).astype(numpy.uint16)
            image = numpy.divide(image,int(rel_ffc_ints[4]))

            pilimage = Image.fromarray(image,'I;16')
            pilimage = pilimage.convert('L')
            pilimage = pilimage.rotate(-90)
            pilimage = pilimage.transpose(Image.FLIP_LEFT_RIGHT)

            # save the result
            name = os.path.basename(file)
            pilimage.save(acq_folder + "647" + name[:-4] + ".tif")

        #save 750ffc image -modified from Hazen's dax_to_png.py
        files = glob.glob(acq_folder + "IRffc*.dax")
        for file in files:
            print "File:", os.path.basename(file)

            # load dax file
            dax_file = daxspereader.DaxReader(file)
            image = dax_file.loadAFrame(0).astype(numpy.uint16)
            image = numpy.divide(image,int(rel_ffc_ints[3]))
            pilimage = Image.fromarray(image,'I;16')
            pilimage = pilimage.convert('L')
            pilimage = pilimage.rotate(-90)
            pilimage = pilimage.transpose(Image.FLIP_LEFT_RIGHT)

            # save the result
            name = os.path.basename(file)
            pilimage.save(acq_folder + "750" + name[:-4] + ".tif")

         #647IRffc
        for i in range(9):
         im = Image.open(acq_folder + '647IRffc_' + str(i) + '.tif')
         crop = im.crop((0,0,256,256))
         imnp = numpy.array(crop)
         imnp = numpy.reshape(imnp,(256,256,1))
         if i == 0:
             imstack = imnp
         else:  
             imstack = numpy.concatenate((imstack, imnp), axis=2)
        avgim = numpy.average(imstack,axis=2)
        
        pilimage = Image.fromarray(avgim)
        ffcIR647np = ndimage.gaussian_filter(pilimage,20)
        ffcIR647np[ffcIR647np == 0] = 1
        ffcIR647mean = numpy.mean(ffcIR647np)
        print ffcIR647mean
        pilimage = Image.fromarray(ffcIR647np)
        pilimage = pilimage.convert('L')

        pilimage.save(acq_folder + "avgIR647ffc.tif")
        
        ###750ffc
        for i in range(9):
         im = Image.open(acq_folder + '750IRffc_' + str(i) + '.tif')
         crop = im.crop((0,256,256,512))
         imnp = numpy.array(crop)
         imnp = numpy.reshape(imnp,(256,256,1))
         if i == 0:
             imstack = imnp
         else:  
             imstack = numpy.concatenate((imstack, imnp), axis=2)
        avgim = numpy.average(imstack,axis=2)
        pilimage = Image.fromarray(avgim)
        ffc750np = ndimage.gaussian_filter(pilimage, sigma=20)
        ffc750np[ffc750np == 0] = 1
        ffc750mean = numpy.mean(ffc750np)
        print ffc750mean
        pilimage = Image.fromarray(ffc750np)
        pilimage = pilimage.convert('L')
        pilimage.save(acq_folder + "avg750ffc.tif")

    files = glob.glob(acq_folder + "IRffc*.dax")
    for i in range(9):
        if len(files)>0:
            im = Image.open(acq_folder + '647IRffc_' + str(i) + '.tif')
            crop = im.crop((0,0,256,256))
            imnp = numpy.array(crop)*ffcIR647mean
            corr = numpy.array(imnp/ffcIR647np)
        else:
            im = Image.open(acq_folder + '647Visffc_' + str(i) + '.tif')            
            crop = im.crop((0,0,256,256))
            imnp = numpy.array(crop)*ffcVis647mean
            corr = numpy.array(imnp/ffcVis647np)
        pilimage = Image.fromarray(corr)
        #pilimage.save(acq_folder + 'dist_corr/647IRffc_' + "%02d" % i + 'ds.tif')
        pilimage = pilimage.resize((2560,2560),Image.BILINEAR)
        pilimage = pilimage.convert('L')
        pilimage = ImageOps.autocontrast(pilimage,cutoff=0)
        pilimage.save(acq_folder + 'dist_corr/647ffc_' + "%02d" % i + '.tif')
        
    return


def align_tile(local_exp, slicenum, rel_conv_ints):
    leseg =  local_exp.split('/')
    rc_exp =  '/'.join(leseg[:-3]) + '/'
    storm_image_scale = int(10)
    acq_folder = local_exp + 'acquisition/'
    analysis_folder = local_exp + 'analysis/'
    ISanalysisfolder = analysis_folder + 'individual_sections/'
    rawstorm = analysis_folder + 'stormpngs/'

    #auto_to_image#
    channels = ["750storm", "647storm","488storm"]
    output_directory = rawstorm
    input_directory = analysis_folder + "bins/"
    image_max = float(256)
    print rel_conv_ints
    rel_conv_ints = map(float,rel_conv_ints.split('_'))
    print rel_conv_ints
    
    for channel in channels:
     image_base = channel
       
     def file_compare(x, y):
         x = int(re.sub(r'[^\d]', r'', x))
         y = int(re.sub(r'[^\d]', r'', y))
         return cmp(x, y)
     if os.path.isfile(input_directory):
        bin_files = [input_directory + channel + '_' + "%03d" % int(slicenum) + '_*alist.bin']
     else:
        bin_files = sorted(glob.glob(input_directory + channel + '_'
                                     + "%03d" % int(slicenum) + '_*alist.bin'), file_compare)
     manual_1st = 0
     index = manual_1st
     for file in bin_files:
         
         # Sort out the file names
     #    index = int(file.split("_")[-2])
     #    out_name = output_directory + image_base + "_7%03d" % index
         name = os.path.basename(file)
         imgname = name[:-9]
         if file == (input_directory + channel + "_0001_alist.bin"):
           out_name = output_directory + imgname
         else:
           out_name = output_directory + imgname

         
         print file
         if os.path.getsize(file)>100000:
             # Create image
             image_i3g = i3togrid.I3GData(file, scale = storm_image_scale)
             if (image_i3g.getNumberMolecules() > 0):
                 print " -> " + out_name
                 if file == (input_directory + channel + "_0001_alist.bin"):
                     index = index
                 else:
                     index += 1

                 image = image_i3g.i3To2DGridAllChannelsMerged()
                 image = numpy.transpose(image).copy()
                 image = image/image_max
                 image[(image > 1.0)] = 1.0
            
                 if (len(sys.argv) == 6):
                     print "  inverting image"
                     image = 1.0 - image
                 arraytoimage.singleColorImage(image, out_name, autoscale = False)
         else:
             print str(os.path.getsize(file)) + "bin file is less than 100kb! bad molecule list"

    convfolder = acq_folder
    #save_all_conv#
    #this is modified from Hazen's dax_to_png.py
    files = glob.glob(convfolder + 'Visconv_' + "%03d" % int(slicenum) + '_*.dax')
    if len(files)>0:        
        for file in files:
            print "File:", os.path.basename(file)

            # load dax file
            dax_file = daxspereader.DaxReader(file)
            image = dax_file.loadAFrame(2).astype(numpy.uint16)
            image = numpy.divide(image,int(rel_conv_ints[0]))

            pilimage = Image.fromarray(image,'I;16')
            pilimage = pilimage.convert('L')
            pilimage = pilimage.rotate(-90)
            pilimage = pilimage.transpose(Image.FLIP_LEFT_RIGHT)

            # save the result
            name = os.path.basename(file)
            idx = name.split('_')
            index = (int(idx[1]))
            pilimage.save(ISanalysisfolder + "%04d" % index + "/rawimages/488" + name[:-4] + ".tif")
            
            #488ffc
            for i in range(9):
                im = Image.open(acq_folder + '488Visffc_' + str(i) + '.tif')
                crop = im.crop((256,256,512,512))
                imnp = numpy.array(crop)
                imnp = numpy.reshape(imnp,(256,256,1))
                if i == 0:
                    imstack = imnp
                else:  
                    imstack = numpy.concatenate((imstack, imnp), axis=2)
            avgim = numpy.average(imstack,axis=2)
            pilimage = Image.fromarray(avgim)
        
            ffc488np = ndimage.gaussian_filter(pilimage,20)
            ffc488np[ffc488np == 0] = 1
            ffc488mean = numpy.mean(ffc488np)
            print ffc488mean
                 
        #this is modified from Hazen's dax_to_png.py
        files = glob.glob(convfolder + 'Visconv_' + "%03d" % int(slicenum) + '_*.dax')
        for file in files:
            print "File:", os.path.basename(file)

            # load dax file
            dax_file = daxspereader.DaxReader(file)
            image = dax_file.loadAFrame(1).astype(numpy.uint16)
            image = numpy.divide(image,int(rel_conv_ints[1]))

            pilimage = Image.fromarray(image,'I;16')
            pilimage = pilimage.convert('L')
            pilimage = pilimage.rotate(-90)
            pilimage = pilimage.transpose(Image.FLIP_LEFT_RIGHT)

            # save the result
            name = os.path.basename(file)
            idx = name.split('_')
            index = (int(idx[1]))
            pilimage.save(ISanalysisfolder + "%04d" % index + "/rawimages/561" + name[:-4] + ".tif")

            #561ffc
            for i in range(9):
             im = Image.open(acq_folder + '561Visffc_' + str(i) + '.tif')
             crop = im.crop((256,0,512,256))
             imnp = numpy.array(crop)
             imnp = numpy.reshape(imnp,(256,256,1))
             if i == 0:
                 imstack = imnp
             else:  
                 imstack = numpy.concatenate((imstack, imnp), axis=2)
            avgim = numpy.average(imstack,axis=2)
            pilimage = Image.fromarray(avgim)
            ffc561np = ndimage.gaussian_filter(pilimage,20)
            ffc561np[ffc561np == 0] = 1
            ffc561mean = numpy.mean(ffc561np)
            print ffc561mean
         
        #this is modified from Hazen's dax_to_png.py
        files = glob.glob(convfolder + 'Visconv_' + "%03d" % int(slicenum) + '_*.dax')
        for file in files:
            print "File:", os.path.basename(file)

            # load dax file
            dax_file = daxspereader.DaxReader(file)
            image = dax_file.loadAFrame(0).astype(numpy.uint16)
            image = numpy.divide(image,int(rel_conv_ints[2]))
            pilimage = Image.fromarray(image,'I;16')
            pilimage = pilimage.convert('L')
            pilimage = pilimage.rotate(-90)
            pilimage = pilimage.transpose(Image.FLIP_LEFT_RIGHT)

            # save the result
            name = os.path.basename(file)
            idx = name.split('_')
            index = (int(idx[1]))
            pilimage.save(ISanalysisfolder + "%04d" % index + "/rawimages/647" + name[:-4] + ".tif")

            #647Visffc
            for i in range(9):
             im = Image.open(acq_folder + '647Visffc_' + str(i) + '.tif')
             crop = im.crop((0,0,256,256))
             imnp = numpy.array(crop)
             imnp = numpy.reshape(imnp,(256,256,1))
             if i == 0:
                 imstack = imnp
             else:  
                 imstack = numpy.concatenate((imstack, imnp), axis=2)
            avgim = numpy.average(imstack,axis=2)
            pilimage = Image.fromarray(avgim)
            ffcVis647np = ndimage.gaussian_filter(pilimage,20)
            ffcVis647np[ffcVis647np == 0] = 1
            ffcVis647mean = numpy.mean(ffcVis647np)
            print ffcVis647mean

        conv_images_per_section = len(glob.glob(convfolder + 'Visconv_' + "%03d" % int(slicenum) + '_*.dax'))
        i = int(slicenum)
        for j in range (0, conv_images_per_section):
            l = "%03d" % int(slicenum)
            k= str(j)
            im = Image.open((ISanalysisfolder + "%04d" % i + "/rawimages/488Visconv_" + l + "_" + k + "_0.tif"))
            im = im.convert('L')
            crop = im.crop((256,256,512,512))
            imnp = numpy.array(crop)*ffc488mean
            corr = numpy.array(imnp/ffc488np)
            pilimage = Image.fromarray(corr)
            pilimage = pilimage.convert('L')
            pilimage.save(ISanalysisfolder + "%04d" % i + "/rawimages/for_matlab/488Visconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif")

            im = Image.open((ISanalysisfolder + "%04d" % i + "/rawimages/561Visconv_" + l + "_" + k + "_0.tif"))
            im = im.convert('L')
            crop = im.crop((256,0,512,256))
            imnp = numpy.array(crop)*ffc561mean
            corr = numpy.array(imnp/ffc561np)
            pilimage = Image.fromarray(corr)
            pilimage = pilimage.convert('L')
            pilimage.save(ISanalysisfolder + "%04d" % i + "/rawimages/for_matlab/561Visconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif")

            im = Image.open((ISanalysisfolder + "%04d" % i + "/rawimages/647Visconv_" + l + "_" + k + "_0.tif"))
            im = im.convert('L')
            crop = im.crop((0,0,256,256))
            imnp = numpy.array(crop)*ffcVis647mean
            corr = numpy.array(imnp/ffcVis647np)
            pilimage = Image.fromarray(corr)
            pilimage = pilimage.convert('L')
            pilimage.save(ISanalysisfolder + "%04d" % i + "/rawimages/for_matlab/647Visconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif")
            
    if len(files)>0:
        #this is modified from Hazen's dax_to_png.py
        files = glob.glob(convfolder + 'IRconv_' + "%03d" % int(slicenum) + '_*.dax')
        for file in files:
            print "File:", os.path.basename(file)

            # load dax file
            dax_file = daxspereader.DaxReader(file)
            image = dax_file.loadAFrame(1).astype(numpy.uint16)
            image = numpy.divide(image,int(rel_conv_ints[4]))

            pilimage = Image.fromarray(image,'I;16')
            pilimage = pilimage.convert('L')
            pilimage = pilimage.rotate(-90)
            pilimage = pilimage.transpose(Image.FLIP_LEFT_RIGHT)

            # save the result
            name = os.path.basename(file)
            idx = name.split('_')
            index = (int(idx[1]))
            pilimage.save(ISanalysisfolder + "%04d" % index + "/rawimages/647" + name[:-4] + ".tif")

             #647IRffc
            for i in range(9):
             im = Image.open(acq_folder + '647IRffc_' + str(i) + '.tif')
             crop = im.crop((0,0,256,256))
             imnp = numpy.array(crop)
             imnp = numpy.reshape(imnp,(256,256,1))
             if i == 0:
                 imstack = imnp
             else:  
                 imstack = numpy.concatenate((imstack, imnp), axis=2)
            avgim = numpy.average(imstack,axis=2)
            
            pilimage = Image.fromarray(avgim)
            ffcIR647np = ndimage.gaussian_filter(pilimage,20)
            ffcIR647np[ffcIR647np == 0] = 1
            ffcIR647mean = numpy.mean(ffcIR647np)
            print ffcIR647mean

        #this is modified from Hazen's dax_to_png.py
        files = glob.glob(convfolder + 'IRconv_' + "%03d" % int(slicenum) + '_*.dax')
        for file in files:
            print "File:", os.path.basename(file)

            # load dax file
            dax_file = daxspereader.DaxReader(file)
            image = dax_file.loadAFrame(0).astype(numpy.uint16)
            image = numpy.divide(image,int(rel_conv_ints[3]))
            pilimage = Image.fromarray(image,'I;16')
            pilimage = pilimage.convert('L')
            pilimage = pilimage.rotate(-90)
            pilimage = pilimage.transpose(Image.FLIP_LEFT_RIGHT)

            # save the result
            name = os.path.basename(file)
            idx = name.split('_')
            index = (int(idx[1]))
            pilimage.save(ISanalysisfolder + "%04d" % index + "/rawimages/750" + name[:-4] + ".tif")
               
            ###750ffc
            for i in range(9):
             im = Image.open(acq_folder + '750IRffc_' + str(i) + '.tif')
             crop = im.crop((0,256,256,512))
             imnp = numpy.array(crop)
             imnp = numpy.reshape(imnp,(256,256,1))
             if i == 0:
                 imstack = imnp
             else:  
                 imstack = numpy.concatenate((imstack, imnp), axis=2)
            avgim = numpy.average(imstack,axis=2)
            pilimage = Image.fromarray(avgim)
            ffc750np = ndimage.gaussian_filter(pilimage, sigma=20)
            ffc750np[ffc750np == 0] = 1
            ffc750mean = numpy.mean(ffc750np)
            print ffc750mean
                    
        #apply ffc to conv images
        conv_images_per_section = len(glob.glob(convfolder + 'IRconv_' + "%03d" % int(slicenum) + '_*.dax'))
        i = int(slicenum)
        for j in range (0, conv_images_per_section):
            l = "%03d" % int(slicenum)
            k= str(j)
            
            im = Image.open((ISanalysisfolder + "%04d" % i + "/rawimages/647IRconv_" + l + "_" + k + "_0.tif"))
            im = im.convert('L')
            crop = im.crop((0,0,256,256))
            imnp = numpy.array(crop)*ffcIR647mean
            corr = numpy.array(imnp/ffcIR647np)
            pilimage = Image.fromarray(corr)
            pilimage = pilimage.convert('L')
            pilimage.save(ISanalysisfolder + "%04d" % i + "/rawimages/for_matlab/647IRconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif")

            im = Image.open((ISanalysisfolder + "%04d" % i + "/rawimages/750IRconv_" + l + "_" + k + "_0.tif"))
            im = im.convert('L')
            crop = im.crop((0,256,256,512))
            imnp = numpy.array(crop)*ffc750mean
            corr = numpy.array(imnp/ffc750np)
            pilimage = Image.fromarray(corr)
            pilimage = pilimage.convert('L')
            pilimage.save(ISanalysisfolder + "%04d" % i + "/rawimages/for_matlab/750IRconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif")

            for p in range(1):
                l = "%03d" % int(slicenum)
                ps = str(p)  
                k= str(j)

                stormfile = (rawstorm + "488storm_" + l + "_" + k + "_" + ps + ".tif")
                if os.path.isfile(stormfile):

                   im = Image.open((stormfile)).convert("L")
                   imeq = ImageOps.equalize(im)
                   npeq = ndimage.gaussian_filter(imeq, sigma=1)
                   ffc488up = ndimage.zoom(ffc488np, 10, order=1)
                   imffc = numpy.multiply(npeq,numpy.log10(ffc488mean))
                   imffc = numpy.divide(imffc,numpy.log10(ffc488up))
                   pilimage = Image.fromarray(imffc)
                   pilimage = pilimage.convert('L')
                   imadj = ImageOps.autocontrast(pilimage)
                   pilimage.save(ISanalysisfolder + "%04d" % i + "/rawimages/for_matlab/488storm_" + "%03d" % i + "_" + "%02d" % j + "_" + ps + ".tif")

                stormfile = (rawstorm + "561storm_" + l + "_" + k + "_" + ps + ".tif")
                if os.path.isfile(stormfile):

                   im = Image.open((stormfile)).convert("L")
                   imeq = ImageOps.equalize(im)
                   npeq = ndimage.gaussian_filter(imeq, sigma=1)
                   ffc561up = ndimage.zoom(ffc561np, 10, order=1)
                   imffc = numpy.multiply(npeq,numpy.log10(ffc561mean))
                   imffc = numpy.divide(imffc,numpy.log10(ffc561up))
                   pilimage = Image.fromarray(imffc)
                   pilimage = pilimage.convert('L')
                   imadj = ImageOps.autocontrast(pilimage)
                   imadj.save(ISanalysisfolder + "%04d" % i + "/rawimages/for_matlab/561storm_" + "%03d" % i + "_" + "%02d" % j + "_" + ps + ".tif")

                stormfile = (rawstorm + "647storm_" + l + "_" + k + "_" + ps + ".tif")
                if os.path.isfile(stormfile):

                   im = Image.open((stormfile)).convert("L")
                   imeq = ImageOps.equalize(im)
                   npeq = ndimage.gaussian_filter(imeq, sigma=1)
                   ffcIR647up = ndimage.zoom(ffcIR647np, 10, order=1)
                   imffc = numpy.multiply(npeq,numpy.log10(ffcIR647mean))
                   imffc = numpy.divide(imffc,numpy.log10(ffcIR647up))
                   pilimage = Image.fromarray(imffc)
                   pilimage = pilimage.convert('L')
                   imadj = ImageOps.autocontrast(pilimage)
                   imadj.save(ISanalysisfolder + "%04d" % i + "/rawimages/for_matlab/647storm_" + "%03d" % i + "_" + "%02d" % j + "_" + ps + ".tif")

                stormfile = (rawstorm + "750storm_" + l + "_" + k + "_" + ps + ".tif")
                if os.path.isfile(stormfile):

                   im = Image.open((stormfile)).convert("L")
                   imeq = ImageOps.equalize(im)
                   npeq = ndimage.gaussian_filter(imeq, sigma=1)
                   ffc750up = ndimage.zoom(ffc750np, 10, order=1)
                   imffc = numpy.multiply(npeq,numpy.log10(ffc750mean))
                   imffc = numpy.divide(imffc,numpy.log10(ffc750up))
                   pilimage = Image.fromarray(imffc)
                   pilimage = pilimage.convert('L')
                   imadj = ImageOps.autocontrast(pilimage)
                   imadj.save(ISanalysisfolder + "%04d" % i + "/rawimages/for_matlab/750storm_" + "%03d" % i + "_" + "%02d" % j + "_" + ps + ".tif")
    return
