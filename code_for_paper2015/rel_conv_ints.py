#determine the relative intensity for conventional images in each channel for each coverslip
import glob
import library.daxspereader as daxspereader
import numpy

def conv_ints(local_exp):
    leseg =  local_exp.split('/')
    rc_exp =  '/'.join(leseg[:-3]) + '/'
    storm_image_scale = int(10)
    acq_folder = local_exp + 'acquisition/'
    analysis_folder = local_exp + 'analysis/'
    ISanalysisfolder = analysis_folder + 'individual_sections/'
    rawstorm = analysis_folder + 'stormpngs/'
    #save_all_conv#
    #this is modified from Hazen's dax_to_png.py
    files = glob.glob(acq_folder + 'Visconv_' + '*.dax')
    ## files = files[0]
    #print files
    if len(files)>0:
        cnt = 0
        aperc_v488 = [0]*len(files)
        aperc_v561 = [0]*len(files)
        aperc_v647 =[0]*len(files)
        for file in files:
            # load dax file
            ## print file
            dax_file = daxspereader.DaxReader(file)
            image = dax_file.loadAFrame(0).astype(numpy.uint16)
            im2 = numpy.split(image,2,0)
            im3 = numpy.split(im2[0],2,1)
            im4 = numpy.split(im2[1],2,1)
            aperc_v647[cnt] = numpy.percentile(im3[0], 99.999)

            image = dax_file.loadAFrame(1).astype(numpy.uint16)
            im2 = numpy.split(image,2,0)
            im3 = numpy.split(im2[0],2,1)
            im4 = numpy.split(im2[1],2,1)
            aperc_v561[cnt] = numpy.percentile(im4[0], 99.999)

            image = dax_file.loadAFrame(2).astype(numpy.uint16)
            im2 = numpy.split(image,2,0)
            im3 = numpy.split(im2[0],2,1) #im3[0]=647channel im3[1]=750
            im4 = numpy.split(im2[1],2,1) #im4[0]=561channel im3[1]=488
            aperc_v488[cnt] = numpy.percentile(im4[1], 99.999)

            cnt = cnt+1

    files = glob.glob(acq_folder + 'IRconv_' + '*.dax')
    ## files = files[0]
    #print files
    if len(files)>0:
        cnt = 0
        aperc_IR750 = [0]*len(files)
        aperc_IR647 = [0]*len(files)
        for file in files:
            # load dax file
            ## print file
            dax_file = daxspereader.DaxReader(file)
            image = dax_file.loadAFrame(0).astype(numpy.uint16)
            im2 = numpy.split(image,2,0)
            im3 = numpy.split(im2[0],2,1)
            im4 = numpy.split(im2[1],2,1)
            aperc_IR750[cnt] = numpy.percentile(im3[1], 99.999)

            image = dax_file.loadAFrame(1).astype(numpy.uint16)
            im2 = numpy.split(image,2,0)
            im3 = numpy.split(im2[0],2,1) #im3[0]=647channel im3[1]=750
            im4 = numpy.split(im2[1],2,1) #im4[0]=561channel im3[1]=488
            aperc_IR647[cnt] = numpy.percentile(im3[0], 99.999)

            cnt = cnt+1
    rel_conv_ints = [0]*5
    
    rel_conv_ints[0] =  numpy.mean(aperc_v488)/256
    rel_conv_ints[1] =  numpy.mean(aperc_v561)/256
    rel_conv_ints[2] =  numpy.mean(aperc_v647)/256
    rel_conv_ints[3] =  numpy.mean(aperc_IR750)/256
    rel_conv_ints[4] =  numpy.mean(aperc_IR647)/256
            
    return (rel_conv_ints)

def ffc_ints(local_exp):
    leseg =  local_exp.split('/')
    rc_exp =  '/'.join(leseg[:-3]) + '/'
    storm_image_scale = int(10)
    acq_folder = local_exp + 'acquisition/'
    analysis_folder = local_exp + 'analysis/'
    ISanalysisfolder = analysis_folder + 'individual_sections/'
    rawstorm = analysis_folder + 'stormpngs/'
    #save_all_conv#
    #this is modified from Hazen's dax_to_png.py
    files = glob.glob(acq_folder + 'Visffc_' + '*.dax')
    ## files = files[0]
    #print files
    if len(files)>0:
        cnt = 0
        aperc_v488 = [0]*len(files)
        aperc_v561 = [0]*len(files)
        aperc_v647 = [0]*len(files)
        for file in files:
            # load dax file
            ## print file
            dax_file = daxspereader.DaxReader(file)
            image = dax_file.loadAFrame(0).astype(numpy.uint16)
            im2 = numpy.split(image,2,0)
            im3 = numpy.split(im2[0],2,1)
            im4 = numpy.split(im2[1],2,1)
            aperc_v647[cnt] = numpy.percentile(im3[0], 99.95)

            image = dax_file.loadAFrame(1).astype(numpy.uint16)
            im2 = numpy.split(image,2,0)
            im3 = numpy.split(im2[0],2,1)
            im4 = numpy.split(im2[1],2,1)
            aperc_v561[cnt] = numpy.percentile(im4[0], 99.95)

            image = dax_file.loadAFrame(2).astype(numpy.uint16)
            im2 = numpy.split(image,2,0)
            im3 = numpy.split(im2[0],2,1) #im3[0]=647channel im3[1]=750
            im4 = numpy.split(im2[1],2,1) #im4[0]=561channel im3[1]=488
            aperc_v488[cnt] = numpy.percentile(im4[1], 99.95)

            cnt = cnt+1

    files = glob.glob(acq_folder + 'IRffc_' + '*.dax')
    ## files = files[0]
    #print files
    if len(files)>0:
        cnt = 0
        aperc_IR750 = [0] * len(files)
        aperc_IR647 = [0]* len(files)
        for file in files:
            # load dax file
            ## print file
            dax_file = daxspereader.DaxReader(file)
            image = dax_file.loadAFrame(0).astype(numpy.uint16)
            im2 = numpy.split(image,2,0)
            im3 = numpy.split(im2[0],2,1)
            im4 = numpy.split(im2[1],2,1)
            aperc_IR750[cnt] = numpy.percentile(im3[1], 99.95)

            image = dax_file.loadAFrame(1).astype(numpy.uint16)
            im2 = numpy.split(image,2,0)
            im3 = numpy.split(im2[0],2,1) #im3[0]=647channel im3[1]=750
            im4 = numpy.split(im2[1],2,1) #im4[0]=561channel im3[1]=488
            aperc_IR647[cnt] = numpy.percentile(im3[0], 99.95)

            cnt = cnt+1
    rel_ffc_ints = [0] * 5
    
    rel_ffc_ints[0] =  numpy.mean(aperc_v488)/256
    rel_ffc_ints[1] =  numpy.mean(aperc_v561)/256
    rel_ffc_ints[2] =  numpy.mean(aperc_v647)/256
    rel_ffc_ints[3] =  numpy.mean(aperc_IR750)/256
    rel_ffc_ints[4] =  numpy.mean(aperc_IR647)/256
            
    return (rel_ffc_ints)
