import os,sys,time
import glob
import shutil
from ij import IJ  
import math
from ij.process import ImageStatistics as IS
from ij import Prefs
from register_virtual_stack import Register_Virtual_Stack_MT
from register_virtual_stack import Transform_Virtual_Stack_MT
import signal


argv = getArgument()

arg_array = argv.split(" ")


local_exp = arg_array[1]
xdim = arg_array[2]
ydim = arg_array[3]
secnum=arg_array[0]
cproc = arg_array[4]
wgaswitch = arg_array[5]
if not local_exp[-1:]=="/":
	local_exp = local_exp + "/"
if not cproc[-1:]=="/":
	cproc = cproc + "/"
f1 = open(cproc + secnum + '_fiji.txt','w+')
print >>f1, str(arg_array)
print >>f1, str(argv)
print >>f1, "starting fiji from the macros folder"

for n in range(1):

	xdim = int(xdim)
	ydim = int(ydim)
	wgaswitch = int(wgaswitch)
	secnum = int(secnum)
	acq_folder = local_exp + 'acquisition/'
	ffc_folder = acq_folder
	analysisfolder = local_exp + "analysis/"
	ISanalysisfolder = analysisfolder + "individual_sections/"
	conv_images_per_section = (xdim * ydim)

	files = glob.glob(ISanalysisfolder + "*")
	files = files[(secnum-1):secnum]

	print >>f1,  "I647_2560"
	distcorr_src = ffc_folder + "dist_corr/distCorr.txt"
	for strseq in [secnum]:
		print >>f1, strseq
		strsequence = "%04d" % strseq
		strseq2 = "%03d" % strseq
		print >>f1, strsequence
		if not os.path.exists(ISanalysisfolder + strsequence + "/aligned/after_dist_corr/"):
			os.mkdir (ISanalysisfolder + strsequence + "/aligned/after_dist_corr/")
		subfiles = glob.glob(ISanalysisfolder + strsequence + "/aligned/after_dist_corr/*")
		print >>f1, subfiles
		for s in subfiles:
			os.remove(s)
		subfiles = glob.glob(ISanalysisfolder + strsequence + "/aligned/*_3.tif")
		print >>f1, subfiles
		for s in subfiles:
			os.remove(s)
		dst = ISanalysisfolder + strsequence + "/aligned/"
		shutil.copy(distcorr_src, dst)

		if os.path.exists(ISanalysisfolder + strsequence + "/aligned/488storm_" +
                   		strseq2 + "_00_0.tif"):
			IJ.run("lenscorrection Distortion Correction", "calibration=" + ISanalysisfolder
                   		+ strsequence + "/aligned/ number_of_images=52 first_image=488storm_" +
                   		strseq2 + "_00_0.tif power_of_polynomial_kernel=3 lambda=0.00010000 " +
                   		"apply_correction_to_images visualize what=load file_name=distCorr.txt " +
                   		"initial_gaussian_blur=1.00 steps_per_scale_octave=7 minimum_image_size=64 " +
                   		"maximum_image_size=400 feature_descriptor_size=8 feature_descriptor_orientation_bins=8 " +
                   		"closest/next_closest_ratio=0.92 maximal_alignment_error=32 inlier_ratio=0.20 " +
                   		"expected_transformation=Rigid target=" + ISanalysisfolder + strsequence +
                   		"/aligned/after_dist_corr/");
		else:
			IJ.run("lenscorrection Distortion Correction", "calibration=" + ISanalysisfolder
                   		+ strsequence + "/aligned/ number_of_images=52 first_image=488storm_" +
                   		strseq2 + "_00_1.tif power_of_polynomial_kernel=3 lambda=0.00010000 " +
                   		"apply_correction_to_images visualize what=load file_name=distCorr.txt " +
                   		"initial_gaussian_blur=1.00 steps_per_scale_octave=7 minimum_image_size=64 " +
                   		"maximum_image_size=400 feature_descriptor_size=8 feature_descriptor_orientation_bins=8 " +
                   		"closest/next_closest_ratio=0.92 maximal_alignment_error=32 inlier_ratio=0.20 " +
                   		"expected_transformation=Rigid target=" + ISanalysisfolder + strsequence +
                   		"/aligned/after_dist_corr/");
		print >>f1, "distcorr for section", file
	
	files = len(files)	
	for i in [secnum]:
		if not os.path.isdir(ISanalysisfolder + "%04d" % i + "/aligned/conv_561/"):
			os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/prestitch/")
			os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/conv_561/")
			os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/conv_561/out/")
			os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/conv_561/mc/")
			os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/storm_merged/")
			os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/storm_merged/out/")
			os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/storm_merged/mc/")
			os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/conv_561adj/")
			os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/conv_561adj/out/")
			os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/conv_561adj/mc/")
			os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/conv_merged/")
			os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/conv_merged/out/")
			os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/conv_merged/mc/")
			#os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/storm_561/")
			#os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/storm_561/out/")
			#os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/storm_561/mc/")
			for j in range (1, (xdim+1)):
				os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/conv_561/c_" + "%01d" % j + "/")
				os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/conv_561adj/c_" + "%01d" % j + "/") 
				os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/conv_merged/c_" + "%01d" % j + "/")
				#os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/storm_561/c_" + "%01d" % j + "/")
				os.mkdir (ISanalysisfolder + "%04d" % i + "/aligned/storm_merged/c_" + "%01d" % j + "/") 
		## this section merges 2 storm images and uses the average
		for j in range (0, conv_images_per_section):
			k= str(j)
			if os.path.isfile((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/488storm_" + "%03d" % i + "_" + "%02d" % j + "_0.tif")):
				imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/488storm_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))
				imp.show()
				IJ.run(imp, "Canvas Size...", "width=2350 height=2350 position=Center zero");
				IJ.run(imp, "Canvas Size...", "width=2560 height=2560 position=Center zero");
				#imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/488storm_" + "%03d" % i + "_" + "%02d" % j + "_1.tif"))
				imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/488storm_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))
				imp.show()
				IJ.run(imp, "Canvas Size...", "width=2350 height=2350 position=Center zero");
				IJ.run(imp, "Canvas Size...", "width=2560 height=2560 position=Center zero");
				IJ.run(imp, "Images to Stack", "name=Stack title=[] use");
				imp = IJ.getImage()  
				IJ.run(imp, "Z Project...", "start=1 stop=2 projection=[Average Intensity]");
				imp = IJ.getImage()  
			else:
				imp = IJ.createImage("Untitled", "8-bit Black", 2560, 2560, 1);
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/488storm_" + "%03d" % i + "_" + "%02d" % j + ".tif"));
			IJ.run("Close All", "");
			
			#if os.path.isfile((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/561storm_" + "%02d" % i + "_" + "%02d" % j + "_0.tif")):
			#	imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/561storm_" + "%02d" % i + "_" + "%02d" % j + "_0.tif"))
			#	imp.show()
			#	IJ.run(imp, "Canvas Size...", "width=2350 height=2350 position=Center zero");
			#	IJ.run(imp, "Canvas Size...", "width=2560 height=2560 position=Center zero");
			#	imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/561storm_" + "%02d" % i + "_" + "%02d" % j + "_1.tif"))
			#	imp.show()
			#	IJ.run(imp, "Canvas Size...", "width=2350 height=2350 position=Center zero");
			#	IJ.run(imp, "Canvas Size...", "width=2560 height=2560 position=Center zero");
			#	IJ.run(imp, "Images to Stack", "name=Stack title=[] use");
			#	imp = IJ.getImage()  
			#	IJ.run(imp, "Z Project...", "start=1 stop=2 projection=[Average Intensity]");
			#	imp = IJ.getImage()  
			#else:
			#	imp = IJ.createImage("Untitled", "8-bit Black", 2560, 2560, 1);
			#IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/561storm_" + "%02d" % i + "_" + "%02d" % j + ".tif"));
			#IJ.run("Close All", ""); 
			
			if os.path.isfile((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/647storm_" + "%03d" % i + "_" + "%02d" % j + "_0.tif")):
				imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/647storm_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))
				imp.show()
				IJ.run(imp, "Canvas Size...", "width=2350 height=2350 position=Center zero");
				IJ.run(imp, "Canvas Size...", "width=2560 height=2560 position=Center zero");
				#imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/647storm_" + "%02d" % i + "_" + "%02d" % j + "_1.tif"))
				imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/647storm_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))
				imp.show()
				IJ.run(imp, "Canvas Size...", "width=2350 height=2350 position=Center zero");
				IJ.run(imp, "Canvas Size...", "width=2560 height=2560 position=Center zero");
				IJ.run(imp, "Images to Stack", "name=Stack title=[] use");
				imp = IJ.getImage()  
				IJ.run(imp, "Z Project...", "start=1 stop=2 projection=[Average Intensity]");
				imp = IJ.getImage()  
			else:
				imp = IJ.createImage("Untitled", "8-bit Black", 2560, 2560, 1);
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/647storm_" + "%03d" % i + "_" + "%02d" % j + ".tif"));
			IJ.run("Close All", ""); 
			
			if os.path.isfile((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/750storm_" + "%03d" % i + "_" + "%02d" % j + "_0.tif")):
				imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/750storm_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))
				imp.show()
				IJ.run(imp, "Canvas Size...", "width=2350 height=2350 position=Center zero");
				IJ.run(imp, "Canvas Size...", "width=2560 height=2560 position=Center zero");
				#imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/750storm_" + "%02d" % i + "_" + "%02d" % j + "_1.tif"))
				imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/750storm_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))
				imp.show()
				IJ.run(imp, "Canvas Size...", "width=2350 height=2350 position=Center zero");
				IJ.run(imp, "Canvas Size...", "width=2560 height=2560 position=Center zero");
				IJ.run(imp, "Images to Stack", "name=Stack title=[] use");
				imp = IJ.getImage()  
				IJ.run(imp, "Z Project...", "start=1 stop=2 projection=[Average Intensity]");
				imp = IJ.getImage()  
			else:
				imp = IJ.createImage("Untitled", "8-bit Black", 2560, 2560, 1);
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/750storm_" + "%03d" % i + "_" + "%02d" % j + ".tif"));
			IJ.run("Close All", "");
			
                        if not os.path.isfile((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/750IRconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif")):
                                imp = IJ.createImage("Untitled", "8-bit Black", 2560, 2560, 1);
                                IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/750IRconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))

			if not os.path.isfile((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/647IRconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif")):
                                imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/647Visconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))
                                IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/647IRconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))

			if not os.path.isfile((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/647Visconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif")):
                                imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/647Irconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))
                                IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/647Visconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))

                        if not os.path.isfile((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/561Visconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif")):
                                imp = IJ.createImage("Untitled", "8-bit Black", 2560, 2560, 1);
                                IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/561Visconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))

                        if not os.path.isfile((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/488Visconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif")):
                                imp = IJ.createImage("Untitled", "8-bit Black", 2560, 2560, 1);
                                IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/488Visconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))

			## this section crops all images the same from 2560 to 2450
			imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/750storm_" + "%03d" % i + "_" + "%02d" % j + ".tif"))
			imp.show()
			imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/647storm_" + "%03d" % i + "_" + "%02d" % j + ".tif"))
			imp.show()
			imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/488storm_" + "%03d" % i + "_" + "%02d" % j + ".tif"))
			imp.show()
			#imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/488storm_" + "%02d" % i + "_" + "%02d" % j + "_01.tif"))
			#imp.show()
			#imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/561storm_" + "%02d" % i + "_" + "%02d" % j + ".tif"))
			#imp.show()
			imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/750IRconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))
			imp.show()  
			imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/647IRconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))
			imp.show() 
			if wgaswitch:
				imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/561Visconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))
				imp.show()  
				imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/488Visconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))
				imp.show()  
			if not wgaswitch:
				imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/488Visconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))
				imp.show()  
				imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/561Visconv_" + "%03d" % i + "_" + "%02d" % j + "_0.tif"))
				imp.show()   
			IJ.run("Images to Stack", "name=Stack title=[] use");
			IJ.run("Select Bounding Box (guess background color)");
			IJ.run("Crop");
			IJ.run("Stack to Images");
			imp = IJ.getImage() 
			IJ.run(imp, "8-bit", "") 
			IJ.run(imp, "Canvas Size...", "width=2450 height=2450 position=Center zero");
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/prestitch/561Visconv_" + "%03d" % i + "_" + "%02d" % j + ".tif"));
			IJ.run(imp, "Enhance Contrast", "saturated=0.4 normalize");
			#IJ.run(imp, "Apply LUT","")
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/prestitch/561Visconv_" + "%03d" % i + "_" + "%02d" % j + "adj.tif"));
			imp.close();

			imp = IJ.getImage()  
			IJ.run(imp, "8-bit", "")
			IJ.run(imp, "Canvas Size...", "width=2450 height=2450 position=Center zero");
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/prestitch/488Visconv_" + "%03d" % i + "_" + "%02d" % j + ".tif"));
			imp.close();
        
			imp = IJ.getImage()  
			IJ.run(imp, "8-bit", "")
			IJ.run(imp, "Canvas Size...", "width=2450 height=2450 position=Center zero");
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/prestitch/647IRconv_" + "%03d" % i + "_" + "%02d" % j + ".tif"));
			imp.close();
        
			imp = IJ.getImage()
			IJ.run(imp, "8-bit", "")
			IJ.run(imp, "Canvas Size...", "width=2450 height=2450 position=Center zero");  
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/prestitch/750IRconv_" + "%03d" % i + "_" + "%02d" % j + ".tif"));
			imp.close();
        
			#imp = IJ.getImage()  
			#IJ.run(imp, "8-bit", "")
			#IJ.run(imp, "Canvas Size...", "width=2450 height=2450 position=Center zero");
			#IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/prestitch/561storm_" + "%02d" % i + "_" + "%02d" % j + ".tif"));
			#imp.close();

			#imp = IJ.getImage()  
			#IJ.run(imp, "8-bit", "")
			#IJ.run(imp, "Canvas Size...", "width=2450 height=2450 position=Center zero");
			#IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/prestitch/488storm_" + "%02d" % i + "_" + "%02d" % j + "_01.tif"));
			#imp.close();
        
			imp = IJ.getImage()  
			IJ.run(imp, "8-bit", "")
			IJ.run(imp, "Canvas Size...", "width=2450 height=2450 position=Center zero");
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/prestitch/488storm_" + "%03d" % i + "_" + "%02d" % j + ".tif"));
			imp.close();
        
			imp = IJ.getImage()  
			IJ.run(imp, "8-bit", "")
			IJ.run(imp, "Canvas Size...", "width=2450 height=2450 position=Center zero");
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/prestitch/647storm_" + "%03d" % i + "_" + "%02d" % j + ".tif"));
			imp.close();
        
			imp = IJ.getImage()  
			IJ.run(imp, "8-bit", "")
			IJ.run(imp, "Canvas Size...", "width=2450 height=2450 position=Center zero");
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/prestitch/750storm_" + "%03d" % i + "_" + "%02d" % j + ".tif"));
			imp.close();   
			IJ.run("Close All", ""); 
			print >>f1, "done with image resizing... starting tiling"
    	#parse out local_exp folder (exp_folder and coverslip number are needed) and make an unaligned_data folder
		leseg =  local_exp.split('/')
		lefolder =  str(leseg[-2:-1])[1:-1]
		coverslip = int(lefolder[-3:-1])
		print >>f1, coverslip
		exp_folder = '/'.join(leseg[:-3])
		print >>f1, exp_folder

		if not os.path.exists(exp_folder + "/unaligned/"):
			os.mkdir (exp_folder + "/unaligned/")
			os.mkdir (exp_folder + "/unaligned/conv_561/")
			os.mkdir (exp_folder + "/unaligned/conv_merged/")
			os.mkdir (exp_folder + "/unaligned/storm_merged/")
		if not os.path.exists(exp_folder + "/unaligned/conv_561_ds/"):
                        os.mkdir (exp_folder + "/unaligned/conv_561_ds/")
			os.mkdir (exp_folder + "/unaligned/conv_merged_ds/")
			os.mkdir (exp_folder + "/unaligned/storm_merged_ds/")
  
		

		for j in range (0, conv_images_per_section):
			cfold = str(int((math.floor(j/ydim))+1))
        	#save merged conv
			imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/prestitch/488Visconv_" + "%03d" % i + "_" + "%02d" % j + ".tif"))
			imp.show()
			IJ.run(imp, "8-bit", "");
			imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/prestitch/647IRconv_" + "%03d" % i + "_" + "%02d" % j + ".tif"))
			imp.show()
			IJ.run(imp, "8-bit", "");
			imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/prestitch/750IRconv_" + "%03d" % i + "_" + "%02d" % j + ".tif"))
			imp.show()
			IJ.run(imp, "8-bit", "");
			IJ.run("Merge Channels...", "red=750IRconv_" + "%03d" % i + "_" + "%02d" % j + ".tif green=647IRconv_" + "%03d" % i + "_" + "%02d" % j + ".tif blue=488Visconv_" + "%03d" % i + "_" + "%02d" % j + ".tif gray=*None*");
			imp = IJ.getImage() 
			IJ.saveAs("Tiff", (exp_folder + "/unaligned/conv_merged/" + "%02d" % coverslip + "%03d" % i + ".tif"))
                        print >>f1, "saving big convimage"
                        imp = IJ.getImage()  
			IJ.run(imp, "Size...", "width=" + str(imp.width/10) + " height=" + str(imp.height/10) + " constrain average interpolation=Bilinear");					
			IJ.saveAs("Tiff", (exp_folder + "/unaligned/conv_merged_ds/" + "%02d" % coverslip + "%03d" % i + ".tif"))
			print >>f1, "saving small convimage"
			IJ.run("Close All", ""); 
	
        	#save conv_561
			imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/prestitch/561Visconv_" + "%03d" % i + "_" + "%02d" % j + ".tif"))
			imp.show()
			IJ.run(imp, "8-bit", "");
			IJ.saveAs("Tiff", (exp_folder + "/unaligned/conv_561/" + "%02d" % coverslip + "%03d" % i + ".tif"))
                        print >>f1, "saving big conv561image"
                        imp = IJ.getImage()  
			IJ.run(imp, "Size...", "width=" + str(imp.width/10) + " height=" + str(imp.height/10) + " constrain average interpolation=Bilinear");					
			IJ.saveAs("Tiff", (exp_folder + "/unaligned/conv_561_ds/" + "%02d" % coverslip + "%03d" % i + ".tif"))
			print >>f1, "saving small conv561image"
			IJ.run("Close All", ""); 
		

        	#save merged storm
			imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/prestitch/488storm_" + "%03d" % i + "_" + "%02d" % j + ".tif"))
			imp.show()
			IJ.run(imp, "8-bit", "");
			imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/prestitch/647storm_" + "%03d" % i + "_" + "%02d" % j + ".tif"))
			imp.show()
			IJ.run(imp, "8-bit", "");
			imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/prestitch/750storm_" + "%03d" % i + "_" + "%02d" % j + ".tif"))
			imp.show()
			IJ.run(imp, "8-bit", "");
			IJ.run("Merge Channels...", "red=750storm_" + "%03d" % i + "_" + "%02d" % j + ".tif green=647storm_" + "%03d" % i + "_" + "%02d" % j + ".tif blue=488storm_" + "%03d" % i + "_" + "%02d" % j + ".tif gray=*None*");
			imp = IJ.getImage() 
			IJ.saveAs("Tiff", (exp_folder + "/unaligned/storm_merged/" + "%02d" % coverslip + "%03d" % i + ".tif"))
                        print >>f1, "saving big stormimage"
                        imp = IJ.getImage()  
			IJ.run(imp, "Size...", "width=" + str(imp.width/10) + " height=" + str(imp.height/10) + " constrain average interpolation=Bilinear");					
			IJ.saveAs("Tiff", (exp_folder + "/unaligned/storm_merged_ds/" + "%02d" % coverslip + "%03d" % i + ".tif"))
                        print >>f1, "saving small stormimage"
                        while not(os.path.exists(exp_folder + "/unaligned/storm_merged_ds/" + "%02d" % coverslip + "%03d" % i + ".tif")):
                                print "why no image???"
                                time.sleep(60)
			IJ.run("Close All", "");
time.sleep(60)
IJ.run("Quit")			 
