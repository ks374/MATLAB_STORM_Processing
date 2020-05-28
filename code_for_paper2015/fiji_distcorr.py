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
                   		strsequence[2:] + "_00_0.tif"):
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
				imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/488storm_" + "%03d" % i + "_" + "%02d" % j + "_1.tif"))
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
				imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/647storm_" + "%03d" % i + "_" + "%02d" % j + "_1.tif"))
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
				imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/after_dist_corr/750storm_" + "%03d" % i + "_" + "%02d" % j + "_1.tif"))
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
			
			imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/prestitch/561Visconv_" + "%03d" % i + "_" + "%02d" % j + "adj.tif"))
			imp.show()
			IJ.run(imp, "8-bit", "");
			image = WindowManager.getCurrentImage()
			imp = IJ.getImage() 
			cfold = str(int((math.floor(j/ydim))+1))
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/conv_561adj/c_" + cfold + "/561Visconv_" + "%03d" % i + "_" + "%02d" % j + ".tif"));
			IJ.run("Close All", ""); 

			target_dir = ISanalysisfolder + "%04d" % i + "/aligned/conv_561adj/out/"
			oldprocs = (glob.glob(target_dir + "*" + ".tif"))
			for op in oldprocs:
				print >>f1, "removing old files"
				os.remove(op)
			t2 = ISanalysisfolder + "%04d" % i + "/aligned/conv_561adj/mc/"
			oldprocs2 = (glob.glob(t2 + "*" + ".tif"))
			for op2 in oldprocs2:
				print >>f1, "removing old mc files"
				os.remove(op2)
		for c in range(1,(xdim+1)):
			cs= str(c)
                           # source directory
			source_dir = ISanalysisfolder + "%04d" % i + "/aligned/conv_561adj/c_" + cs + "/"
                          # output directory
			
                          # transforms directory
			transf_dir = ISanalysisfolder + "%04d" % i + "/aligned/conv_561adj/c_" + cs + "/"
                          # reference image
			file = (glob.glob(source_dir + "*" + ".tif"))
			print >>f1, os.path.basename(file[0])
			reference_name = os.path.basename(file[0])
			# shrinkage option (false)
			use_shrinking_constraint = 0
       
			p = Register_Virtual_Stack_MT.Param()
			p.featuresModelIndex = 0
			p.registrationModelIndex = 0
			p.interpolate = 1
			p.maxEpsilon = 25
			# The "maximum image size":
			p.sift.maxOctaveSize = 400
			p.sift.fdBins = 8
			p.sift.fdSize = 8
			p.sift.initialSigma = 1.0
			p.sift.steps = 6
			p.sift.minOctaveSize  = 64
			# The "inlier ratio":
			p.minInlierRatio = 0.03

               		iterat = 0
                        def signal_handler(signum, frame):
    				raise Exception("Reg_fail!")
			signal.signal(signal.SIGALRM, signal_handler)
			signal.alarm(1200)  
			while iterat<5: 
                                Register_Virtual_Stack_MT.exec(source_dir, target_dir, transf_dir, reference_name, p, use_shrinking_constraint)
                                src_files = glob.glob(source_dir + '*.tif')
                                file1 =  src_files[0].split('/')
                                filelast = src_files[-1].split('/')
                                file1 =  file1[-1]
                                filelast = filelast[-1]
                                print >>f1, file1,filelast
                                if os.path.exists(target_dir + filelast):
                                        print >>f1, "filelast exists"
                                        lastsize = os.path.getsize(target_dir + filelast)
                                        firstsize = os.path.getsize(target_dir + file1)
                                        if firstsize == lastsize:
                                                iterat = 5
                                                print >>f1, "register virtual stack success"
                                        else:
                                                print "files are different sizes"
                                                iterat = iterat+1
                                                print >>f1, "trying register virtual stack again"                   
                                else:
                                        iterat = iterat+1
                                        print >>f1, "trying register virtual stack again"

                        imp = IJ.getImage()  
                        IJ.run(imp, "Z Project...", "start=1 stop=5 projection=[Max Intensity]");
                        imp = IJ.getImage()  
                        IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/conv_561adj/mc/c" + cs + ".tif"));
                        IJ.run("Close All", "");
                        signal.alarm(0)	 					
                        print >>f1, "done with first register virtual stack"
                

		if not xdim==1:  
       		# source directory
			source_dir = ISanalysisfolder + "%04d" % i + "/aligned/conv_561adj/mc/"
       		# output directory
			target_dir = ISanalysisfolder + "%04d" % i + "/aligned/conv_561adj/out/"
       		# transforms directory
			transf_dir = ISanalysisfolder + "%04d" % i + "/aligned/conv_561adj/mc/"
       		# reference image
			file = (glob.glob(source_dir + "*" + ".tif"))
			print >>f1, os.path.basename(file[0])
			reference_name = os.path.basename(file[0])
               		iterat = 0
                        def signal_handler(signum, frame):
    				raise Exception("Reg_fail!")
			signal.signal(signal.SIGALRM, signal_handler)
			signal.alarm(1200)  
			while iterat<5: 
                                Register_Virtual_Stack_MT.exec(source_dir, target_dir, transf_dir, reference_name, p, use_shrinking_constraint)
                                src_files = glob.glob(source_dir + '*.tif')
                                file1 =  src_files[0].split('/')
                                filelast = src_files[-1].split('/')
                                file1 =  file1[-1]
                                filelast = filelast[-1]
                                print >>f1, file1,filelast
                                if os.path.exists(target_dir + filelast):
                                        print  >>f1,"filelast exists"
                                        lastsize = os.path.getsize(target_dir + filelast)
                                        firstsize = os.path.getsize(target_dir + file1)
                                        if firstsize == lastsize:
                                                iterat = 5
                                                print >>f1, "register virtual stack success"
                                        else:
                                                print "files are different sizes"
                                                iterat = iterat+1
                                                print >>f1, "trying register virtual stack again"                   
                                else:
                                        iterat = iterat+1
                                        print >>f1, "trying register virtual stack again"

                        imp = IJ.getImage()  
                        IJ.run(imp, "Z Project...", "start=1 stop=5 projection=[Max Intensity]");
                        imp = IJ.getImage()  
                        IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/conv_561adj/mc/c" + cs + ".tif"));
                        IJ.run("Close All", "");
                        signal.alarm(0)	 					
                        print >>f1, "done with first register virtual stack"
                

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
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/conv_merged/c_" + cfold + "/conv_merged" + "%03d" % i + "_" + "%02d" % j + ".tif"));

			#save storm_561
			#imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/prestitch/561storm_" + "%02d" % i + "_" + "%02d" % j + ".tif"))
			#imp.show()
			#IJ.run(imp, "8-bit", "");
			#IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/storm_561/c_" + cfold + "/561storm" + "%02d" % i + "_" + "%02d" % j + ".tif"));

        	#save storm_488_2nd
			#imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/prestitch/488storm_" + "%02d" % i + "_" + "%02d" % j + "_01.tif"))
			#imp.show()
			#IJ.run(imp, "8-bit", "");
			#IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/storm_488_2nd/c_" + cfold + "/488storm" + "%02d" % i + "_" + "%02d" % j + ".tif"));
	
        	#save conv_561
			imp = IJ.openImage((ISanalysisfolder + "%04d" % i + "/aligned/prestitch/561Visconv_" + "%03d" % i + "_" + "%02d" % j + ".tif"))
			imp.show()
			IJ.run(imp, "8-bit", "");
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/conv_561/c_" + cfold + "/561Visconv" + "%03d" % i + "_" + "%02d" % j + ".tif"));


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
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/storm_merged/c_" + cfold + "/storm_merged" + "%03d" % i + "_" + "%02d" % j + ".tif"));
	
		for c in range(1,(xdim+1)):
			cs= str(c)
			# source directory
			source_dir = ISanalysisfolder + "%04d" % i + "/aligned/conv_merged/c_" + cs + "/"
			# output directory
			target_dir = ISanalysisfolder + "%04d" % i + "/aligned/conv_merged/out/"
			# transforms directory
			transf_dir = ISanalysisfolder + "%04d" % i + "/aligned/conv_561adj/c_" + cs + "/"
        	      
			IJ.run("Transform Virtual Stack Slices", "source=" + source_dir + " output=" + target_dir + " transforms=" + transf_dir + " interpolate");
			imp = IJ.getImage()  
			IJ.run(imp, "Z Project...", "start=1 stop=5 projection=[Max Intensity]");
			imp = IJ.getImage()  
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/conv_merged/mc/c" + cs + ".tif"));
			if xdim==1:
				print >>f1, "saving big convimage"
				IJ.saveAs("Tiff", (exp_folder + "/unaligned/conv_merged/" + "%02d" % coverslip + "%03d" % i + ".tif"))
                                imp = IJ.getImage()  
				IJ.run(imp, "Size...", "width=" + str(imp.width/10) + " height=" + str(imp.height/10) + " constrain average interpolation=Bilinear");					
				print >>f1, "saving little conv image"
				IJ.saveAs("Tiff", (exp_folder + "/unaligned/conv_merged_ds/" + "%02d" % coverslip + "%03d" % i + ".tif"))
                                
			IJ.run("Close All", ""); 
			
			#	# source directory
			#	source_dir = ISanalysisfolder + "%04d" % i + "/aligned/storm_561/c_" + cs + "/"
			#	# output directory
			#	target_dir = ISanalysisfolder + "%04d" % i + "/aligned/storm_561/out/" 
			#	IJ.run("Transform Virtual Stack Slices", "source=" + source_dir + " output=" + target_dir + " transforms=" + transf_dir + " interpolate");
			#	imp = IJ.getImage()  
			#	IJ.run(imp, "Z Project...", "start=1 stop=5 projection=[Max Intensity]");
			#	imp = IJ.getImage()  
			#	IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/storm_561/mc/c" + cs + ".tif"));
			#	IJ.run("Close All", ""); 

			#	source_dir = ISanalysisfolder + "%04d" % i + "/aligned/storm_488_2nd/c_" + cs + "/"
			#	# output directory
			#	target_dir = ISanalysisfolder + "%04d" % i + "/aligned/storm_488_2nd/out/"        
			#	IJ.run("Transform Virtual Stack Slices", "source=" + source_dir + " output=" + target_dir + " transforms=" + transf_dir + " interpolate");
			#	imp = IJ.getImage()  
			#	IJ.run(imp, "Z Project...", "start=1 stop=5 projection=[Max Intensity]");
			#	imp = IJ.getImage()  
			#	IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/storm_488_2nd/mc/c" + cs + ".tif"));
			#	IJ.run("Close All", ""); 

			source_dir = ISanalysisfolder + "%04d" % i + "/aligned/conv_561/c_" + cs + "/"
			# output directory
			target_dir = ISanalysisfolder + "%04d" % i + "/aligned/conv_561/out/"
			IJ.run("Transform Virtual Stack Slices", "source=" + source_dir + " output=" + target_dir + " transforms=" + transf_dir + " interpolate");
			imp = IJ.getImage()  
			IJ.run(imp, "Z Project...", "start=1 stop=5 projection=[Max Intensity]");
			imp = IJ.getImage()  
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/conv_561/mc/c" + cs + ".tif"));
			print >>f1, "saving big conv561image"
			if xdim==1:
				IJ.saveAs("Tiff", (exp_folder + "/unaligned/conv_561/" + "%02d" % coverslip + "%03d" % i + ".tif"))
                                imp = IJ.getImage()  
				IJ.run(imp, "Size...", "width=" + str(imp.width/10) + " height=" + str(imp.height/10) + " constrain average interpolation=Bilinear");					
				IJ.saveAs("Tiff", (exp_folder + "/unaligned/conv_561_ds/" + "%02d" % coverslip + "%03d" % i + ".tif"))
				print >>f1, "saving small561 convimage"
			IJ.run("Close All", "");

			source_dir = ISanalysisfolder + "%04d" % i + "/aligned/storm_merged/c_" + cs + "/"
			# output directory
			target_dir = ISanalysisfolder + "%04d" % i + "/aligned/storm_merged/out/"
			IJ.run("Transform Virtual Stack Slices", "source=" + source_dir + " output=" + target_dir + " transforms=" + transf_dir + " interpolate");
			imp = IJ.getImage()  
			IJ.run(imp, "Z Project...", "start=1 stop=5 projection=[Max Intensity]");
			imp = IJ.getImage()  
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/storm_merged/mc/c" + cs + ".tif"));
			if xdim==1:
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
      	if not xdim==1:   
			# source directory
			source_dir = ISanalysisfolder + "%04d" % i + "/aligned/conv_merged/mc/"
			# output directory
			target_dir = ISanalysisfolder + "%04d" % i + "/aligned/conv_merged/out/"
			# transforms directory
			transf_dir = ISanalysisfolder + "%04d" % i + "/aligned/conv_561adj/mc/"
           
			IJ.run("Transform Virtual Stack Slices", "source=" + source_dir + " output=" + target_dir + " transforms=" + transf_dir + " interpolate");
			imp = IJ.getImage()  
			IJ.run(imp, "Z Project...", "start=1 stop=5 projection=[Max Intensity]");
			imp = IJ.getImage()  
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/conv_merged_" + "%04d" % i + ".tif"));
			IJ.saveAs("Tiff", (exp_folder + "/unaligned/conv_merged/" + "%02d" % coverslip + "%03d" % i + ".tif"))
                        print >>f1, "saving big convimage"
                        imp = IJ.getImage()  
			IJ.run(imp, "Size...", "width=" + str(imp.width/10) + " height=" + str(imp.height/10) + " constrain average interpolation=Bilinear");					
			IJ.saveAs("Tiff", (exp_folder + "/unaligned/conv_merged_ds/" + "%02d" % coverslip + "%03d" % i + ".tif"))
			print >>f1, "saving small convimage"
			IJ.run("Close All", ""); 
	
			## source directory
			#source_dir = ISanalysisfolder + "%04d" % i + "/aligned/storm_561/mc/"
			## output directory
			#target_dir = ISanalysisfolder + "%04d" % i + "/aligned/storm_561/out/"          
			#IJ.run("Transform Virtual Stack Slices", "source=" + source_dir + " output=" + target_dir + " transforms=" + transf_dir + " interpolate");
			#imp = IJ.getImage()  
			#IJ.run(imp, "Z Project...", "start=1 stop=5 projection=[Max Intensity]");
			#imp = IJ.getImage()  
			#IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/storm_561_" + "%04d" % i + ".tif"));
			#IJ.saveAs("Tiff", (rc_store + rc_exp + "storm_561_unaligned/" + "%02d" % coverslip + "%02d" % i + ".tif"))
			#IJ.run("Close All", ""); 
		
			## source directory
			#source_dir = ISanalysisfolder + "%04d" % i + "/aligned/storm_488_2nd/mc/"
			## output directory
			#target_dir = ISanalysisfolder + "%04d" % i + "/aligned/storm_488_2nd/out/"            
			#IJ.run("Transform Virtual Stack Slices", "source=" + source_dir + " output=" + target_dir + " transforms=" + transf_dir + " interpolate");
			#imp = IJ.getImage()  
			#IJ.run(imp, "Z Project...", "start=1 stop=5 projection=[Max Intensity]");
			#imp = IJ.getImage()  
			#IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/storm_488_2nd_" + "%04d" % i + ".tif"));
			#IJ.saveAs("Tiff", (rc_store + rc_exp + "storm_488_2nd_unaligned/" + "%02d" % coverslip + "%02d" % i + ".tif"))
			#IJ.run("Close All", ""); 
		
			# source directory
			source_dir = ISanalysisfolder + "%04d" % i + "/aligned/conv_561/mc/"
			# output directory
			target_dir = ISanalysisfolder + "%04d" % i + "/aligned/conv_561/out/"        
			IJ.run("Transform Virtual Stack Slices", "source=" + source_dir + " output=" + target_dir + " transforms=" + transf_dir + " interpolate");
			imp = IJ.getImage()  
			IJ.run(imp, "Z Project...", "start=1 stop=5 projection=[Max Intensity]");
			imp = IJ.getImage()  
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/conv_561_" + "%04d" % i + ".tif"));
			IJ.saveAs("Tiff", (exp_folder + "/unaligned/conv_561/" + "%02d" % coverslip + "%03d" % i + ".tif"))
                        print >>f1, "saving big conv561image"
                        imp = IJ.getImage()  
			IJ.run(imp, "Size...", "width=" + str(imp.width/10) + " height=" + str(imp.height/10) + " constrain average interpolation=Bilinear");					
			IJ.saveAs("Tiff", (exp_folder + "/unaligned/conv_561_ds/" + "%02d" % coverslip + "%03d" % i + ".tif"))
			print >>f1, "saving small conv561image"
			IJ.run("Close All", ""); 
		
			# source directory
			source_dir = ISanalysisfolder + "%04d" % i + "/aligned/storm_merged/mc/"
			# output directory
			target_dir = ISanalysisfolder + "%04d" % i + "/aligned/storm_merged/out/"
			IJ.run("Transform Virtual Stack Slices", "source=" + source_dir + " output=" + target_dir + " transforms=" + transf_dir + " interpolate");
			imp = IJ.getImage()  
			IJ.run(imp, "Z Project...", "start=1 stop=5 projection=[Max Intensity]");
			imp = IJ.getImage()  
			IJ.saveAs(imp, "Tiff", (ISanalysisfolder + "%04d" % i + "/aligned/storm_merged_" + "%04d" % i + ".tif"));
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
