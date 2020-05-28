import os,sys
import glob
import shutil
from ij import IJ  
from ij.process import ImageStatistics as IS
from ij import Prefs

argv = getArgument()
print str(argv)
arg_array = argv.split(" ")
print str(arg_array)

local_exp = arg_array[0]

acq_folder = local_exp + 'acquisition/'
ffc_folder = acq_folder

#get rid of any output files before we begin to avoid error message
if not os.path.exists(acq_folder + "dist_corr/output/"):
                os.mkdir (acq_folder + "dist_corr/output/")
subfiles = glob.glob(acq_folder + "dist_corr/output/*")
for s in subfiles:
        os.remove(s)



#I647_2560
distcorr_src = ffc_folder + "dist_corr/distCorr.txt"
IJ.run("lenscorrection Distortion Correction", "calibration=" + ffc_folder + 
                "dist_corr/ number_of_images=9 first_image=647ffc_00.tif " +
                "power_of_polynomial_kernel=3 lambda=1.0 apply_correction_to_images what=save " +
        "file_name=distCorr.txt initial_gaussian_blur=1.6 steps_per_scale_octave=4 " +
        "minimum_image_size=64 maximum_image_size=500 feature_descriptor_size=8 " +
        "feature_descriptor_orientation_bins=8 closest/next_closest_ratio=0.92 " +
        "maximal_alignment_error=32 inlier_ratio=0.20 expected_transformation=Rigid target=" +
        ffc_folder + "dist_corr/output/");
IJ.run("Quit")			 

