This folder contains code used in alignment and analysis of data in Sigal et al. Cell, 2015
Analysis is for data collected using STORM-control software Hal4000,Steve,and Dave found at github/ZhuangLab
Analysis requires STORM-analysis software compiled also found at github/ZhuangLab along with the python dependencies described for the STORM-analysis software.  Fiji and Matlab are also required, and additional preexisting matlab functions (provided) should be put on the matlab path.
Parts of the code are written for  job submission for SLURM at Harvard Research computing.
This code has three main parts with brief descriptions as to their function.
1)XY single section alignment
2)Z serial section alignment
3)Neuron and Synapses Analysis





XY single section alignment
uss_align.py - when called from command line, this program launches a gui to guide parameters for single section alignment.  This gui does minimal compute and can be run from login nodes.  For each coverslip run, one bead fitting job will be submitted.  For each physical x-y section run, one molecule fitting and one image processing job will be submitted. Run from folder uss_align.py is in.
Input variables are:
	Storage Drive - path to sever and folder where data is generally stored
	Home Drive - path to home directory
	Experiment Name - 
	Coverslip start/end - if multiple coverslips were acquired, chose those to analyze. Pushing the "Find Coverslips" button will identify those coverslips to be analyzed.  A value of -1 in coverslip end will select all coverslips
	Small/Large x-y dimension - If more than one tile per section is acquired, set the shape (ie. 1x4 or 2x2)
	Redo?	-  Expects three digits either 0 or 1.  If partial or complete analysis for a section is already done, sets if you want to redo.  First digit is bead fitting, second digit is storm molecule list generation, third is image processing.  
	Align all tiles? 0 to run only the first section for each coverslip, 1 to run all sections
	Which Scope - differences are based on the camera used and current choices are "STORM2", "STORM2U" for the ultra camera, and "STORM4" 
	WGA-488? either 0 or 1. switches alignment channel is WGA-488 was used
within code, make note of/change path for partitions, mufit,xml variables, pyimproc variables.
job checking no longer works as currently written with new implementation of squeue function at harvard RC as of 6/2015.

cproc.py - creates folders with dates to store job output and error files for job submission.  
rel_conv_ints.py - determine conventional image intensities of data aquired to use for image normalization
recon_bead.py - launched from uss_align.  Runs bead alignment and generates transforms for chromatic abberation and lens distortion. generates flat field correction images for conventional images.  Should be in same folder as uss_align.py.
add_bead_images.py - sums several frames of bead images to increase signal useful for some channels
image_proc_all.py - converts .dax into tif images, .bin into tif images and generates and applies flat field correction
mufit_analysis.py - interfaces with daostorm fitting software
gen_bead_warp.m - generates chromatic abberation warping from set of fit beads. should be in same folder as uss_align.py
test_bead_warp.m - tests the above generated transform on a second set of beads. should be in same folder as uss_align.py
fiji_distcorr_make.py - generates the lens distortion correction transform.  should be in the macros folder of the fiji install.
recon_tile1.py - launches from uss_align.  performs daostorm molecule fitting for all storm movies for a section.should be in same folder as uss_align.py
recon_tile2.py - launches from uss_align.  processes dax and bin files to images. applies image corrections and tiles overlapping sections. should be in same folder as uss_align.py
image_chrom_align.m -  Applies chromatic abberation correction to conventional and storm images. should be in same folder as uss_align.py
image_chrom_align_visonly.m -  same as above except used if not 750 channel was collected. should be in same folder as uss_align.py
fiji_distcorr.py - Applies distortion correction, and uses SIFT to align overlapping tiles. should be located in fiji macros folder
fiji_distcorr1field.py - same as above, but for when only a single tile per section is acquired.  should be located in fiji macros folder
register_virtual_stack and lenscorrection folders -  these folders contain only small changes to exsisting fiji plugins to suppress desired manual input. should be located in fiji plugins folder







Z serial section alignment
uss_z_align.py - this program is run from the command line and requires prior manual editing of the code for each dataset .  First manually inspect the sections to check for out of focus sections/corupted sections/out of order sections! It calls several sequential scripts to align the serial sections. These are:  1)Copies original data before reordering/removing sections in case an error is made. 2)Reorders sections if necessary and removes "bad" sections. 3) Normalizes WGA signal and creates merged WGA/neuron image for alignment. 4) Performs rigid alignment 5) Rotate and crop the dataset to the desired size. ** this only opens a matlab terminal! you must run crop_datastack.m and manually determine angle of rotation and cropping** 5) elastic alignment of data 6) resave out elastic alignment images without regenerating alignment in case program freezes
 It does not submit jobs but assumes it is being run from a compute node. Usually the requested memory is ~30GB for large datasets.  do not run from a login node.
	parameters to change:
		-Set storage,home,and experiment paths as in to uss_align above.
		- if sections were collected out of order or contain bad sections, change out_of_order variable to 1.  Also provide the correct ordering in starts,storms,and excludes.  Each pair of starts/stops in the list sets the sequential number of images for the proper ordering.

reorder.py - corrects ordering of sections with input provided in uss_z_align
wga_norm_and_thresh.m - normalizes wga intensity across images and combines with neuron signal for alignment.  The neuron signal is included at a higher intensity.  **This is definitely better for rigid alignment, but I'm unsure if it is best for elastic align.  May be better to only use wga signal for elastic align** . should be in same folder as uss_z_align.py
fiji_z_align1.py - performs rigid alignment of dataset.  takes a while to create stack before input is required. Still requires manual input of sift parameters.  those parameters used can be found at the bottom of the uss_z_align script, but may vary from dataset to dataset.  these parameters can be tested using the "Feature Extraction/SIFT" plugin in Fiji. should be in fiji macro folder
crop_datastack.m - as the image stack can become much larger after the rigid alignment and can contain a large amount of blank space, this allows for the rotation and cropping of the dataset to minimize processing time.  angle of rotation and cropping region are determined manually.
fiji_z_align2.py - performs elastic alignment of dataset.   takes a while to create stack before input is required.  Still requires manual input of parameters.  Select elastic registration from pulldown menu. and set first and last sections to the first and last section of the data.  Also necessary to select that the data is pre-aligned. those parameters used can be found at the bottom of the uss_z_align script, but may vary from dataset to dataset.  these parameters can be tested using the "Registration/Elastic/test elastic parameters " plugin in Fiji.   should be in fiji macro folder
fiji_z_align3.py - finishes saving out aligned data if fiji_z_align2 doesn't finish saving images.   should be in fiji macro folder




numel(find(statsGwater(i,1).PixelValues2>(1.2*threshfactorg(1))))
numel(find(statsGwater(i,1).PixelValues(statsGwater(i,1).PixelValues2>threshfactorg(1))>0))



Neuron and Synapses Analysis
Disclaimer- As a whole these scripts and the parameters used are less set and optimized.  The code used for the the large ds cell is being provided here while some changes were made for the receptor datasets.  In its current state, code was run in sections sometimes with manual input.  While the results are robust across the datasets and signals analyzed in this paper, it will not be surprising if these codes require optimization for other types of synaptic signals In addition, some values provided within these scripts may still be dataset specific.  Further improvement of these codes represents ongoing work, so please contact me with questions.
For large datasets, these scripts require a large amound of memory ~250GB for the large ds cell. 
initial whole dataset processing
 It is also important to run the rot_data script first followed by a norm_thresh script before proceeding to further analysis of synapse and neuron signals.  
rot_data_neuronsig.m - this is used to rotate the data into the coordinate system of the IPL.  For whole cells, this is well determined based on the stratification of the neuron.  Also, downsampled images are made which contain only the IPL data, this is used for future signal normalization and synapse identification .
rot_data_manual.m - for incomplete cells that were imaged in cross section, determining the angle of rotation based on the stratification of the cell did not always work.  Instead, here it is possible to manually set the x,y,and z angles to rotate the data into the proper orientation.  
norm_thresh_smoothed- this is used to normalize the signal intensity of conventional and storm images across sections and apply the low pass conventional image filter to the storm data.  the smoothed version is being used for data aquired in the plane of the IPL to all for a rolling histogram across sections.
norm_thresh_mean - same as above, but for data aquired in crosssection of the IPL, assumes a constant histogram
neuron signal processing
parseYFP.m - finds neuron signal and saves out images of selected signal.  Allows for selection of more or less unconnected components of neuron based both on size of signal and distance from main neuron branch
ysurfdist.py - as determining the distance matrix for a large object such as the neuron is too time and memory intensive to do at once, this parses smaller pieces and submits multiple jobs to the cluster for this function.  calls parseyfp_surfdist.m
synapse signal processing
segG.m and segP.m - for gephyrin and presynaptic signals respectively, this is used to find a list of clusters within the signal
synapse_select_G.m and synapse_select_P.m - these are run after the segG+P scripts and are used to separate synaptic from non synaptic populations.
add_to_statsG.m and add_to_statsP.m - these are run after the synapse selection scripts and are used to query the opposing synaptic signal surrounding each synaptic cluster

pairing and on_neuron analysis
synapse_pairing.m - after above synapse signals have been segmented, selected, and the shell analysis has been performed, this is used to generate the centroid nearest neighbor pairing and the two parameter pairing analysis.
geph_on_neuron.m - after all the above analysis is completed, this script is used to query the distance of gephyrin clusters to the surface of the neuron and identify synapses associated with the neuron.
