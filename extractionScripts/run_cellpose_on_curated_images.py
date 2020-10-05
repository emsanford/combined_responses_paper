import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import skimage.io
from cellpose import models, io, plot
import re
import glob
import os
import scipy

class Nucleus:
    """
    lightweight class to represent the data for a given nucleus type
    """
    def __init__(self, rep_name, condition, treatment_duration, antibody_stained, 
                 image_number, mask_number, intensity_value_array,
                 background_intensity_statistic, nuclear_annulus_intensity_statistic):
        self.rep_name = rep_name # "IFT1", "IFT2", etc...
        self.condition = condition  # EtOH, RA, TGFb, or Both
        self.treatment_duration = treatment_duration
        self.antibody_stained = antibody_stained
        self.image_number = image_number
        self.mask_number = mask_number
        self.intensity_value_array = intensity_value_array
        self.mean_intensity = np.mean(intensity_value_array)
        self.median_intensity = np.median(intensity_value_array)
        self.perc90_intensity = np.percentile(intensity_value_array, 90)
        self.perc75_intensity = np.percentile(intensity_value_array, 75)
        self.perc25_intensity = np.percentile(intensity_value_array, 25)
        self.perc10_intensity = np.percentile(intensity_value_array, 10)
        self.background_intensity_statistic = background_intensity_statistic
        self.nuclear_annulus_intensity_statistic = nuclear_annulus_intensity_statistic
        self.normalized_intensity = self.mean_intensity - self.background_intensity_statistic
        self.annulus_subtracted_intensity = self.mean_intensity - self.nuclear_annulus_intensity_statistic
        self.nuc_cytoplasm_ratio = (self.mean_intensity - self.background_intensity_statistic) / (self.nuclear_annulus_intensity_statistic - self.background_intensity_statistic)
        #self.normalized_intensity = self.mean_intensity


def calculate_frac_above_stdev_threshold(list_of_values, n_stdevs):
    """
    calculates fraction of values greater than n_stdevs standard deviations from the median
    """
    value_array = np.array(list_of_values)
    stdev  = np.std(value_array)
    median_val = np.median(value_array)
    n_values = len(value_array)
    frac_above_std_threshold = sum(value_array > (median_val + n_stdevs * stdev)) / float(n_values)
    return(frac_above_std_threshold)

def calculate_frac_above_median_plus_raw_threshold(list_of_values, threshold_interval):
    """
    calculates fraction of values greater than n_stdevs standard deviations from the median
    """
    value_array = np.array(list_of_values)
    median_val = np.median(value_array)
    n_values = len(value_array)
    frac_above_median_plus_threshold = sum(value_array > (median_val + threshold_interval)) / float(n_values)
    return(frac_above_median_plus_threshold)


grand_nucleus_list = []

os.environ['KMP_DUPLICATE_LIB_OK']='True'

# model_type='cyto' or model_type='nuclei'
model = models.Cellpose(model_type='nuclei')

# below define a list of DAPI files to run cellpose on. the folders containing these files should contain 
# similarly named tiff files for the other channels (e.g. cy, alexa, etc.)

files_rep1 = glob.glob("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/*/*/dapi*.tif")
#files_rep2 = glob.glob("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/IF_experiments_for_reviewers/IFT2_pSMAD2_RARA_6-Sep-2020/IFT2_analysis/curated_selection_of_images/*/*/dapi*.tif")
files_rep3 = glob.glob("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/*/*/dapi*.tif")


files =  files_rep3 + files_rep1

# define CHANNELS to run segementation on
# grayscale=0, R=1, G=2, B=3
# channels = [cytoplasm, nucleus]
# if NUCLEUS channel does not exist, set the second channel to 0
channels = [[0,0]]
# IF ALL YOUR IMAGES ARE THE SAME TYPE, you can give a list with 2 elements
# channels = [0,0] # IF YOU HAVE GRAYSCALE
# channels = [2,3] # IF YOU HAVE G=cytoplasm and B=nucleus
# channels = [2,1] # IF YOU HAVE G=cytoplasm and R=nucleus
    

for filename in files:
    print("current image file: {0}".format(filename))
    img = skimage.io.imread(filename)
    
    # if diameter is set to None, the size of the cells is estimated on a per image basis
    # you can set the average cell `diameter` in pixels yourself (recommended)
    # diameter can be a list or a single number for all images
    
    maskfile = filename.replace(".tif", "_cp_masks.png")
    # running cellpose takes the longest time. don't run it again if not necessary!
    if os.path.exists(maskfile):
        masks = skimage.io.imread(maskfile)
    else:
        masks, flows, styles, diams = model.eval(img, diameter=32, channels=channels)
        io.save_to_png(img, masks, flows, filename)
        
    # now make a nucleus for each mask in the maskfile. we will read in the 
    # cy image and collect the intensity values at the nuclei detected in 
    # the dapi image
    antibody_signal_image = skimage.io.imread(filename.replace("dapi", "cy"))
    
    
    regex_for_important_metadata = ".*(IFT[0-9]+)_(.*)_.*timepoint/(.*)_images_(.*)_stained_(.*)/dapi(.*).tif"
    match_object = re.match(regex_for_important_metadata, filename) 
    rep_name  = match_object.group(1)
    condition = match_object.group(3)
    antibody  = match_object.group(4)
    timepoint = match_object.group(5)
    image_num = match_object.group(6)
    
    number_of_nuclei_this_image = np.max(masks)
    
    array_of_background_and_cytoplasm_intensity_values = antibody_signal_image[masks == 0]
    ##percentile
    #background_intensity_statistic = np.percentile(array_of_background_and_cytoplasm_intensity_values, 95)
    #mode
    hist_step_size = 25
    numpy_hist_bg_vals = np.histogram(array_of_background_and_cytoplasm_intensity_values, bins = range(0, 8000, hist_step_size))
    logical_index_max_hist_val = numpy_hist_bg_vals[0] == np.max(numpy_hist_bg_vals[0])
    this_img_mode = numpy_hist_bg_vals[1][:-1][logical_index_max_hist_val] + hist_step_size / 2
    background_intensity_statistic = this_img_mode[0]
    
    # assume annulus around nucleus = cytoplasm, avg value of that
    strel_disk = np.array([[0,0,1,1,1,1,1,0,0],
                           [0,1,1,1,1,1,1,1,0],
                           [1,1,1,1,1,1,1,1,1],
                           [1,1,1,1,1,1,1,1,1],
                           [1,1,1,1,1,1,1,1,1],
                           [1,1,1,1,1,1,1,1,1],
                           [1,1,1,1,1,1,1,1,1],
                           [0,1,1,1,1,1,1,1,0],
                           [0,0,1,1,1,1,1,0,0]])

    for ii in range(1, number_of_nuclei_this_image + 1):
        ring_around_nuc = scipy.ndimage.morphology.binary_dilation(masks == ii, iterations = 2, structure = strel_disk) ^ (masks == ii)
        array_of_juxta_nuclear_intensity_values = antibody_signal_image[ring_around_nuc]
        nuclear_annulus_intensity_statistic = np.mean(array_of_juxta_nuclear_intensity_values)
        
        array_of_intensity_values = antibody_signal_image[masks == ii]
        if len(array_of_intensity_values) == 0:
            print("some mask had no actual intensity values assigned to it")
            continue
        this_nucleus = Nucleus(rep_name, condition, timepoint, antibody, image_num, ii, array_of_intensity_values, 
                               background_intensity_statistic, nuclear_annulus_intensity_statistic)
        grand_nucleus_list.append(this_nucleus)
        
    print("{0} {1} {2} bg intensity avg = {3:0.2f}".format(condition, antibody, timepoint, np.mean(array_of_background_and_cytoplasm_intensity_values)))
    # plt.hist(array_of_background_and_cytoplasm_intensity_values, bins = 200)
    # plt.xlim(0, 10000)
    # plt.title("{0} {1} {2} bg intensity stat = {3:0.2f}".format(condition, antibody, timepoint, background_intensity_statistic))
    # plt.show()
    # fig = plt.figure(figsize=(12,5))
    # plot.show_segmentation(fig, img, masks, flows[0], channels=channels[0])
    # plt.tight_layout()
    # plt.show()
range_list = []
for ii in range(-50, 500):
    range_list.append(ii / 100)


def plot_nucleus_histogram(nuclear_intensities, replicate_name, time_point_name, condition_name, histogram_bins, plot_descriptor):
    plt.hist(nuclear_intensities, bins = histogram_bins)
    plt.title(replicate_name + " " + condition_name + " " + time_point_name + " " + plot_descriptor)
    plt.xlim(-1000, 8000)
    plt.show()


output_filepath = "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/curated_image_nuclei_IF_statistics.tsv"
fout = open(output_filepath, "w")
fout.write('\t'.join(['replicate_ID', 'timepoint', 'condition', 'antibody_stained', 'image_number',
                      'annulus_subtracted_intensity', 'raw_mean_intensity']) + '\n')
for nuc in grand_nucleus_list:
    fout.write('\t'.join(map(str, [nuc.rep_name, nuc.treatment_duration, nuc.condition, 
                          nuc.antibody_stained, nuc.image_number,
                          nuc.annulus_subtracted_intensity, 
                          nuc.mean_intensity])) + '\n')
fout.close()
    
    


# plot_descriptor = "\nannulus-subtracted nuclear intensities"
# histbins = range(-1000, 8000, hist_step_size)
# outlier_removal_lower_boundary = -1000
# outlier_removal_upper_boundary = 8000
# #timepoint_strings = ["40min", "2hr", "4hr", "12hr", "24hr", "72hr"]
# timepoint_strings = ["40min", "2hr", "4hr", "12hr", "24hr", "72hr"]

# for replicate_name in ["IFT1", "IFT3"]:
#     nuclei_subset_1 = [x for x in grand_nucleus_list if x.rep_name == replicate_name]
    
#     for antibody in ["pSMAD2", "RARA"]:
#         nuclei_subset_2 = [x for x in nuclei_subset_1 if x.antibody_stained == antibody]
#         replicate_stat_list = []
#         for time_selection_string in timepoint_strings:
#             # 40 minute timepoint failed in replicate IFT3 due to very high and variable control well IF signal
#             if (time_selection_string == "40min" and replicate_name == "IFT3") or (antibody == "RARA"):
#                 continue
            
#             nuclei_subset_3 = [x for x in nuclei_subset_2 if x.treatment_duration == time_selection_string]
#             for condition_name in ["control", "ra_treated",  "tgfb_treated"]:
#                 nuclei_subset_4 = [x for x in nuclei_subset_3 if x.condition == condition_name]
                

#                 nuclear_intensities = [x.annulus_subtracted_intensity for x in nuclei_subset_4]
#                 plot_nucleus_histogram(nuclear_intensities, replicate_name, time_selection_string, condition_name, histbins, plot_descriptor)
                
#                 nstd_thresh = 2
#                 raw_intensity_thresh = 1000
                
#                 nuclear_intensities = np.array(nuclear_intensities)
#                 nuclear_intensities = nuclear_intensities[np.logical_and(nuclear_intensities > outlier_removal_lower_boundary, nuclear_intensities < outlier_removal_upper_boundary)]
                
#                 frac_above_std_thresh = calculate_frac_above_stdev_threshold(nuclear_intensities, nstd_thresh)
#                 frac_above_median_plus_fixed_thresh = calculate_frac_above_median_plus_raw_threshold(nuclear_intensities, raw_intensity_thresh)
#                 print("{0}\t{1}\t{2}   \tmedian: {3:.01f}\tstdev: {4:.01f}\tpct > 2 stdev: {5:.03f}\tmax: {6:.01f}\tmin: {7:.01f}\t{8:.3f}".format(replicate_name, time_selection_string,
#                                                                                               condition_name, np.median(nuclear_intensities),
#                                                                                               np.std(nuclear_intensities), frac_above_std_thresh,
#                                                                                               np.max(nuclear_intensities),
#                                                                                               np.min(nuclear_intensities),
#                                                                                               frac_above_median_plus_fixed_thresh))
#                 if antibody == 'pSMAD2':
#                     replicate_stat_list.append(frac_above_median_plus_fixed_thresh)
#                 elif antibody == 'RARA':
#                     replicate_stat_list.append(np.mean(nuclear_intensities))
                    
        
#         replicate_stat_list = np.array(replicate_stat_list)
#         n_stat_values  = len(replicate_stat_list)
#         control_stat_values = replicate_stat_list[range(0, n_stat_values, 3)]
#         ra_stat_values = replicate_stat_list[range(1, n_stat_values, 3)]
#         tgfb_stat_values = replicate_stat_list[range(2, n_stat_values, 3)]
#         x_axis_vals = range(1, 1 + len(timepoint_strings))
#         plt.plot(x_axis_vals, control_stat_values, "k")
#         plt.plot(x_axis_vals, ra_stat_values, "r")
#         plt.plot(x_axis_vals, tgfb_stat_values, "b")
#         if antibody == 'pSMAD2':
#             plt.title(replicate_name + " " + antibody + " fraction nuclei > median intensity + 1000")
#         elif antibody == 'RARA':
#             plt.title(replicate_name + " " + antibody + " mean nuclear intensity")
        
        
#         ctrl_patch = mpatches.Patch(color='k', label='etoh control')
#         ra_patch = mpatches.Patch(color='r', label='retinoic acid')
#         tgfb_patch = mpatches.Patch(color='b', label='tgf beta')
#         plt.legend(handles=[ctrl_patch, ra_patch, tgfb_patch])
        
#         plt.show()
    




















            
            