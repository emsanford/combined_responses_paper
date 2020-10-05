% This script is an example that makes image overlays. Will try and
% annotate all the relevant decisions and so forth.
% The paths below are based on the current working directory being the location of this script, 
% if it's not you can set the working directory by doing cd('~/Dropbox(Rajlab)/Projects/example_project/Paper/plotScripts/imageOverlays/'). 
writeDirectory_pSMAD2 = '../extractedData/immunofluorescence_image_data/replicate_1_cropped_images/pSMAD2_IF/';
writeDirectory_RARA   = '../extractedData/immunofluorescence_image_data/replicate_2_cropped_images/RARA_IF/';

cropped_image_width  = 511;
cropped_image_height = 511;

image_condition_order = [repelem("control", 6) repelem("ra_treated", 6) repelem("tgfb_treated", 6)];
image_timepoints      = ["40 min" "2 hr" "4 hr" "12 hr" "24 hr" "72  hr" "40 min" "2 hr" "4 hr" "12 hr" "24 hr" "72  hr" "40 min" "2 hr" "4 hr" "12 hr" "24 hr" "72  hr"];

% six paths for each selected image, ordered from early to later time point
time_ordered_paths_control      =  ["/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_40_min_timepoint/control_images_pSMAD2_stained_40min/cy026.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_2_hour_timepoint/control_images_pSMAD2_stained_2hr/cy023.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_4_hour_timepoint/control_images_pSMAD2_stained_4hr/cy021.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_12_hour_timepoint/control_images_pSMAD2_stained_12hr/cy017.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_24_hour_timepoint/control_images_pSMAD2_stained_24hr/cy008.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_72_hour_timepoint/control_images_pSMAD2_stained_72hr/cy025.tif"];
crop_starting_points_control    = [220 225;
                                   511 511;
                                   482 420;
                                     5  91;
                                   280  60;
                                   283 511];

time_ordered_paths_ra_treated   =  ["/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_40_min_timepoint/ra_treated_images_pSMAD2_stained_40min/cy020.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_2_hour_timepoint/ra_treated_images_pSMAD2_stained_2hr/cy009.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_4_hour_timepoint/ra_treated_images_pSMAD2_stained_4hr/cy003.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_12_hour_timepoint/ra_treated_images_pSMAD2_stained_12hr/cy027.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_24_hour_timepoint/ra_treated_images_pSMAD2_stained_24hr/cy028.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_72_hour_timepoint/ra_treated_images_pSMAD2_stained_72hr/cy026.tif"];
crop_starting_points_ra_treated = [511 511;
                                     1 511;
                                   129 251;
                                   381   4;
                                   285 115;
                                   497  80];

time_ordered_paths_tgfb_treated =  ["/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_40_min_timepoint/tgfb_treated_images_pSMAD2_stained_40min/cy030.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_2_hour_timepoint/tgfb_treated_images_pSMAD2_stained_2hr/cy016.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_4_hour_timepoint/tgfb_treated_images_pSMAD2_stained_4hr/cy021.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_12_hour_timepoint/tgfb_treated_images_pSMAD2_stained_12hr/cy010.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_24_hour_timepoint/tgfb_treated_images_pSMAD2_stained_24hr/cy009.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_1_curated_images/IFT1_72_hour_timepoint/tgfb_treated_images_pSMAD2_stained_72hr/cy018.tif"];
crop_starting_points_tgfb_treated = [503 410;
                                     497   1;
                                     245 106;
                                     500 152;
                                     395 366;
                                     326 133];

combined_image_filepaths = [time_ordered_paths_control;
                            time_ordered_paths_ra_treated;
                            time_ordered_paths_tgfb_treated];
combined_crop_starting_points = [crop_starting_points_control; 
                                 crop_starting_points_ra_treated;
                                 crop_starting_points_tgfb_treated];
                             
contrast_min = 1000;
contrast_max = 9000;

um_per_pixel = 665.60 /1024; % from 20x objective specs from IanD
number_um_in_scale_bar = 50;

                             
for ii = 1:length(image_condition_order)
    image_path     = combined_image_filepaths(ii);
    this_image     = imread(image_path);
    crop_rectangle = [combined_crop_starting_points(ii, 1) combined_crop_starting_points(ii, 2) cropped_image_width cropped_image_height];
    %disp(ii);
    %disp(crop_rectangle);
    cropped_image = imcrop(this_image, crop_rectangle);
    
    if ii == 1
        cropped_image = add_scale_bar(cropped_image, um_per_pixel, number_um_in_scale_bar);
    end
    
    split_image_path = split(image_path, '/');
    image_number_string = split_image_path(end);
    output_filepath = strjoin([writeDirectory_pSMAD2 image_condition_order(ii) '_' image_timepoints(ii) '_' image_number_string], "");
    imwrite(im2uint16(scale(cropped_image, [contrast_min contrast_max])), output_filepath);

end


%%% now do the same thing but for RARA images
time_ordered_paths_control      =  ["/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_40_min_timepoint/control_images_RARA_stained_40min/cy027.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_2_hour_timepoint/control_images_RARA_stained_2hr/cy015.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_4_hour_timepoint/control_images_RARA_stained_4hr/cy019.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_12_hour_timepoint/control_images_RARA_stained_12hr/cy036.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_24_hour_timepoint/control_images_RARA_stained_24hr/cy035.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_72_hour_timepoint/control_images_RARA_stained_72hr/cy009.tif"];
crop_starting_points_control    = [424 242;
                                   511 484;
                                    50 300;
                                     1 511;
                                   001 001;
                                   332 391];

time_ordered_paths_ra_treated   =  ["/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_40_min_timepoint/ra_treated_images_RARA_stained_40min/cy018.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_2_hour_timepoint/ra_treated_images_RARA_stained_2hr/cy014.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_4_hour_timepoint/ra_treated_images_RARA_stained_4hr/cy031.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_12_hour_timepoint/ra_treated_images_RARA_stained_12hr/cy035.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_24_hour_timepoint/ra_treated_images_RARA_stained_24hr/cy023.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_72_hour_timepoint/ra_treated_images_RARA_stained_72hr/cy022.tif"];
crop_starting_points_ra_treated = [423 511;
                                   250 422;
                                   351 511;
                                   213 266;
                                   112  33;
                                   167 269];

time_ordered_paths_tgfb_treated =  ["/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_40_min_timepoint/tgfb_treated_images_RARA_stained_40min/cy027.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_2_hour_timepoint/tgfb_treated_images_RARA_stained_2hr/cy016.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_4_hour_timepoint/tgfb_treated_images_RARA_stained_4hr/cy009.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_12_hour_timepoint/tgfb_treated_images_RARA_stained_12hr/cy033.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_24_hour_timepoint/tgfb_treated_images_RARA_stained_24hr/cy025.tif";
                                    "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/immunofluorescence_image_data/replicate_2_curated_images/IFT3_72_hour_timepoint/tgfb_treated_images_RARA_stained_72hr/cy025.tif"];
crop_starting_points_tgfb_treated = [438 79;
                                     186 202;
                                     290 511;
                                     001 001;
                                       4 343;
                                       9 448];

combined_image_filepaths = [time_ordered_paths_control;
                            time_ordered_paths_ra_treated;
                            time_ordered_paths_tgfb_treated];
combined_crop_starting_points = [crop_starting_points_control; 
                                 crop_starting_points_ra_treated;
                                 crop_starting_points_tgfb_treated];
                             
contrast_min = 1000;
contrast_max = 20000;
                             
for ii = 1:length(image_condition_order)
    image_path     = combined_image_filepaths(ii);
    this_image     = imread(image_path);
    crop_rectangle = [combined_crop_starting_points(ii, 1) combined_crop_starting_points(ii, 2) cropped_image_width cropped_image_height];
    %disp(ii);
    %disp(crop_rectangle);
    cropped_image = imcrop(this_image, crop_rectangle);
    
    if ii == 1
        cropped_image = add_scale_bar(cropped_image, um_per_pixel, number_um_in_scale_bar);
    end
    
    split_image_path = split(image_path, '/');
    image_number_string = split_image_path(end);
    output_filepath = strjoin([writeDirectory_RARA image_condition_order(ii) '_' image_timepoints(ii) '_' image_number_string], "");
    imwrite(im2uint16(scale(cropped_image, [contrast_min contrast_max])), output_filepath);

end


function image_scale_bar = add_scale_bar(image, um_per_pixel, width_um)

    width_pixels = width_um / um_per_pixel;
    height_pixels = 8;
    
    offset = 20;
    start_x = offset;
    start_y = size(image, 1) - 10 - height_pixels;
    
    image_scale_bar = insertShape(image, 'FilledRectangle', [start_x, start_y, width_pixels, height_pixels], 'Color', 'white', 'Opacity', 1);

end
