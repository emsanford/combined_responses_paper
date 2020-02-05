import os
import glob
import re

class ParamSet:
	# a parameter set will define the parameters we're testing in the sensitivity analysis
	# these include:
	#    x: the size for which to merge neighboring differential peaks together
	#    x: the minimum normalized coverage for allowing a peak to be considered differential
	#    x: the minimum fold change for allowing a peak to be considered differential
	#    x: 
	def __init__(self, base_directory, peak_merge_distance = 250, peak_width = 150, diffpeak_algorithm_min_normalized_fragments = 10, 
				 diffpeak_algorithm_min_fold_change = 1.1, min_control_TPM = 0, fixed_vs_variable_width = "FixedWidth"):
		self.base_directory      						 = base_directory           # this should be the analysis root directory, which contains folders like "extractedData" and "plotScripts"
		self.peak_merge_distance 						 = peak_merge_distance
		self.peak_width          				 		 = peak_width
		self.fixed_vs_variable_width 					 = fixed_vs_variable_width
		self.diffpeak_algorithm_min_normalized_fragments = diffpeak_algorithm_min_normalized_fragments
		self.diffpeak_algorithm_min_fold_change          = diffpeak_algorithm_min_fold_change
		self.min_control_TPM         					 = min_control_TPM
		param_summary_string = "mergedist{0}_peakwidth{1}_minNormFrags{2}_minFoldChange{3}".format(peak_merge_distance, 
																								   peak_width,
																								   diffpeak_algorithm_min_normalized_fragments, 
																								   diffpeak_algorithm_min_fold_change)
		# input fastq files (paired end)
		self.unmerged_differential_atac_peaks_bed_file = '"{0}"'.format(os.sep.join([base_directory, "extractedData", "differentialAtacPeaks_forIGV.bed"]))
		self.merged_differential_atac_peaks_file = '"{0}"'.format(os.sep.join([base_directory, "plots", "sensitivity_analysis", "differentialAtacPeaks_mergedist{0}.bed".format(peak_merge_distance)]))
		self.diffpeak_algorithm_output_file      = '"{0}"'.format(os.sep.join([base_directory, "plots", "sensitivity_analysis", "differentialAtacPeaks_{0}.tsv".format(param_summary_string)]))
		self.annotated_diffpeaks_output_file     = '"{0}"'.format(os.sep.join([base_directory, "plots", "sensitivity_analysis", "differentialAtacPeaks_{0}.annotated.tsv".format(param_summary_string)]))
		self.upregulated_diffpeaks_output_file   = '"{0}"'.format(os.sep.join([base_directory, "plots", "sensitivity_analysis", "differentialAtacPeaks_{0}.annotated.upregulated.tsv".format(param_summary_string)]))
		self.peak_fdr_grid_plot_path_prefix      = '"{0}"'.format(os.sep.join([base_directory, "plots", "sensitivity_analysis", "output_plots", "fdr_grid_{0}".format(param_summary_string)]))
		self.integration_histogram_path_prefix   = '"{0}"'.format(os.sep.join([base_directory, "plots", "sensitivity_analysis", "output_plots", "cval_histogram_{0}".format(param_summary_string)]))
		self.upreg_peak_cats_bed_file_prefix     = '"{0}"'.format(os.sep.join([base_directory, "plots", "sensitivity_analysis", "upregulated_peaks_{0}".format(param_summary_string)]))
		self.joined_gene_and_peak_table_file     = '"{0}"'.format(os.sep.join([base_directory, "plots", "sensitivity_analysis", "upregJoinedPeakGeneTib_{0}.tsv".format(param_summary_string)]))
		self.peaks_near_genes_plots_path_prefix  = '"{0}"'.format(os.sep.join([base_directory, "plots", "sensitivity_analysis", "output_plots", "peaks_near_genes_{0}".format(param_summary_string)]))


		self.path_to_createFragCountMatx = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "createAtacFragmentCountMatrix.R"]))
		self.atac_frag_count_rds_file    = '"{0}"'.format(os.sep.join([base_directory, "plots", "sensitivity_analysis", "atac_fragment_counts_mergedist{0}.rds".format(peak_merge_distance)]))

		# list of filepaths, we need to add quotes around them to pass as command line arguments because dropbox added a space to its own folder and we're using os.system to run commands
		self.path_to_fdrBasedDiffPeakCalling = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "runEmpiricalFdrBasedDifferentialPeakSelection.R"]))
		self.path_to_peak_annotation_script  = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "addIntegrationMetricsAndMotifMatchesToDiffPeaks.R"]))
		self.path_to_create_upregulated_peaks_script  = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "createMasterSetOfUpregulatedPeaks.R"]))
		self.path_to_peak_integration_category_histograms_script = '"{0}"'.format(os.sep.join([base_directory, "plotScripts", "peakIntegrationSummaryPieChartsAndHistograms.R"]))
		self.path_to_make_bed_files_for_each_category = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "createBedFilesForPeakIntegrationCategories.R"]))
		self.path_to_chromVAR_motif_analysis_script = '"{0}"'.format(os.sep.join([base_directory, "plotScripts", "chromVarIndSignalMotifEnrichments.R"]))
		self.path_to_join_peak_to_gene_tib = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "joinNearbyPeaksToGenes.R"]))
		self.path_to_make_peak_near_gene_analysis_plots = '"{0}"'.format(os.sep.join([base_directory, "plotScripts", "makeAdjacentSuperadditivePeakAssocModeOfIntegrationPlots.R"]))


	def __str__(self):
		obj_string = "peak_merge_distance = {0}\n" \
					 "peak_width = {1}\n" \
					 "diffpeak_algorithm_min_normalized_fragments = {2}\n" \
					 "diffpeak_algorithm_min_fold_change = {3}\n" \
					 "fixed_vs_variable_width = {4}".format(self.peak_merge_distance, 
					 										self.peak_width, 
					 										self.diffpeak_algorithm_min_normalized_fragments, 
					 										self.diffpeak_algorithm_min_fold_change,
					 										self.fixed_vs_variable_width)
		return(obj_string)

def run_command(cmd, run_scripts_without_arguments):
	if run_scripts_without_arguments:
		match_obj = re.match('(Rscript \".*?\").*', cmd)
		cmd2 = match_obj.group(1)
		print(cmd2)
		os.system(cmd2)
	else:
		print(cmd)
		os.system(cmd)


def main(param_obj, run_scripts_without_arguments = False):
	# make a merged peak file with the specified merge distance (note, we could in theory use the master peak list for this rather than the original FDR 0.01 peak list)
	if run_scripts_without_arguments or not os.path.exists(param_obj.merged_differential_atac_peaks_file[1:-1]):
		if run_scripts_without_arguments:
			merge_infile           = param_obj.unmerged_differential_atac_peaks_bed_file
			merge_outfile          = '"{0}"'.format(os.sep.join([param_obj.base_directory, "extractedData", "differentialAtacPeaks_merged_forIGV.bed"]))
			merge_distance_default = 250
			cmd = "bedtools merge -d {0} -i {1} > {2}".format(param_obj.peak_merge_distance, merge_infile, merge_outfile)
			print(cmd)
			os.system(cmd)
		else:
			cmd = "bedtools merge -d {0} -i {1} > {2}".format(param_obj.peak_merge_distance, param_obj.unmerged_differential_atac_peaks_bed_file, param_obj.merged_differential_atac_peaks_file)
			run_command(cmd, run_scripts_without_arguments)

	# make a new fragment count matrix for the specified peak file
	if run_scripts_without_arguments or not os.path.exists(param_obj.atac_frag_count_rds_file[1:-1]):
		cmd = "Rscript {0} {1} {2} {3} {4}".format(param_obj.path_to_createFragCountMatx,
										   param_obj.merged_differential_atac_peaks_file,
										   param_obj.peak_width, 
										   param_obj.atac_frag_count_rds_file, 
										   "VariableWidth")

		run_command(cmd, run_scripts_without_arguments)


	#option for long scripts:
		# if not os.path.exists(output file)
	# run a new FDR-based differential peak analysis on the new peak set 
	if run_scripts_without_arguments or not os.path.exists(param_obj.diffpeak_algorithm_output_file[1:-1]):
		cmd = "Rscript {0} {1} {2} {3} {4} {5}".format(param_obj.path_to_fdrBasedDiffPeakCalling, 
													   param_obj.diffpeak_algorithm_min_normalized_fragments, 
													   param_obj.diffpeak_algorithm_min_fold_change,
													   param_obj.atac_frag_count_rds_file, 
													   param_obj.diffpeak_algorithm_output_file,
													   param_obj.peak_fdr_grid_plot_path_prefix)
		run_command(cmd, run_scripts_without_arguments)

	# annotate the new peak set with motif matches and additive / multiplicative predictions
	if run_scripts_without_arguments or not os.path.exists(param_obj.annotated_diffpeaks_output_file[1:-1]):
		cmd = "Rscript {0} {1} {2}".format(param_obj.path_to_peak_annotation_script, 
										   param_obj.diffpeak_algorithm_output_file, 
										   param_obj.annotated_diffpeaks_output_file)
		run_command(cmd, run_scripts_without_arguments)

	# create the set of upregulated peaks from the annotated peak set
	if run_scripts_without_arguments or not os.path.exists(param_obj.upregulated_diffpeaks_output_file[1:-1]):
		cmd = "Rscript {0} {1} {2}".format(param_obj.path_to_create_upregulated_peaks_script, 
										   param_obj.annotated_diffpeaks_output_file, 
										   param_obj.upregulated_diffpeaks_output_file)
		run_command(cmd, run_scripts_without_arguments)

	# make the peak signal integration histogram and pie charts for this parameter set
	if run_scripts_without_arguments or not os.path.exists(param_obj.integration_histogram_path_prefix[1:-1] + "_med_dose.svg"):
		cmd = "Rscript {0} {1} {2}".format(param_obj.path_to_peak_integration_category_histograms_script, 
										   param_obj.upregulated_diffpeaks_output_file, 
										   param_obj.integration_histogram_path_prefix)
		run_command(cmd, run_scripts_without_arguments)

	# make bed files for sub-additive, additive, and super-additive peaks
	upregulated_peaks_files = glob.glob(param_obj.upreg_peak_cats_bed_file_prefix[1:-1] + "*.bed")
	if run_scripts_without_arguments or len(upregulated_peaks_files) == 0:
		cmd = "Rscript {0} {1} {2}".format(param_obj.path_to_make_bed_files_for_each_category, 
										   param_obj.upregulated_diffpeaks_output_file, 
										   param_obj.upreg_peak_cats_bed_file_prefix)
		run_command(cmd, run_scripts_without_arguments)

	# make new fragment counts files for each peakset
	upregulated_peaks_fragCount_r_objects = glob.glob(param_obj.upreg_peak_cats_bed_file_prefix[1:-1] + "*.rds")
	if run_scripts_without_arguments or len(upregulated_peaks_fragCount_r_objects) == 0:
		for upreg_peak_file in upregulated_peaks_files:
			output_filename = upreg_peak_file + ".rds"
			cmd = "Rscript {0} {1} {2} {3} {4}".format(param_obj.path_to_createFragCountMatx,
											   '"{0}"'.format(upreg_peak_file),
											   param_obj.peak_width, 
											   '"{0}"'.format(output_filename), 
											   "VariableWidth")
			run_command(cmd, run_scripts_without_arguments)

	# use chromVAR to do motif analysis at the different categories of upregulated peaks
	# outputPlotPrefix, "raw_dev_score_by_tf_name.svg"
	if run_scripts_without_arguments or not os.path.exists(param_obj.atac_frag_count_rds_file[1:-1] + "_motifPlots_raw_dev_score_by_tf_name.svg"):
		upregulated_peaks_fragCount_r_objects = glob.glob(param_obj.upreg_peak_cats_bed_file_prefix[1:-1] + "*.rds")
		upregulated_peaks_fragCount_r_objects.append(param_obj.atac_frag_count_rds_file[1:-1])
		for peak_r_object in upregulated_peaks_fragCount_r_objects:
			cmd = "Rscript {0} {1} {2}".format(param_obj.path_to_chromVAR_motif_analysis_script, 
										   '"{0}"'.format(peak_r_object), 
										   '"{0}"'.format(peak_r_object + "_motifPlots_"))
			run_command(cmd, run_scripts_without_arguments)

	# are the super-additive peaks more likely to be close to super-multiplicative genes?
	path_to_upregulated_gene_table = '"{0}"'.format(r"/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv")
	# make new joined peak tib
	if run_scripts_without_arguments or not os.path.exists(param_obj.joined_gene_and_peak_table_file[1:-1]):
		cmd = "Rscript {0} {1} {2} {3}".format(param_obj.path_to_join_peak_to_gene_tib,
											   path_to_upregulated_gene_table,
											   param_obj.upregulated_diffpeaks_output_file,
											   param_obj.joined_gene_and_peak_table_file)
		run_command(cmd, run_scripts_without_arguments)
	# use the new joined peak tib to make the "special peak types near genes" plots
	peaks_near_genes_analysis_plots = glob.glob(param_obj.peaks_near_genes_plots_path_prefix[1:-1] + "*.svg")
	if run_scripts_without_arguments or len(peaks_near_genes_analysis_plots) == 0:
		cmd = "Rscript {0} {1} {2} {3} {4} {5}".format(param_obj.path_to_make_peak_near_gene_analysis_plots,
													   path_to_upregulated_gene_table,
													   param_obj.joined_gene_and_peak_table_file,
													   param_obj.min_control_TPM,
													   param_obj.peaks_near_genes_plots_path_prefix + "superadditivePeaks_",
													   "superadditive")
		run_command(cmd, run_scripts_without_arguments)

		cmd = "Rscript {0} {1} {2} {3} {4} {5}".format(param_obj.path_to_make_peak_near_gene_analysis_plots,
													   path_to_upregulated_gene_table,
													   param_obj.joined_gene_and_peak_table_file,
													   param_obj.min_control_TPM,
													   param_obj.peaks_near_genes_plots_path_prefix + "additivePeaks_",
													   "additive")
		run_command(cmd, run_scripts_without_arguments)

		cmd = "Rscript {0} {1} {2} {3} {4} {5}".format(param_obj.path_to_make_peak_near_gene_analysis_plots,
													   path_to_upregulated_gene_table,
													   param_obj.joined_gene_and_peak_table_file,
													   param_obj.min_control_TPM,
													   param_obj.peaks_near_genes_plots_path_prefix + "subadditivePeaks_",
													   "subadditive")
		run_command(cmd, run_scripts_without_arguments)

	

if __name__ == '__main__':
	basedir = r"/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4"
	fixed_vs_variable_width = "VariableWidth"           # choose from "FixedWidth", "VariableWidth"
	peak_merge_distances = [0, 250, 375, 500]           # also try [0, 150, 250, 500]
	peak_width           = 150                          # also try 225 later, for motif analysis
	diffpeak_algorithm_min_normalized_fragments_list = [10, 30]   # do 10, 20, 30, 40
	diffpeak_algorithm_min_fold_change_list          = [1.1, 1.5]  # do 1,1, 1.5, 2, 4
	param_object_list = []

	for pmd in peak_merge_distances:
		for ii in [0, 1]:
			param_object_list.append(ParamSet(basedir, 
											  peak_merge_distance = pmd, 
											  peak_width = 150,
											  diffpeak_algorithm_min_normalized_fragments = diffpeak_algorithm_min_normalized_fragments_list[ii],
											  diffpeak_algorithm_min_fold_change = diffpeak_algorithm_min_fold_change_list[ii], 
											  fixed_vs_variable_width = fixed_vs_variable_width))

	param_obj = ParamSet(basedir, 
						 peak_merge_distance = 250, 
						 peak_width = 150,
						 diffpeak_algorithm_min_normalized_fragments = 30,
						 diffpeak_algorithm_min_fold_change = 1.5, 
						 fixed_vs_variable_width = fixed_vs_variable_width,
						 min_control_TPM = 1)

	param_object_list = [param_obj]  # uncomment this line to do a parameter sweep
	run_scripts_without_arguments = True

	for param_obj in param_object_list:
		if not run_scripts_without_arguments:
			print("Running parameter set:\n{0}".format(str(param_obj)))
		main(param_obj, run_scripts_without_arguments = run_scripts_without_arguments)


