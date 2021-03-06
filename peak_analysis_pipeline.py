import os
import glob
import re
import time
from pyprojroot import here


class ParamSet:
	def __init__(self, base_directory, extractedDataDir, plotsDir, peak_merge_distance = 250, initial_peak_width = 150, 
				 initial_diffpeak_algorithm_min_normalized_fragments = 10, initial_diffpeak_algorithm_min_fold_change = 1.1, 
				 final_diffpeak_algorithm_min_normalized_fragments = 30, final_diffpeak_algorithm_min_fold_change = 1.5,
				 addPredFcDiffMin_integration_histogram = 0.0, minNormFragCtDiff_integration_histogram = 0,
				 use_default_param_string = False):
		self.base_directory      						 	     = base_directory  # this should be the analysis root directory, which contains folders like "extractedData" and "plotScripts"
		self.sample_metadata_file								 = '"{0}"'.format(base_directory + os.sep + "sampleMetadata_SI2-SI4.txt")
		self.extractedDataDir 							 	     = extractedDataDir
		self.plotsDir       							 	     = plotsDir
		self.peak_merge_distance 						 		 = peak_merge_distance
		self.initial_peak_width          				 		 = initial_peak_width
		self.initial_diffpeak_algorithm_min_normalized_fragments = initial_diffpeak_algorithm_min_normalized_fragments
		self.initial_diffpeak_algorithm_min_fold_change          = initial_diffpeak_algorithm_min_fold_change		
		self.final_diffpeak_algorithm_min_normalized_fragments   = final_diffpeak_algorithm_min_normalized_fragments
		self.final_diffpeak_algorithm_min_fold_change            = final_diffpeak_algorithm_min_fold_change
		self.addPredFcDiffMin_integration_histogram              = addPredFcDiffMin_integration_histogram
		self.minNormFragCtDiff_integration_histogram             = minNormFragCtDiff_integration_histogram
		param_summary_string = "mergedist{0}_peakwidth{1}_minNormFrags{2}_minFoldChange{3}".format(peak_merge_distance, 
																								   initial_peak_width,
																								   final_diffpeak_algorithm_min_normalized_fragments, 
																								   final_diffpeak_algorithm_min_fold_change)
		if use_default_param_string:
			param_summary_string = "defaultParams"

		# paths to files that should already exist (from runnin gene analysis pipeline)
		self.path_to_annotated_gene_table   = '"{0}"'.format(extractedDataDir + os.sep + "DeSeqOutputAllConds.annotated.tsv")
		self.path_to_upregulated_gene_table = '"{0}"'.format(extractedDataDir + os.sep + "DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv")


		# paths to folders containing output files
		self.consensus_peak_file_dir                             = '"{0}"'.format(os.sep.join([extractedDataDir, "consensusPeakFiles"]))
		self.venn_diagrams_directory                             = os.sep.join([plotsDir, "venn_diagrams"])
		self.integration_summary_plots_dir                       = os.sep.join([plotsDir, "peak_integration_summary_plots"])
		self.integration_summary_null_distribution_plots_dir     = os.sep.join([plotsDir, "peak_integration_summary_plots", "null_distributions"])
		self.peaks_near_genes_plots_dir                          = os.sep.join([plotsDir, "peaks_near_gene_types_plots"]) 
		self.motif_analysis_plots_dir                            = os.sep.join([plotsDir, "motif_analysis_plots"]) 
		self.subdirectories 							 		 = [plotsDir, extractedDataDir, self.consensus_peak_file_dir[1:-1], 
																	self.integration_summary_plots_dir, self.integration_summary_null_distribution_plots_dir,
																	self.venn_diagrams_directory, self.peaks_near_genes_plots_dir, self.motif_analysis_plots_dir]

		# paths to input and output files
		self.merged_consensus_peak_file                         = '"{0}"'.format(os.sep.join([extractedDataDir, "mergedConsensusMacs2PeakFilesAllConditions.bed"]))
		self.unmerged_differential_atac_peaks_bed_file          = '"{0}"'.format(os.sep.join([base_directory, "extractedData", "initialDifferentialAtacPeaksWidth{0}_minNlCovg{1}_minFc{2}.bed".format(initial_peak_width, initial_diffpeak_algorithm_min_normalized_fragments, initial_diffpeak_algorithm_min_fold_change)]))
		self.merged_differential_atac_peaks_bed_file            = '"{0}"'.format(os.sep.join([extractedDataDir, "differentialAtacPeaks_mergedist{0}.bed".format(peak_merge_distance)]))
		self.final_differential_atac_peaks_bed_file             = '"{0}"'.format(os.sep.join([extractedDataDir, "differentialAtacPeaks_final_{0}.bed".format(param_summary_string)]))
		self.atac_fragment_count_file_merged_consensus_peak_set = '"{0}"'.format(os.sep.join([extractedDataDir, "mergedConsensusPeakSetAtacFragmentCounts.rds"]))
		self.merged_differential_atac_frag_count_rds_file       = '"{0}"'.format(os.sep.join([extractedDataDir, "merged_diffPeaks_fragment_counts_mergedist{0}.rds".format(peak_merge_distance)]))
		self.final_merged_differential_atac_frag_count_rds_file = '"{0}"'.format(os.sep.join([extractedDataDir, "final_diffPeaks_fragment_counts_{0}.rds".format(param_summary_string)]))
		self.initial_diffpeak_algorithm_output_file             = '"{0}"'.format(os.sep.join([extractedDataDir, "initialDifferentialAtacPeaksWidth{0}_minNlCovg{1}_minFc{2}.tsv".format(initial_peak_width, initial_diffpeak_algorithm_min_normalized_fragments, initial_diffpeak_algorithm_min_fold_change)]))
		self.final_diffpeak_algorithm_output_file               = '"{0}"'.format(os.sep.join([extractedDataDir, "differentialAtacPeaks_{0}.tsv".format(param_summary_string)]))
		self.most_variable_motifs_file                          = '"{0}"'.format(os.sep.join([extractedDataDir, "mostVariableMotifs_{0}.rds".format(param_summary_string)]))
		self.annotated_diffpeaks_output_file                    = '"{0}"'.format(os.sep.join([extractedDataDir, "differentialAtacPeaks_{0}.annotated.tsv".format(param_summary_string)]))
		self.upregulated_diffpeaks_output_file                  = '"{0}"'.format(os.sep.join([extractedDataDir, "differentialAtacPeaks_{0}.annotated.upregulated.tsv".format(param_summary_string)]))
		self.initial_peak_fdr_grid_plot_path_prefix             = '"{0}"'.format(os.sep.join([plotsDir, "initial_fdr_grid_{0}".format(param_summary_string)]))
		self.final_peak_fdr_grid_plot_path_prefix               = '"{0}"'.format(os.sep.join([plotsDir, "final_fdr_grid_{0}".format(param_summary_string)]))
		self.upreg_peak_cats_bed_file_prefix                    = '"{0}"'.format(os.sep.join([extractedDataDir, "upregulated_peaks_{0}".format(param_summary_string)]))
		self.joined_gene_and_upreg_peak_table_file              = '"{0}"'.format(os.sep.join([extractedDataDir, "upregGenesUpregPeaksJoinedTib_{0}.tsv".format(param_summary_string)]))
		self.joined_gene_and_all_peak_table_file                = '"{0}"'.format(os.sep.join([extractedDataDir, "upregGenesAllPeaksJoinedTib_{0}.tsv".format(param_summary_string)]))
		self.peaks_near_genes_plots_path_prefix                 = '"{0}"'.format(os.sep.join([plotsDir, "peaks_near_genes_{0}".format(param_summary_string)]))


		# list of filepaths to scripts, we need to add quotes around them to pass as command line arguments because dropbox added a space to its own folder and we're using os.system to run commands
		self.path_to_makeConsensusPeakFilesForEachCond           = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "makeConsensusPeakFilesForEachCondition.R"]))
		self.path_to_mergeConsensusPeaksOfSeveralConditions      = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "mergeConsensusPeaksOfSeveralConditions.R"]))
		self.path_to_createFragCountMatx                         = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "createAtacFragmentCountMatrix.R"]))
		self.path_to_fdrBasedDiffPeakCalling                     = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "runEmpiricalFdrBasedDifferentialPeakSelection.R"]))
		self.path_to_createBedFileFromDiffpeakTable              = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "createBedFileFromDiffpeakTable.R"]))
		self.path_to_makeMostVariableMotifSet                    = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "makeMostVariableMotifSet.R"]))
		self.path_to_peak_annotation_script                      = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "addIntegrationMetricsAndMotifMatchesToDiffPeaks.R"]))
		self.path_to_create_upregulated_peaks_script             = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "createMasterSetOfUpregulatedPeaks.R"]))
		self.path_to_peak_integration_category_histograms_script = '"{0}"'.format(os.sep.join([base_directory, "plotScripts", "peakIntegrationSummaryPieChartsAndHistograms.R"]))
		self.path_to_diff_peak_and_gene_venn_diagram_script      = '"{0}"'.format(os.sep.join([base_directory, "plotScripts", "make_DiffPeakAndDiffGeneVennDiagrams.R"]))
		self.path_to_makeNullDistributionCorDvalue               = '"{0}"'.format(os.sep.join([base_directory, "plotScripts", "makeNullDistributionCorDvalue.R"]))
		self.path_to_make_bed_files_for_each_category            = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "createBedFilesForPeakIntegrationCategories.R"]))
		self.path_to_join_peak_to_gene_tib                       = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "joinNearbyPeaksToGenes.R"]))
		self.path_to_make_peak_near_gene_analysis_plots          = '"{0}"'.format(os.sep.join([base_directory, "plotScripts", "makeAdjacentPeakTypeAssocToGeneTypeModeOfIntegrationPlots.R"]))
		self.path_to_motif_analysis_plots                        = '"{0}"'.format(os.sep.join([base_directory, "plotScripts", "makeMotifAnalysisPlots.R"]))
		self.path_to_motif_size_and_density_analysis_plots       = '"{0}"'.format(os.sep.join([base_directory, "plotScripts", "avgNumMotifsPerPeakType.R"]))
		self.path_to_dual_motif_match_plot_script                = '"{0}"'.format(os.sep.join([base_directory, "plotScripts", "obsVsExpectedDualMotifRateByPeakType.R"]))
		self.path_to_supplemental_motif_analysis_plots           = '"{0}"'.format(os.sep.join([base_directory, "plotScripts", "makeSupplementalMotifAnalysisPlots.R"]))



	def __str__(self):
		obj_string = "peak_merge_distance = {0}\n" \
					 "initial_peak_width = {1}\n" \
					 "final_diffpeak_algorithm_min_normalized_fragments = {2}\n" \
					 "final_diffpeak_algorithm_min_fold_change = {3}\n".format(self.peak_merge_distance, 
					 										self.initial_peak_width, 
					 										self.final_diffpeak_algorithm_min_normalized_fragments, 
					 										self.final_diffpeak_algorithm_min_fold_change)
		return(obj_string)

def run_command(cmd):
	time.sleep(1)  # but in this buffer to help avoid the problem where a new process tries to run before its output file is fully written
	print(cmd)
	os.system(cmd)


def main(param_obj, run_all_steps = False):
	# make sub-directories if they don't already exist
	for dirname in param_obj.subdirectories:
		if not os.path.exists(dirname):
			os.mkdir(dirname)

	# make consensus peaks file
	consensus_peak_files = glob.glob(param_obj.consensus_peak_file_dir[1:-1] + "/*.bed")
	if run_all_steps or len(consensus_peak_files) != 12:
		cmd = 'Rscript {0} {1} {2}'.format(param_obj.path_to_makeConsensusPeakFilesForEachCond,
										   param_obj.sample_metadata_file,
										   param_obj.consensus_peak_file_dir)
		run_command(cmd)

	# merge consensus peak files
	if run_all_steps or not os.path.exists(param_obj.merged_consensus_peak_file[1:-1]):
		cmd = 'Rscript {0} {1} {2}'.format(param_obj.path_to_mergeConsensusPeaksOfSeveralConditions,
										   param_obj.consensus_peak_file_dir,
										   param_obj.merged_consensus_peak_file)
		run_command(cmd)

	# make atac frag count file for consensus peak file
	if run_all_steps or not os.path.exists(param_obj.atac_fragment_count_file_merged_consensus_peak_set[1:-1]):
		cmd = "Rscript {0} {1} {2} {3} {4}".format(param_obj.path_to_createFragCountMatx,
										   param_obj.merged_consensus_peak_file,
										   param_obj.initial_peak_width, 
										   param_obj.atac_fragment_count_file_merged_consensus_peak_set, 
										   "FixedWidth")

		run_command(cmd)

	# run diff peak algorithm with initial FDR-based parameters
	if run_all_steps or not os.path.exists(param_obj.initial_diffpeak_algorithm_output_file[1:-1]):
		cmd = "Rscript {0} {1} {2} {3} {4} {5}".format(param_obj.path_to_fdrBasedDiffPeakCalling, 
													   param_obj.initial_diffpeak_algorithm_min_normalized_fragments, 
													   param_obj.initial_diffpeak_algorithm_min_fold_change,
													   param_obj.atac_fragment_count_file_merged_consensus_peak_set, 
													   param_obj.initial_diffpeak_algorithm_output_file,
													   param_obj.initial_peak_fdr_grid_plot_path_prefix)
		run_command(cmd)

	# make a bed file from the differential peak file
	if run_all_steps or not os.path.exists(param_obj.unmerged_differential_atac_peaks_bed_file[1:-1]):
		cmd = 'Rscript {0} {1} {2}'.format(param_obj.path_to_createBedFileFromDiffpeakTable,
										   param_obj.initial_diffpeak_algorithm_output_file,
										   param_obj.unmerged_differential_atac_peaks_bed_file)
		run_command(cmd)


	# make a merged peak file with the specified merge distance (note, we could in theory use the master peak list for this rather than the original FDR 0.01 peak list)
	if run_all_steps or not os.path.exists(param_obj.merged_differential_atac_peaks_bed_file[1:-1]):
		cmd = "bedtools merge -d {0} -i {1} > {2}".format(param_obj.peak_merge_distance, param_obj.unmerged_differential_atac_peaks_bed_file, param_obj.merged_differential_atac_peaks_bed_file)
		run_command(cmd)

	# make a new fragment count matrix for the merged initial differential peak file
	if run_all_steps or not os.path.exists(param_obj.merged_differential_atac_frag_count_rds_file[1:-1]):
		cmd = "Rscript {0} {1} {2} {3} {4}".format(param_obj.path_to_createFragCountMatx,
										   param_obj.merged_differential_atac_peaks_bed_file,
										   param_obj.initial_peak_width, 
										   param_obj.merged_differential_atac_frag_count_rds_file, 
										   "VariableWidth")

		run_command(cmd)

	# run a new FDR-based differential peak analysis on the new peak set 
	if run_all_steps or not os.path.exists(param_obj.final_diffpeak_algorithm_output_file[1:-1]):
		cmd = "Rscript {0} {1} {2} {3} {4} {5}".format(param_obj.path_to_fdrBasedDiffPeakCalling, 
													   param_obj.final_diffpeak_algorithm_min_normalized_fragments, 
													   param_obj.final_diffpeak_algorithm_min_fold_change,
													   param_obj.merged_differential_atac_frag_count_rds_file, 
													   param_obj.final_diffpeak_algorithm_output_file,
													   param_obj.final_peak_fdr_grid_plot_path_prefix)
		run_command(cmd)

	# make a bed file from the final differential peak file
	if run_all_steps or not os.path.exists(param_obj.final_differential_atac_peaks_bed_file[1:-1]):
		cmd = 'Rscript {0} {1} {2}'.format(param_obj.path_to_createBedFileFromDiffpeakTable,
										   param_obj.final_diffpeak_algorithm_output_file,
										   param_obj.final_differential_atac_peaks_bed_file)
		run_command(cmd)

	# create a final fragmentCount rds file for which to do motif analysis
	if run_all_steps or not os.path.exists(param_obj.final_merged_differential_atac_frag_count_rds_file[1:-1]):
		cmd = "Rscript {0} {1} {2} {3} {4}".format(param_obj.path_to_createFragCountMatx, 
												   param_obj.final_differential_atac_peaks_bed_file,
												   param_obj.initial_peak_width,
												   param_obj.final_merged_differential_atac_frag_count_rds_file,
												   "VariableWidth")
		run_command(cmd)

	# make the list of the most variable motifs in the differential peak set
	if run_all_steps or not os.path.exists(param_obj.most_variable_motifs_file[1:-1]):
		cmd = "Rscript {0} {1} {2}".format(param_obj.path_to_makeMostVariableMotifSet, 
										   param_obj.final_merged_differential_atac_frag_count_rds_file,
										   param_obj.most_variable_motifs_file)
		run_command(cmd)

	# annotate the new peak set with motif matches and additive / multiplicative predictions
	if run_all_steps or not os.path.exists(param_obj.annotated_diffpeaks_output_file[1:-1]):
		cmd = "Rscript {0} {1} {2} {3}".format(param_obj.path_to_peak_annotation_script, 
											   param_obj.final_diffpeak_algorithm_output_file,
											   param_obj.most_variable_motifs_file,
											   param_obj.annotated_diffpeaks_output_file)
		run_command(cmd)

	# create the set of upregulated peaks from the annotated peak set
	if run_all_steps or not os.path.exists(param_obj.upregulated_diffpeaks_output_file[1:-1]):
		cmd = "Rscript {0} {1} {2}".format(param_obj.path_to_create_upregulated_peaks_script, 
										   param_obj.annotated_diffpeaks_output_file, 
										   param_obj.upregulated_diffpeaks_output_file)
		run_command(cmd)

	# make the peak signal integration histogram and pie charts for this parameter set
	pie_chart_output_files = glob.glob(param_obj.integration_summary_plots_dir + os.sep + "*pie_chart*.svg")
	if run_all_steps or len(pie_chart_output_files) == 0:
		cmd = "Rscript {0} {1} {2}".format(param_obj.path_to_peak_integration_category_histograms_script, 
										   param_obj.upregulated_diffpeaks_output_file, 
										   '"{0}"'.format(param_obj.integration_summary_plots_dir))
		run_command(cmd)

	# make the venn diagrams of number of differential peaks and genes for each signal treatment
	venn_diagram_paths = glob.glob(param_obj.venn_diagrams_directory + os.sep + '*.svg')
	if run_all_steps or len(venn_diagram_paths) == 0:
		cmd = "Rscript {0} {1} {2} {3}".format(param_obj.path_to_diff_peak_and_gene_venn_diagram_script, 
											   param_obj.path_to_annotated_gene_table, 
											   param_obj.annotated_diffpeaks_output_file,
											   '"{0}{1}"'.format(param_obj.venn_diagrams_directory, os.sep))
		run_command(cmd)

	integration_summary_null_histogram_paths = glob.glob(param_obj.integration_summary_null_distribution_plots_dir + '/*.svg')
	if run_all_steps or len(integration_summary_null_histogram_paths) == 0:
		cmd = 'Rscript {0} {1} {2} {3} {4} {5}'.format(param_obj.path_to_makeNullDistributionCorDvalue,
													   param_obj.upregulated_diffpeaks_output_file,
													   param_obj.addPredFcDiffMin_integration_histogram,
													   param_obj.minNormFragCtDiff_integration_histogram,
													   "peaks",
													   '"{0}"'.format(param_obj.integration_summary_null_distribution_plots_dir))
		run_command(cmd)

	# make bed files for sub-additive, additive, and super-additive peaks
	upregulated_peaks_files = glob.glob(param_obj.upreg_peak_cats_bed_file_prefix[1:-1] + "*.bed")
	if run_all_steps or len(upregulated_peaks_files) == 0:
		cmd = "Rscript {0} {1} {2}".format(param_obj.path_to_make_bed_files_for_each_category, 
										   param_obj.upregulated_diffpeaks_output_file, 
										   param_obj.upreg_peak_cats_bed_file_prefix)
		run_command(cmd)

	# make new fragment counts files for each peakset
	upregulated_peaks_fragCount_r_objects = glob.glob(param_obj.upreg_peak_cats_bed_file_prefix[1:-1] + "*.rds")
	if run_all_steps or len(upregulated_peaks_fragCount_r_objects) == 0:
		for upreg_peak_file in upregulated_peaks_files:
			output_filename = upreg_peak_file + ".rds"
			cmd = "Rscript {0} {1} {2} {3} {4}".format(param_obj.path_to_createFragCountMatx,
											   '"{0}"'.format(upreg_peak_file),
											   param_obj.initial_peak_width, 
											   '"{0}"'.format(output_filename), 
											   "VariableWidth")
			run_command(cmd)

	# are the super-additive peaks more likely to be close to super-multiplicative genes?
	# make new joined peak tib
	if run_all_steps or not os.path.exists(param_obj.joined_gene_and_upreg_peak_table_file[1:-1]):
		cmd = "Rscript {0} {1} {2} {3}".format(param_obj.path_to_join_peak_to_gene_tib,
											   param_obj.path_to_upregulated_gene_table,
											   param_obj.upregulated_diffpeaks_output_file,
											   param_obj.joined_gene_and_upreg_peak_table_file)
		run_command(cmd)

		cmd = "Rscript {0} {1} {2} {3}".format(param_obj.path_to_join_peak_to_gene_tib,
											   param_obj.path_to_upregulated_gene_table,
											   param_obj.annotated_diffpeaks_output_file,
											   param_obj.joined_gene_and_all_peak_table_file)
		run_command(cmd)

	# use the new joined peak tib to make the "special peak types near genes" plots
	peaks_near_genes_analysis_plots = glob.glob(param_obj.peaks_near_genes_plots_dir + os.sep + "*.svg")
	if run_all_steps or len(peaks_near_genes_analysis_plots) == 0:
		cmd = "Rscript {0} {1} {2} {3} {4}".format(param_obj.path_to_make_peak_near_gene_analysis_plots,
											   param_obj.path_to_upregulated_gene_table,
											   param_obj.joined_gene_and_upreg_peak_table_file,
											   param_obj.joined_gene_and_all_peak_table_file,
											   '"{0}{1}"'.format(param_obj.peaks_near_genes_plots_dir, os.sep))
		run_command(cmd)

	# make the motif analysis plots
	motif_analysis_plots = glob.glob(param_obj.motif_analysis_plots_dir + os.sep + "*.svg")
	if run_all_steps or len(motif_analysis_plots) == 0:
		cmd = "Rscript {0} {1} {2} {3} {4}".format(param_obj.path_to_motif_analysis_plots,
												   param_obj.final_merged_differential_atac_frag_count_rds_file,
												   param_obj.upregulated_diffpeaks_output_file,
												   param_obj.most_variable_motifs_file,
												   '"{0}{1}"'.format(param_obj.motif_analysis_plots_dir, os.sep))
		run_command(cmd)

	# make the motif size and motif density by peak type (sub-additive, additive, super-additive) plots
	motif_analysis_plots = glob.glob(param_obj.motif_analysis_plots_dir + os.sep + "*.svg")
	if run_all_steps or len(motif_analysis_plots) != 5:
		cmd = "Rscript {0} {1} {2} {3}".format(param_obj.path_to_motif_size_and_density_analysis_plots,
											   param_obj.upregulated_diffpeaks_output_file,
											   param_obj.most_variable_motifs_file,
											   '"{0}{1}"'.format(param_obj.motif_analysis_plots_dir, os.sep))
		run_command(cmd)

	# make the expected vs. observed dual motif match rate plot
	motif_analysis_plots = glob.glob(param_obj.motif_analysis_plots_dir + os.sep + "*.svg")
	if run_all_steps or len(motif_analysis_plots) != 6:
		cmd = "Rscript {0} {1} {2}".format(param_obj.path_to_dual_motif_match_plot_script,
										   param_obj.upregulated_diffpeaks_output_file,
										   '"{0}{1}"'.format(param_obj.motif_analysis_plots_dir, os.sep))
		run_command(cmd)

	# make the supplemental motif analysis plot focusing on canonical TFs activated by each signal
	motif_analysis_plots = glob.glob(param_obj.motif_analysis_plots_dir + os.sep + "*.svg")
	if run_all_steps or len(motif_analysis_plots) != 10:
		cmd = "Rscript {0} {1} {2}".format(param_obj.path_to_supplemental_motif_analysis_plots,
										   param_obj.final_merged_differential_atac_frag_count_rds_file,
										   '"{0}{1}"'.format(param_obj.motif_analysis_plots_dir, os.sep))
		run_command(cmd)


if __name__ == '__main__':
	do_parameter_sweep = False
	run_all_steps = False

	if do_parameter_sweep:
		# basedir = str(here())
		# extractedDataDir = os.sep.join([basedir, "extractedData", "sensitivityAnalysis"])
		# plotsDir = os.sep.join([basedir, "plots", "sensitivityAnalysis"])

		# peak_merge_distances = [0, 250, 375, 500]           # also try [0, 150, 250, 500]
		# initial_peak_width           = 150                          #  these are centered at macs2 summits
		# final_diffpeak_algorithm_min_normalized_fragments_list = [10, 30]   # do 10, 20, 30, 40
		# final_diffpeak_algorithm_min_fold_change_list          = [1.1, 1.5]  # do 1,1, 1.5, 2, 4
		# param_object_list = []

		# for pmd in peak_merge_distances:
		# 	for ii in [0, 1]:
		# 		param_object_list.append(ParamSet(basedir, 
		# 										  extractedDataDir,
		# 										  plotsDir,
		# 										  peak_merge_distance = pmd, 
		# 										  initial_peak_width = 150,
		# 										  final_diffpeak_algorithm_min_normalized_fragments = final_diffpeak_algorithm_min_normalized_fragments_list[ii],
		# 										  final_diffpeak_algorithm_min_fold_change = final_diffpeak_algorithm_min_fold_change_list[ii])

		# for param_obj in param_object_list:
		# 	print("Running parameter set:\n{0}".format(str(param_obj)))
		# 	main(param_obj)
		pass
	else:
		# this parameter object stores the default parameters that we chose after parameter sweeps
		basedir = str(here())
		extractedDataDir = os.sep.join([basedir, "extractedData"])
		plotsDir = os.sep.join([basedir, "plots"])

		param_obj = ParamSet(basedir, 
							 extractedDataDir,
							 plotsDir,
							 peak_merge_distance = 250, 
							 initial_peak_width = 150,
							 initial_diffpeak_algorithm_min_normalized_fragments = 10,
							 initial_diffpeak_algorithm_min_fold_change = 1.1, 
							 final_diffpeak_algorithm_min_normalized_fragments = 30,
							 final_diffpeak_algorithm_min_fold_change = 1.5)
		main(param_obj, run_all_steps = run_all_steps)

	




