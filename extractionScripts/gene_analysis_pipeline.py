import os
import glob
import re
import time
from pyprojroot import here


class ParamSet:
	def __init__(self, base_directory, extractedDataDir, plotsDir, refsDir, deseq_qval = 0.05, deseq_min_log2_fold_change = 2, use_default_param_string = False):
		self.base_directory      	 	= base_directory  # this should be the analysis root directory, which contains folders like "extractedData" and "plotScripts"
		self.sample_metadata_file		= '"{0}"'.format(base_directory + os.sep + "sampleMetadata_SI2-SI4.txt")
		self.extractedDataDir 			= extractedDataDir
		self.plotsDir       			= plotsDir
		self.refsDir                    = refsDir
		self.deseq_qval 				= deseq_qval
		self.deseq_min_log2_fold_change = deseq_min_log2_fold_change

		param_summary_string = "qval{0}_minlfc{1}".format(deseq_qval, deseq_min_log2_fold_change)
		if use_default_param_string:
			param_summary_string = "defaultParams"

		# directories that need to exist for saving results
		self.rna_seq_matrix_dir = os.sep.join([extractedDataDir, "rnaSeqMatrixFormatted"])
		self.subdirectories = [plotsDir, extractedDataDir, refsDir, self.rna_seq_matrix_dir]

		# paths to input and output files
		self.rnaseq_pipeline_counts_output_file      = '"{0}"'.format(os.sep.join([extractedDataDir, "si2-si4_RNA-seq-pipeline-output-counts.tsv"]))
		self.gene_counts_file_with_normalized_values = '"{0}"'.format(os.sep.join([extractedDataDir, "si2-si4_RNA-seq-pipeline-output-normalized.tsv"]))
		self.gene_counts_matrix_for_deseq            = '"{0}"'.format(os.sep.join([self.rna_seq_matrix_dir, "counts.RNA-seq-matrix.min-count-filtered.rds"]))
		self.deseq_output_table                      = '"{0}"'.format(os.sep.join([extractedDataDir, "DeSeqOutputAllConds.tsv"]))
		self.annotated_deseq_output_table            = '"{0}"'.format(os.sep.join([extractedDataDir, "DeSeqOutputAllConds.annotated.tsv"]))
		self.upregulated_genes_table                 = '"{0}"'.format(os.sep.join([extractedDataDir, "DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv"]))


		# paths to reference files
		self.hg38_gtf_file                      = '"{0}"'.format(os.sep.join([refsDir, "hg38.gtf"]))  # Homo_sapiens.GRCh38.90.gtf.gz
		self.ensg_to_hgnc_symbol_mapping        = '"{0}"'.format(os.sep.join([refsDir, "EnsgHgncSymbolMapping.tsv"]))  # we use the gene name from the gtf file if there's no HGNC symbol for it
		self.ensg_to_hg38_canonical_tss_mapping = '"{0}"'.format(os.sep.join([refsDir, "EnsgToHg38CanonicalTssMapping.tsv"])) 


		# paths to scripts, we need to add quotes around them to pass as command line arguments because dropbox added a space to its own folder and we're using os.system to run commands
		self.path_to_normalizePipelineCountsOutputAndAddGeneSymbol = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "normalizePipelineCountsOutputAndAddGeneSymbol.R"]))
		self.makeGeneExpressionMatrixWithMinCounts                 = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "makeGeneExpressionMatrixWithMinCounts.R"]))
		self.path_to_runDESeqOnConditionSet                        = '"{0}"'.format(os.sep.join([base_directory, "extractionScripts", "runDESeqOnConditionSet.R"]))
		self.path_to_zzzzzzz          = '"{0}"'.format(os.sep.join([base_directory, "plotScripts", "makeAdjacentSuperadditivePeakAssocModeOfIntegrationPlots.R"]))



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

	if run_all_steps or not os.path.exists(param_obj.gene_counts_file_with_normalized_values[1:-1]):
		cmd = 'Rscript {0} {1} {2} {3} {4} {5}'.format(param_obj.path_to_normalizePipelineCountsOutputAndAddGeneSymbol,
													   param_obj.rnaseq_pipeline_counts_output_file,
													   param_obj.hg38_gtf_file,
													   param_obj.ensg_to_hgnc_symbol_mapping,
													   param_obj.ensg_to_hg38_canonical_tss_mapping,
													   param_obj.gene_counts_file_with_normalized_values)
		run_command(cmd)

	matrix_output_files = glob.glob(param_obj.rna_seq_matrix_dir + '/*.rds')
	if run_all_steps or not len(matrix_output_files) == 3:
		cmd = 'Rscript {0} {1} {2}'.format(param_obj.makeGeneExpressionMatrixWithMinCounts,
										   param_obj.gene_counts_file_with_normalized_values,
										   '"{0}"'.format(param_obj.rna_seq_matrix_dir))
		run_command(cmd)
	
	if run_all_steps or not os.path.exists(param_obj.deseq_output_table[1:-1]):
		cmd = 'Rscript {0} {1} {2}'.format(param_obj.path_to_runDESeqOnConditionSet,
										   param_obj.gene_counts_matrix_for_deseq,
										   param_obj.sample_metadata_file,
										   param_obj.deseq_output_table)
		run_command(cmd)


	

if __name__ == '__main__':
	do_parameter_sweep = False
	run_all_steps      = False

	if do_parameter_sweep:
		pass
	else:
		# this parameter object stores the default parameters that we chose after parameter sweeps
		basedir = str(here())
		extractedDataDir = os.sep.join([basedir, "extractedData"])
		plotsDir = os.sep.join([basedir, "plots"])
		refsDir  = os.sep.join([basedir, "refs"])

		param_obj = ParamSet(basedir, 
							 extractedDataDir,
							 plotsDir,
							 refsDir,
							 deseq_qval = 0.05, 
							 deseq_min_log2_fold_change = 2, 
							 use_default_param_string = True)
		main(param_obj, run_all_steps = run_all_steps)

	




