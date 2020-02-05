import glob
import os
import re
import numpy as np
from csv import DictReader

sample_metadata_file = "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/sampleMetadata_SI2-SI4.txt"
chrom_sizes_file = "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/refs/hg38.chrom.sizes"
path_to_bigWigToBedGraph = "/Users/emsanford/anaconda2/bin/bigWigToBedGraph"  # installed with: "conda install -c bioconda/label/cf201901 ucsc-bigwigtobedgraph"
path_to_bedGraphToBigWig = "/Users/emsanford/anaconda2/bin/bedGraphToBigWig"  # installed with: "conda install -c bioconda/label/cf201901 ucsc-bedgraphtobigwig"
output_dir_prefix = "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/extractedData/IGV_tracks_avgFragmentCountsPerMillion/"

sampleCond_bedGraph_file_dict = {}
nConditions = 12

# this step converts the macs2 output bigwig file back into a bedgraph file
dr = DictReader(open(sample_metadata_file), delimiter='\t')
for dl in dr:
	sample_number = int(dl["sampleID"][0:2])
	if 1 <= sample_number <= 36:
		# print(dl["ATAC_analysisDir"])
		glob_pattern = dl["ATAC_analysisDir"] + "macs2*/*.bigWig"
		# print(glob_pattern)
		glob_results = glob.glob(glob_pattern)
		assert(len(glob_results) == 1)
		this_bigWig_file = glob_results[0]
		outputBedGraphName = this_bigWig_file + '.ToBedGraph.bdg'

		# save the output bedgraph file in a dictionary, to be used later to average together the results from each replicate
		if not dl["condition"] in sampleCond_bedGraph_file_dict:
			sampleCond_bedGraph_file_dict[dl["condition"]] = [outputBedGraphName]
		else:
			sampleCond_bedGraph_file_dict[dl["condition"]].append(outputBedGraphName)

		cmd = "{0} {1} {2}".format(path_to_bigWigToBedGraph, this_bigWig_file, outputBedGraphName)
		# print(cmd)
		# os.system(cmd)  # comment out this line to skip the actual conversation of bigWig files to bedgraph files

# this step takes the bedgraph files for each condition and merges them into a single file
for condition in sampleCond_bedGraph_file_dict.keys():
	# this step is done with bedtools v2.27.1
	nreps = 3
	assert(len(sampleCond_bedGraph_file_dict[condition]) == nreps)
	output_file_name = output_dir_prefix + condition + ".unionBedGraphFile.bdg"
	cmd = 'bedtools unionbedg -i "{0}" > "{1}"'.format('" "'.join(sampleCond_bedGraph_file_dict[condition]), output_file_name)
	# print(cmd)
	# os.system(cmd)

# this step iterates over each of the merged bedgraph files, computes an average fragment counts per million, then produces a new single bedgraph file
unionbedg_file_list = glob.glob(output_dir_prefix + "*.unionBedGraphFile.bdg")
assert(len(unionbedg_file_list) == nConditions)
for filepath in unionbedg_file_list:
	regexmatchobj = re.match(".*/(.*).unionBedGraphFile.bdg", filepath)
	condname = regexmatchobj.group(1)
	output_filepath = output_dir_prefix + condname + ".averagedBedgraphFile.bdg"
	print("calculating average values for {0}".format(condname))
	fout = open(output_filepath, "w")
	with open(filepath) as f:
		for line in f:
			linefields = line.strip().split('\t')
			fragments_per_million_each_sample = [float(x) for x in linefields[-3:]]
			mean_fragment_count = np.mean(fragments_per_million_each_sample)
			output_list_this_line = linefields[:3] + ["{0:0.08f}".format(mean_fragment_count)]
			output_string = '\t'.join(output_list_this_line) + '\n'
			fout.write(output_string)
	fout.close()
	

# this step takes the output of the last step and converts it back into a bigwig file, (for fast IGV viewing)
averagedBedgraphFileList = glob.glob(output_dir_prefix + "*.averagedBedgraphFile.bdg")
assert(len(averagedBedgraphFileList) == nConditions)
for filepath in averagedBedgraphFileList:
	regexmatchobj = re.match(".*/(.*).averagedBedgraphFile.bdg", filepath)
	condname = regexmatchobj.group(1)
	output_filepath = output_dir_prefix + condname + ".atacFragmentsPerMillion.bigWig"
	cmd = "{0} {1} {2} {3}".format(path_to_bedGraphToBigWig, '"{0}"'.format(filepath), '"{0}"'.format(chrom_sizes_file), '"{0}"'.format(output_filepath))
	print(cmd)
	os.system(cmd)














