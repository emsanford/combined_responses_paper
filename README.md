This repository includes all analysis scripts used for the eLife publication, "Gene regulation gravitates towards either addition or multiplication when combining the effects of two signals", located at https://elifesciences.org/articles/59388

To generate the figures and intermediate data files associated with the publication, 

1. Download the raw and processed RNA-seq and ATAC-seq data files from https://www.dropbox.com/sh/fhx7huyhhtf8fux/AACKW5Bd7k34uy6Rrk3k0WZ4a?dl=0&lst=
2. Update the file paths in sampleMetadata_SI2-SI4.txt
3. Run the gene_analysis_pipeline.py script from the command line (requires python 2.7, R 3.5, and installation of associated packages imported in downstream scripts).
4. Run the peak_analysis_pipeline.py from the command line.
