Cistrome-wide association studies (CWAS)

This repository contains scripts for the cistrome-wide association studies described in https://www.biorxiv.org/content/10.1101/2021.05.10.443466v1. The concept of CWAS is to identify chromatin (cistrome) features that are genetically associated with a trait of interest.

The basic steps are:
1. Use reference epigenome data from many individuals to model chromatin activity as a function of nearby SNP genotypes
2. Associate predicted chromatin activity with a trait of interest using gwas summary statistics

Please see https://www.biorxiv.org/content/10.1101/2021.05.10.443466v1 for more information

This analysis is incorporated into a snakemake workflow in cwas.snakefile. The analysis can be performed in a conda environment whose recipe is provided in env/cwas_env.yml

To run the analysis:
1. Create and activate the conda environment specified by env/cwas_env.yml
2. Update the config.yaml file with the locations of reference data files. These include, for each sample, (a) a bam file with epigenomic data (eg, H3K27ac ChIP-seq), (b) a bed file with the coordinates of peaks from the epigenomic data and (c) a vcf file with PHASED genotypes. For this project, phasing was performed with Eagle2 using the Sanger Imputation Service (https://imputation.sanger.ac.uk/). Example config files are provided in run_files/
3. Download the needed data / repositories. This includes: 
(a) The WASP pipeline for mitigating mapping bias for reads covering heterozygous SNPs (https://github.com/bmvdgeijn/WASP). The location of this directory must be specified in the config file. 
(b) Download GATK3 package (package-archive_gatk_GenomeAnalysisTK-3.8-0-ge9d806836.tar from https://software.broadinstitute.org/gatk/download/archive. Then expand the .tar file and activate it using "gatk3-register {path to gatk3 jar file}".
(c) LDREF data from 1000 genomes (https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2). This goes in the LDREF folder. 
(d) GWAS summary statistics files, with file names *.sumstats.gz
These files go in gwas_data/ and have the following format:

SNP	A1	A2	N	Z
rs7899632	A	G	140254	-2.98765
rs3750595	A	C	140254	2.98765
rs10786405	T	C	140254	2.90123

4. Activate the cwas conda environment (ie, with "conda activate cwas")

5. Run the snakemake workflow. To submit to a compute cluster with slurm workload manager, use sbatch submit.sh. Output will populate in the analysis folder

