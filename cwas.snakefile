configfile: "config.yaml"

import pandas as pd
import yaml
import os

def get_chr():
	return ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
              "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
              "chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX", "chrY"]

def get_bams(wildcards):
	return config["bams"][wildcards.sample]

def get_vcfs(wildcards):
	return config["vcfs"][wildcards.sample]

def get_beds(wildcards):
        return config["beds"][wildcards.sample]

def local_param_targets(wildcards):
	ls = []
	for sample in config["bams"]:
		ls.append("analysis/stratas/%s/%s.local.params" % (sample,sample))
	return ls

def global_param_targets(wildcards):
	ls = []
	for sample in config["bams"]:
		ls.append("analysis/stratas/%s/%s.global.params" % (sample,sample))
	return ls

def counts_targets(wildcards):
	ls = []
	for sample in config["bams"]:
		ls.append("analysis/stratas/%s/%s.counts.mat" % (sample,sample))
	return ls

def quan_targets(wildcards):
	ls = []
	for sample in config["bams"]:
		ls.append("analysis/qtl/%s/%s.exon.rpkm.bed" % (sample,sample))
	return ls

def vcf_list(wildcards):
	ls = []
	for sample in config["vcfs"]:
		ls.append("vcf_links/%s.vcf.gz" % (sample))
	return ls

# create targets for "all" rule
def get_targets(wildcards):
	ls = []
	ls.append("analysis/qtl/cis/permute-significant.tsv.gz") 
	ls.append("analysis/qtl/cis/permute-check_beta_approx-neglog10.png")
	ls.append("analysis/qtl/summary/qtl_ai_overlap.txt")
	ls.append("analysis/qtl/summary/beta_vs_mu.pdf")
	ls.append("analysis/qtl/summary/peaks/combined.sig.bed")
	ls.append("analysis/qtl/summary/snps/combined.txt") 
	ls.append("analysis/enrich/results/enrichment.pdf")
	ls.append("analysis/motifs/ai_vs_motifs/motif1.ai.pdf")
	ls.append("analysis/motifs/ai_vs_motifs/motif2.ai.pdf")
	ls.append("analysis/motifs/ai_vs_motifs/motif3.ai.pdf")
	ls.append("analysis/motifs/ai_vs_motifs/motif4.ai.pdf")
	ls.append("analysis/motifs/ai_vs_motifs/motif5.ai.pdf")
	ls.append("analysis/motifs/ai_vs_motifs/motif6.ai.pdf")
	ls.append("analysis/motifs/ai_vs_motifs/motifX.ai.pdf")
	ls.append("analysis/predict/results/predict.all.txt")
#	ls.append("analysis/hsq/hsq.summary.txt")
	ls.append("analysis/fusion/ProstateCancer_Meta_Schumacher2018.nodup/plots/manhattan.pdf")
#	ls.append("analysis/fusion/ProstateCancer_Meta_Schumacher2018.nodup/plots/gwas_twas_cwas_plot.pdf")
	ls.append("analysis/fusion/ProstateCancer_Meta_Schumacher2018.nodup/merged/summary.txt")
	ls.append("analysis/fusion/ProstateCancer_Meta_Schumacher2018.nodup/postprocess/best/conditioned.summary.txt")
	ls.append("analysis/fusion/ADT.gwas.chi2/merged/sig/cwas.sig.gwas.restricted.txt")
	ls.append("analysis/fusion/bph/merged/summary.txt")
	ls.append("analysis/fusion/bph/postprocess/best/conditioned.summary.txt")
	ls.append("analysis/fusion/compare_gwas/compare.gwas.pdf")
	return ls

rule all:
	input:
		get_targets

# --- rules for data preparation and read remapping at heterozygous SNPs

rule bams_links_and_index:
	"""create links to bam files and index bam files"""
	input:
		bam = get_bams
	output:
		bam_ln = "bam_links/{sample}.bam",
		ix = "bam_links/{sample}.bam.bai"
	shell:
		"ln -s $PWD/{input.bam} {output.bam_ln} && "
		"samtools index {output.bam_ln}"

rule vcf_link_and_index:
	"""subset vcf files to only SNPs and create vcf index"""
	input:
		get_vcfs
	output:
		vcf = "vcf_links/{sample}.vcf.gz",
		csi = "vcf_links/{sample}.vcf.gz.csi",
		tbi = "vcf_links/{sample}.vcf.gz.tbi",
	shell:
		""" bcftools view -i 'TYPE="snp"' -Oz {input} > {output.vcf}
		bcftools index {output.vcf}
		tabix {output.vcf} """

rule split_vcf_by_chr:
	"""split vcf file with SNPs by chromosome"""
	input:
		"vcf_links/{sample}.vcf.gz"
	output:
		["analysis/wasp/variants/{sample}/" + "%s.vcf.gz" % s for s in get_chr()]
	shell:
		"""for c in {{1..22}} X Y MT; do bcftools view -Oz -r $c {input} > analysis/wasp/variants/{wildcards.sample}/chr${{c}}.vcf.gz; done"""

rule prep_snps:
	"""create list of snps in the format expected by wasp"""
	input:
		["analysis/wasp/variants/{sample}/" + "%s.vcf.gz" % s for s in get_chr()]
	output:
		snps = ["analysis/wasp/variants/{sample}/" + "%s.snps.txt.gz" % s for s in get_chr()]	
	shell:
		"sh {config[wasp_dir]}/mapping/extract_vcf_snps.sh analysis/wasp/variants/{wildcards.sample} analysis/wasp/variants/{wildcards.sample}"

rule find_intersecting_snps:
	"""find reads that intersect SNPs for remapping by WASP"""
	input:
		bam = "bam_links/{sample}.bam",
		snps = ["analysis/wasp/variants/{sample}/" + "%s.snps.txt.gz" % s for s in get_chr()]
	output:
		fastq1 = "analysis/wasp/find_intersecting_snps/{sample}/{sample}.remap.fq1.gz",
		keep_bam = temp("analysis/wasp/find_intersecting_snps/{sample}/{sample}.keep.bam"),
		remap_bam = temp("analysis/wasp/find_intersecting_snps/{sample}/{sample}.to.remap.bam"),
		log = "analysis/wasp/find_intersecting_snps/{sample}/{sample}.log"
	params:
		seq_type = config["seq_type"]
	shell:
		""" if [[ "{params.seq_type}" == "single" ]]
		then
			python {config[wasp_dir]}/mapping/find_intersecting_snps.py --is_sorted --output_dir analysis/wasp/find_intersecting_snps/{wildcards.sample} --snp_dir analysis/wasp/variants/{wildcards.sample} {input.bam} > {output.log}
			mv analysis/wasp/find_intersecting_snps/{wildcards.sample}/{wildcards.sample}.remap.fq.gz {output.fastq1}
		elif [[ "{params.seq_type}" == "paired" ]]
		then
			python {config[wasp_dir]}/mapping/find_intersecting_snps.py  --is_paired_end --is_sorted --output_dir analysis/wasp/find_intersecting_snps/{wildcards.sample} --snp_dir analysis/wasp/variants/{wildcards.sample} {input.bam} > {output.log}
		fi """

# NOTE: will need to update this if using other aligners
rule remap_at_variants:
	"""remap reads with alternate allele"""
	input:
		fastq1 = "analysis/wasp/find_intersecting_snps/{sample}/{sample}.remap.fq1.gz"
	output:
		bam = temp("analysis/wasp/remap/{sample}/{sample}.bam"),
		log = "analysis/wasp/find_intersecting_snps/logs/{sample}.remap.log"
	params:
		seq_type = config["seq_type"]
	shell:
		""" if [[ "{params.seq_type}" == "single" ]]
		then
			bwa mem {config[bwa_index]} {input.fastq1} | samtools view -b > {output.bam} 2> {output.log}
		elif [[ "{params.seq_type}" == "paired" ]]	
		then
			bwa mem {config[bwa_index]} {input.fastq1} analysis/wasp/find_intersecting_snps/{wildcards.sample}/{wildcards.sample}.remap.fq2.gz  | samtools view -b > {output.bam} 2> {output.log}
		fi
		"""

rule sort_and_index_bam:
	"""sort and index bam generated by remapping step"""
	input:
		"analysis/wasp/remap/{sample}/{sample}.bam"
	output:
		temp("analysis/wasp/remap/{sample}/{sample}.sorted.bam"),
		temp("analysis/wasp/remap/{sample}/{sample}.sorted.bam.bai")
	shell:
		"samtools sort -o {output[0]} {input}; "
		"samtools index {output[0]}"

rule filter_remapped_reads:
	"""filter reads from second mapping step"""
	input:
		to_remap_bam = "analysis/wasp/find_intersecting_snps/{sample}/{sample}.to.remap.bam",
		remap_bam = "analysis/wasp/remap/{sample}/{sample}.sorted.bam",
		bai = "analysis/wasp/remap/{sample}/{sample}.sorted.bam.bai"
	output:
		keep_bam = temp("analysis/wasp/filter_remapped_reads/{sample}/{sample}.keep.bam"),
		log = "analysis/wasp/filter_remapped_reads/{sample}/{sample}.log"
	shell:
		"python {config[wasp_dir]}/mapping/filter_remapped_reads.py "
		"  {input.to_remap_bam} {input.remap_bam} {output.keep_bam} > {output.log}"

    
rule merge_bams:
	"""merge 'keep' BAM files from mapping steps 1 and 2, then sort and index"""
	input:
		keep1 = "analysis/wasp/find_intersecting_snps/{sample}/{sample}.keep.bam",
		keep2 = "analysis/wasp/filter_remapped_reads/{sample}/{sample}.keep.bam"
	output:
		merge = temp("analysis/wasp/merge/{sample}/{sample}.keep.merge.bam"), 
		sort = temp("analysis/wasp/merge/{sample}/{sample}.keep.merge.sort.bam") 
	shell:
		"samtools merge {output.merge} {input.keep1} {input.keep2}; "
		"samtools sort -o {output.sort} {output.merge}; "
		"samtools index {output.sort}"

rule rmdup:
	"""remove duplicate reads"""
	input:
		"analysis/wasp/merge/{sample}/{sample}.keep.merge.sort.bam"
	output:
		rmdup = temp("analysis/wasp/rmdup/{sample}/{sample}.keep.merge.rmdup.bam"),
		sort = temp("analysis/wasp/rmdup/{sample}/{sample}.keep.merge.rmdup.sort.bam")
	params:
		seq_type = config["seq_type"]
	shell:
		""" if [[ "{params.seq_type}" == "single" ]]
		then
			python {config[wasp_dir]}/mapping/rmdup.py {input} {output.rmdup}
		elif [[ "{params.seq_type}" == "paired" ]]
		then
			python {config[wasp_dir]}/mapping/rmdup_pe.py {input} {output.rmdup}
		fi
		samtools sort -o {output.sort} {output.rmdup}
		samtools index {output.sort} """

rule prep_bam_for_asereader:
	"""prepare bam for analysis by ASEReadCounter"""
	input:
		"analysis/wasp/rmdup/{sample}/{sample}.keep.merge.rmdup.sort.bam"
	output:
		bam = temp("analysis/asereader/{sample}/{sample}.remapped.bam"),
		bai = "analysis/asereader/{sample}/{sample}.remapped.bam.bai"
	params:
		genome = config["dict"]
	shell:
		"sh scripts/prep_bam_for_asereader.sh {input} {output.bam} {wildcards.sample} analysis/asereader/{wildcards.sample} {params.genome}"

# --- rules for running stratAS to identify allelic imbalance

rule allele_counts:
	"""get allele-specific read counts for SNPs"""
	input:
		bam = "analysis/asereader/{sample}/{sample}.remapped.bam",
		bai = "analysis/asereader/{sample}/{sample}.remapped.bam.bai",
		vcf = "vcf_links/{sample}.vcf.gz"
	output:
		"analysis/asereader/{sample}/{sample}.allele.counts.tsv"
	params:
		genome = config["genome"]
	shell:
		"sh scripts/ase_readcounter.sh {input.bam} {input.vcf} {params.genome} {output}"

rule consensus_peaks:
	"""identify consensus peaks from the provided bed files"""
	input:
		config["beds"].values()
	output:
		"analysis/consensus_peaks/consensus.peaks.stratas.tsv"
	params:
		n = config["min_samples_for_peak"],
		chr_sizes = config["chr_sizes"]
	shell:
		""" cat {input} | sort -k1,1 -k2,2n > analysis/consensus_peaks/cat.bed
		sh scripts/consensus.peaks.sh analysis/consensus_peaks/cat.bed {params.chr_sizes} {params.n}
		"""

rule prep_peaks:
	"""prepare files with peak coordinates for use by copywriteR"""
	input:
		get_beds	
	output:
		os.getcwd() + "/analysis/cna/{sample}/{sample}.peaks.bed" # need full path name for copywriter bed file
	shell:
		"grep -v chrM {input} | sed 's/chr//' > {output}"

rule copywriter:
	"""use copywriteR to estimate copy number from off-target reads"""
	input:
		bam = "analysis/asereader/{sample}/{sample}.remapped.bam",
		bed = os.getcwd() + "/analysis/cna/{sample}/{sample}.peaks.bed" # need full path name for copywriter bed file
	output:
		"analysis/cna/{sample}/{sample}.igv"
	shell:
		# must remove the directory if it already exists to prevent copywriter from throwing an error
		"""if [ -d "analysis/cna/{wildcards.sample}/CNAprofiles" ]; then rm -Rf analysis/cna/{wildcards.sample}/CNAprofiles; fi
		Rscript scripts/copywriter.R analysis/asereader/{wildcards.sample} {wildcards.sample} analysis/cna/{wildcards.sample} {input.bed}"""

rule stratas_cna_file:
	"""prepare copy number alteration (CNA) file for use with stratAS"""
	input:
		"analysis/cna/{sample}/{sample}.igv"
	output:
		"analysis/cna/{sample}/{sample}.stratas.cna"
	shell:
		"printf 'CHR\tP0\tP1\tCNV\n' > analysis/cna/{wildcards.sample}/tmp.header.txt && tail -n +27 {input} | cut -f1-3,5 >> analysis/cna/{wildcards.sample}/tmp.header.txt && mv analysis/cna/{wildcards.sample}/tmp.header.txt {output}"
# tail -n +27 is to remove the headers/comments plus first 24 rows of the .igv file, which are extraneous

rule stratas_input:
	"""generate input files needed by stratAS"""
	input:
		vcf = "vcf_links/{sample}.vcf.gz",
		counts = "analysis/asereader/{sample}/{sample}.allele.counts.tsv",
		cna = "analysis/cna/{sample}/{sample}.stratas.cna"
	output:
		"analysis/stratas/{sample}/{sample}.local.params",
		"analysis/stratas/{sample}/{sample}.global.params",
		"analysis/stratas/{sample}/{sample}.counts.mat"
	shell:
		"sh scripts/prep_stratas_input.sh {wildcards.sample} {input.vcf} {input.counts} {input.cna}"

rule merge_params:
	"""merge stratAS parameter files for all samples"""
	input:
		local_params = local_param_targets,
		global_params = global_param_targets
	output:
		local_params = "analysis/stratasrun/all.local.params",
		global_params = "analysis/stratasrun/all.global.params",
		plot1 = "analysis/stratasrun/local.params.pdf",
		plot2 = "analysis/stratasrun/global.params.pdf"
	shell:
		""" printf 'ID\tPHI\tMU\tN\n' > {output.global_params}
		cat {input.global_params} | grep -v $'^ID\tPHI' >> {output.global_params}
		printf 'ID\tCHR\tP0\tP1\tCNV\tPHI\tMU\tN\n' > {output.local_params}
		cat {input.local_params} | grep -v $'^ID\tCHR' >> {output.local_params} 
		Rscript scripts/plot_params.R {output.local_params} {output.global_params} {output.plot1} {output.plot2}"""

rule merge_counts:
	"""merge stratAS-formatted read count files for all samples"""
	input:
		counts_targets
	output:
		"analysis/stratasrun/all.counts.mat"
	shell:
		"""cut -f1-5 $(echo {input} | sed 's/ .*//') -d ' ' > analysis/stratasrun/out.tmp 
		for f in {input}
		do 
			echo merging $f
			cut -d ' ' -f6-9 $f | paste -d ' ' analysis/stratasrun/out.tmp - > analysis/stratasrun/out2.tmp 
			mv analysis/stratasrun/out2.tmp analysis/stratasrun/out.tmp
		done
		grep -v CHR analysis/stratasrun/out.tmp > {output} 
		rm analysis/stratasrun/out.tmp """

rule split_counts:
	"""split up counts files so they can be run in parallel"""
	input:
		"analysis/stratasrun/all.counts.mat"
	output:
		["analysis/stratasrun/split/counts.%02d" % s for s in range(0,50)]
	shell:
		"split -d -n l/50 {input} analysis/stratasrun/split/counts."

rule run_stratas:
	"""identify allelic imbalance with stratAS"""
	input:
		peaks = "analysis/consensus_peaks/consensus.peaks.stratas.tsv",
		local_params = "analysis/stratasrun/all.local.params",
		global_params = "analysis/stratasrun/all.global.params",
		split = "analysis/stratasrun/split/counts.{split}"
	output:
		"analysis/stratasrun/results/results.{split}"
	params:
		samples = config["samples"]
	shell:
		"Rscript stratAS/stratas.R --input {input.split} --global_param {input.global_params} --local_param {input.local_params} --peaks {input.peaks} --samples {params.samples} --max_rho 0.2 --window -1 --indiv TRUE --min_cov 0 --seed 123 > analysis/stratasrun/results/results.{wildcards.split}"

rule merge_results:
	"""merge stratAS output"""
	input:
		["analysis/stratasrun/results/results.%02d" % s for s in range(0,50)]
	output:
		"analysis/stratasrun/results/results.all.txt"
	shell:
		""" head -n 1 $(echo {input} | sed 's/ .*//') > {output}
		cat {input} | grep -v CHR >> {output} """

rule summarize_results:
	"""summarize allelic imbalance results"""
	input:
		"analysis/stratasrun/results/results.all.txt",
	output:
		"analysis/summary/sig.BOTH.txt",
		"analysis/summary/bed/sig.BOTH.peak.bed",
		plot = "analysis/summary/allele_counts_plot.pdf"
	params:
		samples = config["samples"]
	shell:
		""" Rscript scripts/summarize_results.R {input} {params.samples} {output.plot}
		cp {input.peaks} analysis/summary/bed """

# --- rules for identifying chromatin QTLs

rule bed_to_gtf:
	"""convert bed file with peaks to gtf format"""
	input:
		"analysis/consensus_peaks/consensus.peaks.stratas.tsv"
	output:
		"analysis/qtl/peaks/peaks.gtf"
	shell:
		""" sh scripts/bed_to_gtf.sh {input} {output} """
	
rule qtl_quan:
	"""get normalized read counts at peaks"""
	input:
		bam = "analysis/asereader/{sample}/{sample}.remapped.bam",
		bai = "analysis/asereader/{sample}/{sample}.remapped.bam.bai",
		gtf = "analysis/qtl/peaks/peaks.gtf"
	output:
		"analysis/qtl/{sample}/{sample}.exon.rpkm.bed"
	params:
		qtl = config["qtltools"]
	shell:
		""" {params.qtl} quan --bam {input.bam}  --gtf {input.gtf} --rpkm --filter-mismatch 5 --filter-mismatch-total 5 --filter-mapping-quality 30 --out-prefix analysis/qtl/{wildcards.sample}/{wildcards.sample} --no-hash
		sed -i 's/analysis.*/{wildcards.sample}/' {output} """

rule merge_quan:
	"""merge normalized read counts at peaks across all samples"""
	input:
		quan_targets
	output:
		"analysis/qtl/quan/rpkm.merged.bed.gz"
	shell:
		""" cut -f1-6 $(echo {input} | sed 's/ .*//') > analysis/qtl/quan/out.tmp 
		for f in {input}
		do
			echo merging $f
			cut -f7 $f | paste analysis/qtl/quan/out.tmp - > analysis/qtl/quan/out2.tmp
			mv analysis/qtl/quan/out2.tmp analysis/qtl/quan/out.tmp
		done
		# this requires combined rpkm > 10 across samples to prevent error from gsl gamm.c:
		awk 'BEGIN{{FS="\t"}} NR==1{{print $0}} NR>1 {{for(i=7; i<=NF; i++) t+=$i; if(t>10) {{print $0}}; t=0}}' \
		  analysis/qtl/quan/out.tmp > analysis/qtl/quan/rpkm.merged.bed # rm empty rows
		bgzip analysis/qtl/quan/rpkm.merged.bed
		tabix analysis/qtl/quan/rpkm.merged.bed.gz """

rule qtl_pca:
	"""pca of normalized read counts at peaks"""
	input:
		"analysis/qtl/quan/rpkm.merged.bed.gz"
	output:
		"analysis/qtl/quan/covariate.pca"
	params:
		qtl = config["qtltools"]
	shell:
		""" {params.qtl} pca --bed analysis/qtl/quan/rpkm.merged.bed.gz --scale --center --out analysis/qtl/quan/covariate """

# may want to update this - in theory it could try to run before indexing of vcf files finishes
rule qtl_prep_vcf:
	"""prepare vcf files for qtl analysis"""
	input:
		vcf_list
	output:
		"analysis/qtl/cis/all.samples.vcf.gz"
	shell:
		""" bcftools merge {input} -Oz --force-samples -o {output} 
		echo {input} | sed 's/vcf_links\///g' | \
		sed 's/.vcf.gz /\\n/g' | sed 's/.vcf.gz//' > list.tmp
		bcftools reheader {output} -s list.tmp > analysis/qtl/cis/vcf.tmp.gz
		mv analysis/qtl/cis/vcf.tmp.gz {output} 
		tabix {output}
		rm list.tmp  """

# adapted from snakemake-qtl (https://github.com/letaylor/snakemake-qtl)
rule qtl_cis:
	"""identify chromatin QTLs"""
	input:
		vcf = "analysis/qtl/cis/all.samples.vcf.gz",
		bed = "analysis/qtl/quan/rpkm.merged.bed.gz",
		cov = "analysis/predict/covar.mat" #has just first 6 PCs
	output:
		temp('analysis/qtl/cis/nominals.{j_cur}_{j_total}.txt')
	params:
		qtl = config["qtltools"]
	shell:
		""" {params.qtl} cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov}\
		--nominal 1 --normal --window 25000 --chunk {wildcards.j_cur} {wildcards.j_total} --out {output} """

# adapted from snakemake-qtl (https://github.com/letaylor/snakemake-qtl)
rule qtl_cis_concat:
	"""merge cQTL results"""
	input:
		expand(
			'analysis/qtl/cis/nominals.{j_cur}_{j_total}.txt', 
			j_cur = range(1, 51), 
			j_total = 50
		)
	output:
		'analysis/qtl/cis/nominals.tsv.gz'
	shell:
		# add the header to the output file
		'echo "pheno_id '
			'pheno_chr pheno_start pheno_end pheno_strand '
			'n_proximal_var pheno_var_dist '
			'var_id var_chr var_start var_end '
			'p_nominal beta top_proximal_var" | '
			'sed s/" "/"\t"/g | gzip -c > {output}; '
            
		# cat swarm files to make the output file
		'cat {input} | sed s/" "/"\t"/g | gzip -c >> {output}; '

# adapted from snakemake-qtl (https://github.com/letaylor/snakemake-qtl)
rule qtl_cis_permute:
	"""identify chromatin QTLs and calculate permutaion-based p-values"""
	input:
		vcf = "analysis/qtl/cis/all.samples.vcf.gz",
		bed = "analysis/qtl/quan/rpkm.merged.bed.gz",
		cov = "analysis/predict/covar.mat" #has just first 6 PCs
	output:
		temp('analysis/qtl/cis/permute.{j_cur}_{j_total}.txt')
	params:
		qtl = config["qtltools"]
	shell:
		# NOTE: increased window size to 1e6 because phenotypes with only 0 or 1 variants cause an error with gsl gamm.c
		""" {params.qtl} cis --permute 1000 --vcf {input.vcf} --bed {input.bed} --cov {input.cov} \
		--normal --window 1000000 --seed 123 --chunk {wildcards.j_cur} {wildcards.j_total} --out {output} """

# adapted from snakemake-qtl (https://github.com/letaylor/snakemake-qtl)
rule qtl_cis_permute_concat:
	"""merge chromatin QTL permutation results"""
	input:
		expand(
			'analysis/qtl/cis/permute.{j_cur}_{j_total}.txt', 
			j_cur = range(1, 51), 
			j_total = 50
		)
	output:
		'analysis/qtl/cis/permute.tsv.gz'
	shell:
		# add the header to the output file        
		'echo "pheno_id '
		'pheno_chr pheno_start pheno_end pheno_strand '
		'n_proximal_var pheno_var_dist '
		'var_id var_chr var_start var_end '
		'p_degree_freedom '
		'dummy '
		'beta_dist_par1 beta_dist_par2 '
		'p_nominal beta '
		'p_empirical p_adjust_beta_dist" | '
		'sed s/" "/"\t"/g | gzip -c > {output}; '
            
		# cat swarm files to make the output file
		'cat {input} | sed s/" "/"\t"/g | gzip -c >> {output}; '

# adapted from snakemake-qtl (https://github.com/letaylor/snakemake-qtl)
rule qtl_cis_permute_fdr:
	"""FDR adjustment of permutation-based p-values for chromatin QTLs"""
	input:
		'analysis/qtl/cis/permute.tsv.gz'
	output:
		'analysis/qtl/cis/permute-significant.tsv.gz',
		'analysis/qtl/cis/permute-thresholds.txt.gz'
	shell:
		'Rscript scripts/qtltools-runFDR_cis.R {input} 0.05 analysis/qtl/cis/`basename {input} .tsv.gz`'

# adapted from snakemake-qtl (https://github.com/letaylor/snakemake-qtl)
rule qtl_check_beta_approx:
	"""check the beta approximation used by qtltools"""
	input:
		'analysis/qtl/cis/permute.tsv.gz'
	output:
		'analysis/qtl/cis/permute-check_beta_approx.png',
		'analysis/qtl/cis/permute-check_beta_approx-neglog10.png'
	shell:
		'Rscript scripts/qtltools-check_beta_approx.R {input} '
		'analysis/qtl/cis/`basename {input} .tsv.gz`-check_beta_approx'

rule qtl_snps_peaks_distances:
	"""plot histogram of distance to top SNP for significant qtl peaks"""
	input:
		"analysis/qtl/cis/permute-significant.tsv.gz"
	output:
		"analysis/qtl/summary/snp.dist.pdf"
	shell:
		""" zcat {input} | sed '1d' | awk '{{print $1 "\t" $2 "\t" $7 "\t" $20}}' | \
		awk 'BEGIN{{FS="_";OFS="\t"}}{{print $1,$2,$3,$4}}' > analysis/qtl/summary/snp.dist.tmp
		awk '{{print $4 "\t" $2 "\t" $3}}' analysis/qtl/summary/snp.dist.tmp | \
		bedtools intersect -u -v - -b analysis/qtl/cis/all.samples.vcf.gz > analysis/qtl/summary/cpeaks.no.snp.bed
		Rscript scripts/snp.qtl.dist.R analysis/qtl/summary/snp.dist.tmp """

rule qtl_ai_overlap:
	"""test overlap of qtl peaks and allelically imbalanced peaks"""
	input: 
		qtl = "analysis/qtl/summary/peaks/qtl.in.peaks.sig.bed",
		ai = "analysis/summary/bed/sig.BOTH.peak.bed"
	output:
		"analysis/qtl/summary/qtl_ai_overlap.txt"
	shell:
		""" wc -l {input.ai} > {output}
		wc -l {input.qtl} >> {output}
		bedtools intersect -u -a {input.ai} -b {input.qtl} | wc -l >> {output} """

rule qtl_snp_in_peak:
	"""get qtl p-values for snps and the peaks they fall within"""
	input:
		"analysis/qtl/cis/nominals.tsv.gz"
	output:
		"analysis/qtl/summary/qtl.snps.in.peaks.txt"
	shell:
		""" printf 'peak\tstart\tend\tchr\tsnp.pos\tpval\tbeta\ttop\n' > {output}
		zcat {input} | sed '1d' | \
			awk 'BEGIN{{OFS="\t"}}{{print $1, $2, $10, $12, $13, $14}}' | \
			awk 'BEGIN{{FS="_";OFS="\t"}}{{print $1,$2,$3}}' | \
			awk 'BEGIN{{OFS="\t"}}{{if($5 >= $2 && $5 <= $3) {{print $0}}}}' >> \
			{output} """

rule qtl_snp_in_peak_fdr:
	"""get q-values for qtl snps in peaks"""
	input:
		"analysis/qtl/summary/qtl.snps.in.peaks.txt"
	output:
		snp = "analysis/qtl/summary/snps/qtl.in.peaks.sig.txt",
		bed ="analysis/qtl/summary/peaks/qtl.in.peaks.sig.bed"
	shell:
		" Rscript scripts/qtl.in.peaks.fdr.R {input} {output.snp} {output.bed} "

rule beta_vs_mu:
	"""compare beta from qtls to mu from allelically imbalanced peaks"""
	input:
		ai = "analysis/summary/sig.BOTH.txt",
		beta = "analysis/qtl/summary/qtl.snps.in.peaks.txt"
	output:
		plot = "analysis/qtl/summary/beta_vs_mu.pdf"
	shell:
		" Rscript scripts/mu_vs_beta.R {input.ai} {input.beta} {output.plot} "

rule combined_test:
	"""combine qtl and ai p-values into one significance test"""
	input:
		ai = "analysis/stratasrun/results/results.all.txt",
		qtl = "analysis/qtl/summary/qtl.snps.in.peaks.txt"
	output:
		comb = "analysis/qtl/summary/snps/combined.txt",
		snp = "analysis/qtl/summary/snps/combined.sig.txt",
		bed = "analysis/qtl/summary/peaks/combined.sig.bed"
	shell:
		" Rscript scripts/combined.test.R {input.ai} {input.qtl} {output.comb} {output.snp} {output.bed} "

# --- rules for testing genetically determined chromatin for eQTL enrichment

#test groups of peaks for GWAS SNP enrichment:
rule enrich_fg:
	input:
		fg="analysis/qtl/summary/peaks/combined.sig.bed",
		targets = "{features}.bed"
	output:
		"analysis/enrich/results/{features}/fg.enrich.txt"
	shell:
		" sh scripts/enrich.sh {input.fg} {input.targets} {output} " 


rule enrich_permute:
	"""test groups of peaks for eQTL SNP enrichment"""
	input:
		bg = "analysis/consensus_peaks/consensus.peaks.stratas.tsv",
		fg = "analysis/qtl/summary/peaks/combined.sig.bed",
		targets = "{features}.bed"
	output:
		"analysis/enrich/perm/{features}/group.{i}.bg.enrich"
	shell:
		" sh scripts/enrich.permute.sh {wildcards.i} {input.bg} {input.fg} {input.targets} analysis/enrich/perm/{wildcards.features} " 

rule enrich_random_background:
	"""test random genomic intervals for eQTL SNP enrichment"""
	input:
		fg = "analysis/qtl/summary/peaks/combined.sig.bed",
		targets = "{features}.bed"
	output:
		"analysis/enrich/perm/{features}/rand/group.{i}.bg.enrich"
	shell:
		""" # concat multiple shuffle results to provide a large sample size for enrich.perm.sh to choose from:
		printf "" > analysis/enrich/perm/{wildcards.features}/rand/tmp.{wildcards.i}
		for i in $(seq 1 50); do
			bedtools shuffle -seed $i -chrom -noOverlapping -i {input.fg} -g config["chr_sizes"] >> analysis/enrich/perm/{wildcards.features}/rand/tmp.{wildcards.i}
		done
		sh scripts/enrich.permute.sh {wildcards.i} analysis/enrich/perm/{wildcards.features}/rand/tmp.{wildcards.i} {input.fg} {input.targets} analysis/enrich/perm/{wildcards.features}/rand
		rm analysis/enrich/perm/{wildcards.features}/rand/tmp.{wildcards.i} """

rule enrich_stats:
	"""summarize eQTL enrichment stats"""
	input:
		perm = ["analysis/enrich/perm/{features}/group.%d.bg.enrich" % s for s in range(1,51)],
		perm_rand = ["analysis/enrich/perm/{features}/rand/group.%d.bg.enrich" % s for s in range(1,51)],	
		enrich = "analysis/enrich/results/{features}/fg.enrich.txt"
	output:
		"analysis/enrich/results/{features}/enrichment.stats.txt"
	shell:
		""" fg_enrich=`sed '1d' {input.enrich} | cut -f4`
		hits=`cat {input.perm} | \
			awk -v n=$fg_enrich 'BEGIN{{count=0}} $1>n{{count++}}END{{print count}}'`
		tot=`cat {input.perm} | wc -l`
		bg_enrich=`cat {input.perm} | awk -f scripts/avg.awk`
		rel_enrich=`echo $fg_enrich $bg_enrich | awk '{{print $1/$2}}'`
		pval=`echo $hits $tot | awk '{{print $1/$2}}'`
		hits_rand=`cat {input.perm_rand} | \
                        awk -v n=$fg_enrich 'BEGIN{{count=0}} $1>n{{count++}}END{{print count}}'`
                tot_rand=`cat {input.perm_rand} | wc -l`                                                           
                bg_enrich_rand=`cat {input.perm_rand} | awk -f scripts/avg.awk`                                    
                rel_enrich_rand=`echo $fg_enrich $bg_enrich_rand | awk '{{print $1/$2}}'`                          
                pval_rand=`echo $hits_rand $tot_rand | awk '{{print $1/$2}}'` 

		printf "bg_enrich\trel_enrich\tpval\tbg_enrich_rand\trel_enrich_rand\tpval_rand\n%s\t%s\t%s\t%s\t%s\t%s\n" $bg_enrich $rel_enrich $pval $bg_enrich_rand $rel_enrich_rand $pval_rand| \
			paste {input.enrich} - > {output}"""

rule plot_enrichment:
	"""plot eQTL enrichment stats"""
	input:
		"analysis/enrich/results/qtl/Prostate.qtl.hg19/enrichment.stats.txt",
		"analysis/enrich/results/qtl/Whole_Blood.qtl.hg19/enrichment.stats.txt",
		"analysis/enrich/results/qtl/Pancreas.qtl.hg19/enrichment.stats.txt",
		"analysis/enrich/results/qtl/Esophagus_Mucosa.qtl.hg19/enrichment.stats.txt",
		"analysis/enrich/results/qtl/Lung.qtl.hg19/enrichment.stats.txt",
		"analysis/enrich/results/qtl/Liver.qtl.hg19/enrichment.stats.txt",
		"analysis/enrich/results/qtl/Brain_Cortex.qtl.hg19/enrichment.stats.txt",
		"analysis/enrich/results/qtl/Thyroid.qtl.hg19/enrichment.stats.txt",
		"analysis/enrich/results/qtl/Stomach.qtl.hg19/enrichment.stats.txt",
		"analysis/enrich/results/qtl/Spleen.qtl.hg19/enrichment.stats.txt",
		"analysis/enrich/results/qtl/Breast_Mammary_Tissue.qtl.hg19/enrichment.stats.txt",
		"analysis/enrich/results/qtl/Artery_Aorta.qtl.hg19/enrichment.stats.txt",
		"analysis/enrich/results/qtl/Kidney_Cortex.qtl.hg19/enrichment.stats.txt",
		"analysis/enrich/results/qtl/Skin_Sun_Exposed_Lower_leg.qtl.hg19/enrichment.stats.txt"
	output:
		summary = "analysis/enrich/results/enrichment.summary.txt",
		plot = "analysis/enrich/results/enrichment.pdf"
	shell:
		"""grep -v fg {input} | sed 's/.*results\///' | sed 's/\/enrichment.stats.txt:/\t/' > {output.summary}
		Rscript scripts/plot_enrichment.R {output.summary} {output.plot} """

# --- rules for motif analysis in allelically imbalanced peaks

rule peak_motifs:
	"""find top de novo motifs with homer"""
	input:
		"analysis/consensus_peaks/consensus.peaks.stratas.tsv"		
	output:
		expand("analysis/motifs/homer/homerResults/motif{i}.motif", i=range(1,7))
	shell:
		""" sed '1d' {input} | shuf -n 10000 --random-source={input} > analysis/motifs/tmp.motif.bed 
		findMotifsGenome.pl analysis/motifs/tmp.motif.bed /home/scb20/ref/hg19/Homo_sapiens_assembly19.fasta analysis/motifs/homer -size given -p 4 -mask -seqlogo -preparsedDir analysis/motifs/homer 
		rm analysis/motifs/tmp.motif.bed"""


rule ai_vs_motifs:
	"""plot allele fraction vs motif score"""
	input:
		"analysis/motifs/homer/homerResults/motif{i}.motif",
		ai = "analysis/summary/sig.BOTH.txt",
		counts = "analysis/stratasrun/all.counts.mat"
	output:
		"analysis/motifs/ai_vs_motifs/motif{i}.ai.pdf"
	params:
		genome = config["genome"]
	shell:
		""" sh scripts/motifs.sh {input.ai} {input.counts} 20 analysis/motifs/ai_vs_motifs/snps_in_motifs/motif{wildcards.i}/ {params.genome} analysis/motifs/homer/homerResults/motif{wildcards.i}.motif
		Rscript scripts/motif.vs.AF.R analysis/motifs/ai_vs_motifs/snps_in_motifs/motif{wildcards.i}/motifs.by.allele {input.ai} {output}"""

rule control_motif:
	"""plot allele fraction vs motif score for control motif"""
	input:
		config["control_motif"]
	output:
		"analysis/motifs/homer/homerResults/motifX.motif"
	shell:
		"cp {input} {output}"

# --- rules for building genetic models of chromatin

rule prep_stratas_predict:
	"""prepare inputs for stratAS model building step"""
	input:
		mat = "analysis/qtl/quan/rpkm.merged.bed.gz",
		covar = "analysis/qtl/quan/covariate.pca"
		
	output:
		mat = "analysis/predict/total.mat",
		covar = "analysis/predict/covar.mat"
	shell:
		""" zcat {input.mat} | cut -f1,4 | \
		  awk 'BEGIN{{FS="_";OFS="\t"}}{{print $1, $2, $3}}' | \
		  awk 'BEGIN{{FS="\t"}}{{print "chr" $1 ":" $3 "-" $4}}' > analysis/predict/peaks.tmp
		zcat {input.mat} | cut -f7- > analysis/predict/vals.tmp
		paste analysis/predict/peaks.tmp analysis/predict/vals.tmp | sed 's/chr#chr:-\t//' > {output.mat}

		# keep only first 6 covariates:
		head -n 7 {input.covar} > {output.covar}

		rm analysis/predict/vals.tmp analysis/predict/peaks.tmp """

rule run_stratas_predict:
	"""generate chromatin ~ SNP models with stratAS"""
	input:
		peaks = "analysis/consensus_peaks/consensus.peaks.stratas.tsv",
		local_params = "analysis/stratasrun/all.local.params",
		global_params = "analysis/stratasrun/all.global.params",
		split = "analysis/stratasrun/split/counts.{split}",
		mat = "analysis/predict/total.mat",
		covar = "analysis/predict/covar.mat"
	output:
		"analysis/predict/results/predict.{split}"
	params:
		samples = config["samples"]
	shell:
		""" 
		mkdir -p analysis/predict/WEIGHTS
		Rscript stratAS/stratas.R \
			--input {input.split} \
			--samples {params.samples} \
			--peaks {input.peaks} \
			--global_param {input.global_params} \
			--local_param {input.local_params} \
			--max_rho 0.2 \
			--min_cov 1 \
			--window 25000 \
			--predict_snps "LDREF/hm3.pos" \
			--total_matrix {input.mat} \
			--covar {input.covar} \
			--predict --predict_only \
			--weights_out "analysis/predict/WEIGHTS" \
			--seed 123 > analysis/predict/results/predict.{wildcards.split}
		"""

rule merge_predict_results:
	"""merge chromatin models"""
	input:
		["analysis/predict/results/predict.%02d" % s for s in range(0,50)]
	output:
		"analysis/predict/results/predict.all.txt"
	shell:
		" cat misc/predict.head {input} > {output} "

rule hsq:
	"""calculate cis-SNP heritability for chromatin QTLs"""
	input:
		vcf = "analysis/qtl/cis/all.samples.vcf.gz",
		cov = "analysis/predict/covar.mat",
		rpkm = "analysis/qtl/quan/rpkm.merged.bed.gz",
	output:
		"analysis/hsq/hsq/{chr}.hsq"
	params:
		win = "5e5",
		LDREF = "LDREF",
		hsq = "hsq"
	shell:
		"sh hsq/compute_weights.sh {input.vcf} {input.cov} {input.rpkm} {params.LDREF} {wildcards.chr} {params.win} analysis/hsq {params.hsq}"

rule merge_hsq:
	"""merge cis-SNP heritability results for chromatin QTLs"""
	input:
		["analysis/hsq/hsq/%d.hsq" % s for s in range(1,22)]
	output:
		summary = "analysis/hsq/hsq.summary.txt"
	shell:
		"cat {input} > {output.summary} "

# --- rules for testing chromatin model - trait associations

rule prep_fusion:
	"""prepare input files for FUSION"""
	input:
		"analysis/predict/results/predict.all.txt"
	output:
		"analysis/fusion/weights.pos"
	shell:
		""" ls analysis/predict/WEIGHTS/ | grep RDat > analysis/fusion/weights.tmp
			sed 's/chr//' analysis/fusion/weights.tmp | \
			sed 's/.wgt.RDat//' | sed 's/:/\t/' | \
			sed 's/-/\t/' > analysis/fusion/coords.tmp
		paste analysis/fusion/weights.tmp \
			<(sed 's/.wgt.RDat//' analysis/fusion/weights.tmp) \
			analysis/fusion/coords.tmp > analysis/fusion/comb.tmp
		printf "WGT\tID\tCHR\tP0\tP1\n" > {output}
		cat analysis/fusion/comb.tmp >> {output}
		rm analysis/fusion/coords.tmp analysis/fusion/weights.tmp analysis/fusion/comb.tmp """

rule prep_gwas_bed:
	"""create bed files with gwas-significant SNPs and their regions"""
	input:
		"gwas_data/{gwas}.sumstats.gz"
	output:
		gwas_snps = "analysis/fusion/{gwas}/sig.snps.bed",
		gwas_bed = "analysis/fusion/{gwas}/gwas.bed"
	shell:
		""" Rscript scripts/sumstats.to.bed.R {input} {output.gwas_snps}
		bedtools slop -i {output.gwas_snps} -b 1000000 -g misc/hg19.genome | sort -k1,1 -k2,2n | bedtools merge -i - > {output.gwas_bed} """ 

rule fusion:
	"""run FUSION (chromatin model - GWAS association)"""
	input:
		weights = "analysis/fusion/weights.pos",
		sumstats = "gwas_data/{gwas}.sumstats.gz"	
	output:
		"analysis/fusion/{gwas}/meta.lasso.{chr}.dat",
		"analysis/fusion/{gwas}/meta.lasso.as.{chr}.dat",
		"analysis/fusion/{gwas}/meta.lasso.plasma.{chr}.dat",
		"analysis/fusion/{gwas}/meta.top1.as.{chr}.dat",
		"analysis/fusion/{gwas}/meta.top1.qtl.{chr}.dat",
		"analysis/fusion/{gwas}/meta.top1.{chr}.dat"
	shell:
		"""for mod in lasso lasso.as lasso.plasma top1.as top1.qtl top1
		do
			Rscript fusion/FUSION.assoc_test.R \
			--sumstats {input.sumstats} \
			--out analysis/fusion/{wildcards.gwas}/meta.$mod.{wildcards.chr}.dat \
			--weights {input.weights} \
			--weights_dir analysis/predict/WEIGHTS/ \
			--ref_ld_chr ./LDREF/1000G.EUR. \
			--chr {wildcards.chr} \
			--force_model $mod
		done """

rule merge_fusion:
	"""merge results from FUSION"""
	input:
		lasso = ["analysis/fusion/{gwas}/meta.lasso.%d.dat" % s for s in range(1,23)],
		lasso_as = ["analysis/fusion/{gwas}/meta.lasso.as.%d.dat" % s for s in range(1,23)],
		lasso_plasma = ["analysis/fusion/{gwas}/meta.lasso.plasma.%d.dat" % s for s in range(1,23)],
		top1_as = ["analysis/fusion/{gwas}/meta.top1.as.%d.dat" % s for s in range(1,23)],
		top1_qtl = ["analysis/fusion/{gwas}/meta.top1.qtl.%d.dat" % s for s in range(1,23)],
		top1 = ["analysis/fusion/{gwas}/meta.top1.%d.dat" % s for s in range(1,23)]
	output:
		lasso = "analysis/fusion/{gwas}/merged/fusion.lasso.txt",
		lasso_as = "analysis/fusion/{gwas}/merged/fusion.lasso.as.txt",
		lasso_plasma = "analysis/fusion/{gwas}/merged/fusion.lasso.plasma.txt",
		top1_as = "analysis/fusion/{gwas}/merged/fusion.top1.as.txt",
		top1_qtl = "analysis/fusion/{gwas}/merged/fusion.top1.qtl.txt",
		top1 = "analysis/fusion/{gwas}/merged/fusion.top1.txt"
	shell:
		""" HEADER=analysis/fusion/{wildcards.gwas}/merged/header.tmp
		head -n 1 $(echo {input} | sed 's/ .*//') > $HEADER
		cp $HEADER {output.lasso}
		cp $HEADER {output.lasso_as}
		cp $HEADER {output.lasso_plasma}
		cp $HEADER {output.top1_as}
		cp $HEADER {output.top1_qtl}
		cp $HEADER {output.top1}
		cat {input.lasso} | grep -v PANEL >> {output.lasso}
		cat {input.lasso_as} | grep -v PANEL >> {output.lasso_as}
		cat {input.lasso_plasma} | grep -v PANEL >> {output.lasso_plasma}
		cat {input.top1_as} | grep -v PANEL >> {output.top1_as}
		cat {input.top1_qtl} | grep -v PANEL >> {output.top1_qtl}
		cat {input.top1} | grep -v PANEL >> {output.top1} 
		rm $HEADER """

rule fusion_sig:
	"""prune chromatin-trait associations by significance"""
	input:
		"analysis/fusion/{gwas}/merged/fusion.lasso.txt",
		"analysis/fusion/{gwas}/merged/fusion.lasso.as.txt",
		"analysis/fusion/{gwas}/merged/fusion.lasso.plasma.txt",
		"analysis/fusion/{gwas}/merged/fusion.top1.as.txt",
		"analysis/fusion/{gwas}/merged/fusion.top1.qtl.txt",
		"analysis/fusion/{gwas}/merged/fusion.top1.txt",
		"analysis/fusion/{gwas}/gwas.bed"
	output:
		"analysis/fusion/{gwas}/merged/sig/cwas.sig.lenient.txt",
		"analysis/fusion/{gwas}/merged/sig/cwas.sig.stringent.txt",
		"analysis/fusion/{gwas}/merged/sig/cwas.sig.best.txt",
		"analysis/fusion/{gwas}/merged/sig/cwas.sig.gwas.txt",
		"analysis/fusion/{gwas}/merged/sig/cwas.sig.gwas.restricted.txt",
                "analysis/fusion/{gwas}/merged/cv/cwas.cv.lenient.txt",
                "analysis/fusion/{gwas}/merged/cv/cwas.cv.stringent.txt",
                "analysis/fusion/{gwas}/merged/cv/cwas.cv.best.txt",
                "analysis/fusion/{gwas}/merged/cv/cwas.cv.gwas.txt",
                "analysis/fusion/{gwas}/merged/cv/cwas.cv.gwas.restricted.txt"
	params:
		path = "analysis/fusion/{gwas}/merged"
	shell:
		"Rscript scripts/sig.models.R {params.path}"

rule manhattan_plot:
	"""create manhattan plots of chromatin-trait associations"""
	input:
		"analysis/fusion/{gwas}/merged/sig/cwas.sig.lenient.txt",
		"analysis/fusion/{gwas}/merged/sig/cwas.sig.stringent.txt",
		"analysis/fusion/{gwas}/merged/sig/cwas.sig.best.txt",
		"analysis/fusion/{gwas}/merged/sig/cwas.sig.gwas.txt",
		"analysis/fusion/{gwas}/merged/cv/cwas.cv.lenient.txt",
		"analysis/fusion/{gwas}/merged/cv/cwas.cv.stringent.txt",
		"analysis/fusion/{gwas}/merged/cv/cwas.cv.gwas.txt"
	output:
		manhattan = "analysis/fusion/{gwas}/plots/manhattan.pdf",
	params:
		path = "analysis/fusion/{gwas}/merged",
		zzplot_base = "analysis/fusion/{gwas}/plots/zzplot"
	shell:
		"Rscript scripts/manhattan.R {params.path} {output.manhattan} {params.zzplot_base}"

rule annotate_cwas:
	"""summarize and annotate significant cwas results"""
	input:
		cwas_tested = "analysis/fusion/{gwas}/merged/cv/cwas.cv.best.txt", 
		cwas_sig = "analysis/fusion/{gwas}/merged/sig/cwas.sig.best.txt",
		twas = "twas_data/ONCOARRAY.TWAS.txt"
	output:
		summary = "analysis/fusion/{gwas}/merged/summary.txt"
	params:
		tss = config["tss"],
		hichip = config["HiChIP"],
		peaks1 = config["peaks1"],
		peaks2 = config["peaks2"]
	shell:
		"Rscript scripts/annotate_cwas.R {input.cwas_tested} {input.cwas_sig} {input.twas} {output.summary} {params.tss} {params.hichip} {params.peaks1} {params.peaks2}"

rule fusion_post_process:
	"""run post-processing FUSION steps, including conditioning analysis"""
	input:
		"analysis/fusion/{gwas}/merged/sig/cwas.sig.{strategy}.txt"
	output:
		"analysis/fusion/{gwas}/postprocess/{strategy}/done.tmp"
	shell:
		""" for chr in {{1..22}}; do
			echo chromosome $chr
			Rscript fusion/FUSION.post_process.R \
				--input analysis/fusion/{wildcards.gwas}/merged/sig/cwas.sig.{wildcards.strategy}.txt \
				--out analysis/fusion/{wildcards.gwas}/postprocess/{wildcards.strategy}/chr${{chr}} \
				--sumstats gwas_data/{wildcards.gwas}.sumstats.gz \
				--ref_ld_chr ./LDREF/1000G.EUR. \
				--plot TRUE \
				--plot_eqtl TRUE \
				--plot_corr TRUE \
				--eqtl_model top1.qtl \
				--save_loci TRUE \
				--report TRUE \
				--plot_individual TRUE \
				--chr $chr
			done
		touch {output} """

rule summarize_conditioning:
	"""summarize FUSION conditioning analysis"""
	input:
		"analysis/fusion/{gwas}/postprocess/{strategy}/done.tmp"
	output:
		"analysis/fusion/{gwas}/postprocess/{strategy}/conditioned.summary.txt"
	shell:
		"sh scripts/summarize_conditioning.sh analysis/fusion/{wildcards.gwas}/postprocess/{wildcards.strategy} {output}"

rule compare_gwas:
	"""compare z-scores for chromatin associations with a variety of gwas traits"""
	input:
                "analysis/fusion/ProstateCancer_Meta_Schumacher2018.nodup/merged/cv/cwas.cv.best.txt",
                "analysis/fusion/UKB_460K.biochemistry_Testosterone_Male/merged/cv/cwas.cv.best.txt",
                "analysis/fusion/UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED/merged/cv/cwas.cv.best.txt",
                "analysis/fusion/UKB_460K.disease_PSORIASIS/merged/cv/cwas.cv.best.txt",
                "analysis/fusion/UKB_460K.disease_DERMATOLOGY/merged/cv/cwas.cv.best.txt",
                "analysis/fusion/UKB_460K.disease_AID_ALL/merged/cv/cwas.cv.best.txt",
                "analysis/fusion/UKB_460K.disease_ASTHMA_DIAGNOSED/merged/cv/cwas.cv.best.txt",
                "analysis/fusion/UKB_460K.disease_CARDIOVASCULAR/merged/cv/cwas.cv.best.txt",
                "analysis/fusion/UKB_460K.disease_HI_CHOL_SELF_REP/merged/cv/cwas.cv.best.txt",
                "analysis/fusion/UKB_460K.disease_HYPERTENSION_DIAGNOSED/merged/cv/cwas.cv.best.txt",
                "analysis/fusion/UKB_460K.disease_RESPIRATORY_ENT/merged/cv/cwas.cv.best.txt",
                "analysis/fusion/UKB_460K.disease_T2D/merged/cv/cwas.cv.best.txt",
		"analysis/fusion/UKB_460K.disease_HYPOTHYROIDISM_SELF_REP/merged/cv/cwas.cv.best.txt",
		"analysis/fusion/UKB_460K.biochemistry_Testosterone_Male/merged/cv/cwas.cv.best.txt",
		"analysis/fusion/UKB_460K.blood_EOSINOPHIL_COUNT/merged/cv/cwas.cv.best.txt",
		"analysis/fusion/UKB_460K.blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT/merged/cv/cwas.cv.best.txt",
		"analysis/fusion/UKB_460K.blood_LYMPHOCYTE_COUNT/merged/cv/cwas.cv.best.txt",
		"analysis/fusion/UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN/merged/cv/cwas.cv.best.txt",
		"analysis/fusion/UKB_460K.blood_MEAN_PLATELET_VOL/merged/cv/cwas.cv.best.txt",
		"analysis/fusion/UKB_460K.blood_MEAN_SPHERED_CELL_VOL/merged/cv/cwas.cv.best.txt",
		"analysis/fusion/UKB_460K.blood_MONOCYTE_COUNT/merged/cv/cwas.cv.best.txt",
		"analysis/fusion/UKB_460K.blood_PLATELET_COUNT/merged/cv/cwas.cv.best.txt",
		"analysis/fusion/UKB_460K.blood_PLATELET_DISTRIB_WIDTH/merged/cv/cwas.cv.best.txt",
		"analysis/fusion/UKB_460K.blood_RBC_DISTRIB_WIDTH/merged/cv/cwas.cv.best.txt",
		"analysis/fusion/UKB_460K.blood_RED_COUNT/merged/cv/cwas.cv.best.txt",
		"analysis/fusion/UKB_460K.blood_WHITE_COUNT/merged/cv/cwas.cv.best.txt"
	output:
		"analysis/fusion/compare_gwas/compare.gwas.pdf"
	shell:
		"Rscript scripts/compare.gwas.R {input}"
