#!/usr/bin/env -S snakemake --cores 1 --quiet --snakefile

import tempfile
import csv
from pathlib import Path
import glob

# Site-specific config, most changes would be made in that file rather than this one.
configfile: "config.yaml"

analysis_dir = config["analysis_dir"]
bcl_dir = config["bcl_dir"]
output_dir = config["output_dir"]
samplesheet_dir = config["samplesheet_dir"]
temp_dir = config["temp_dir"]
sequencer = config["sequencer"]
refgenome = config["refgenome"]
ref_exon_annotations = config["ref_exon_annotations"]


# Convenience method to strip lines starting with a hashmark from a CSV file being read (assumed to be comments).
def decomment(csvfile):
    for row in csvfile:
        raw = row.split('#')[0].strip()
        if raw: yield raw

# The format of the paired DNA samples file is tab-delimited, with subject ID in the zeroth column, tumor sample name in first column, then matched normal in the second. 
# Additional columns are metadata: Boolean indicating if some tumor is expected in the normal sample, and a number representing the tumor site as per PCGR guidelines (see generate_pcgr.sh for details).
dna_paired_samples = {}
with open(config["dna_paired_samples_tsv"], 'r') as data_in:
	tsv_file = csv.reader(decomment(data_in), delimiter="\t")
	for line in tsv_file:
		if len(line[0]) > 35 or len(line[0]) < 6:
			raise NameError("Subject names are limited to between 6 and 35 characters by PCGR reporting requirements, please revise "+line[0])
		dna_paired_samples_key = "_".join((line[0],line[1],line[2]))
		dna_paired_samples[dna_paired_samples_key] = [line[0],line[3],line[4]]

# The format of the paired RNA samples file is tab-delimited, with tumor RNA sample name in first column, then matched normal *DNA* in the second. Additional columns are ignored.
rna_paired_samples = []
#with open(config["rna_paired_samples_tsv"], 'r') as data_in:
#	tsv_file = csv.reader(data_in, delimiter="\t")
#	for line in tsv_file:
#		rna_paired_samples.append("-".join((line[0],line[1],line[2])))

rule all:
	# Generate a Dragen small nucleotide variant somatic mutation analysis for each pair for tumor-normal listed in the pairs file.
	# Also generate a PCGR user-friendly variant report, which will require the germline variant calls as well in antoher rule.
	input:
		expand("{output_dir}/{sample_combo}.dna.somatic.hard-filtered.vcf.gz", output_dir=output_dir, sample_combo=list(dna_paired_samples.keys())),
		expand("{output_dir}/{sample_combo}.dna.somatic.cnv.vcf.gz", output_dir=output_dir, sample_combo=list(dna_paired_samples.keys())),
		expand("{output_dir}/{sample_combo}.dna.somatic.sv.vcf.gz", output_dir=output_dir, sample_combo=list(dna_paired_samples.keys())),
		expand(output_dir+"/pcgr/{sample_combo}/{subject}.pcgr_acmg.grch38.flexdb.html", zip, subject=[tuple[0] for tuple in dna_paired_samples.values()], sample_combo=list(dna_paired_samples.keys()))
		#expand("{output_dir}/{sample_combo}.rna.quant.sf", output_dir=output_dir, sample_combo=rna_paired_samples)
		#expand("{output_dir}/{sample_combo}.rna.fusion_candidates.features.csv", output_dir=output_dir, sample_combo=rna_paired_samples)

def get_samplesheet(wildcards):
	return samplesheet_dir+'/'+wildcards.run+'.csv'

def get_tumor_site(wildcards):
	return dna_paired_samples["_".join((wildcards.subject, wildcards.tumor, wildcards.normal))][2]

def get_file(sample_name, suffix):
	# Code that returns a list of bam files for based on *sample_name* 
	bams = sorted(glob.glob('{analysis_dir}/secondary/{sequencer}/*/'+sample_name+suffix))
	# If the file with the given suffix hasn't been generated yet, find which sequencing run has it
	if not bams:
		for path in Path(samplesheet_dir).iterdir():
			if path.is_file() and Path(path).suffix == '.csv':
				# print("Reading file " + str(path))
				current_file = open(path, "r")
				csv_file = csv.reader(current_file)
				sample_name_index = -1
				sample_id_index = -1
                        	for line in csv_file:
					if 'Sample_Name' in line and 'Sample_ID' in line: #it's the header for the sample list
                                        	sample_name_index = line.index('Sample_Name')
                                        	sample_id_index = line.index('Sample_ID')
					elif sample_name_index != -1 and sample_id_index != -1 and line[sample_name_index] == sample_name and (len(line) < 11 or not "RNA" in line[10]):
                                        	return analysis_dir+'/secondary/'+sequencer+'/'+os.path.splitext(os.path.basename(path))[0]+'/'+sample_name+suffix;
				if sample_name_index == -1:
					print("Missing Sample_Name column in Illumina samplesheet "+str(path)+", skipping")
				if sample_id_index == -1:
					print("Missing Sample_ID column in Illumina samplesheet "+str(path)+", skipping")
			
def get_normal_bam(wildcards):
	return get_file(wildcards.normal, ".bam")

def get_tumor_bam(wildcards):
	return get_file(wildcards.tumor, ".bam")

def get_tumor_rna_bam(wildcards):
	return get_file(wildcards.tumor, ".rna_spliced.bam")

# Germline cancer susceptibility reporting
rule generate_cpsr_json:
	input:
		germline_snv_vcf = output_dir+'/{subject}_{normal}.dna.germline.hard-filtered.vcf.gz'
	output:
		output_dir+'/{subject}_{normal}.cpsr.grch38.json.gz'
	shell:
		"cpsr --input_vcf {input.germline_snv_vcf} --pcgr_dir . --output_dir {output_dir} --genome_assembly grch38 --panel_id 0 --sample_id {wildcards.subject}_{wildcards.normal} --secondary_findings --classify_all --maf_upper_threshold 0.2 --force_overwrite"

# Personal Cancer Genome Report (user-friendly triaged variants self-contained Web page)
# The Web page name is limited by PCGR to 35 characters, so we have the full subject-tumor-normal unique combo in the dir name only.
rule generate_pcgr_html:
	input:
		cpsr_json = output_dir+'/{subject}_{normal}.cpsr.grch38.json.gz',
		somatic_snv_vcf = output_dir+'/{subject}_{tumor}_{normal}.dna.somatic.hard-filtered.vcf.gz',
		somatic_cnv_vcf = output_dir+'/{subject}_{tumor}_{normal}.dna.somatic.cnv.vcf.gz'
	output:
		output_dir+'/pcgr/{subject}_{tumor}_{normal}/{subject}.pcgr_acmg.grch38.flexdb.html'
	run:
		tumor_site = get_tumor_site(wildcards)

		# Somatic small nucleotide variants reformatting
		tf=tempfile.NamedTemporaryFile(suffix=".vcf")
		SNVFILE=tf.name
		shell("gzip -cd {input.somatic_snv_vcf} | perl -pe 'if(/^##INFO=<ID=DP,/){{print \"##INFO=<ID=TDP,Number=1,Type=Integer,Description=\\\"Read depth of alternative allele in the tumor\\\">\\n##INFO=<ID=TVAF,Number=1,Type=Float,Description=\\\"Alternative allele proportion of reads in the tumor\\\">\\n\"}}($tdp, $tvaf) = /\\t[01][\/|][01]:\d+.?\d*:\d+,(\d+):([0-9]+\.[0-9]*):\S+?$/;s/\\tDP=/\\tTDP=$tdp;TVAF=$tvaf;DP=/; s/;SOMATIC//' > {SNVFILE}")
		# Only keeping the original to avoid FileNotFoundError when temp file automatically cleaned up by Snakemake after rule application.
		shell("bgzip -c {SNVFILE} > {SNVFILE}.gz")
		shell("tabix {SNVFILE}.gz")

		# Somatic copy number variants reformatting
		tf2=tempfile.NamedTemporaryFile(suffix=".txt")
		CNAFILE=tf2.name
		shell("gzip -cd {input.somatic_cnv_vcf} | perl -ane 'if($.==1){{print \"Chromosome\\tStart\\tEnd\\tSegment_Mean\\n\"}}next if /^#/; ($end) = /END=(\d+)/; @d = split /:/, $F[$#F]; print join(\"\\t\", $F[0], $F[1], $end, $d[5]-2),\"\\n\"' > {CNAFILE}")

		# Generate PCGR report with all these data
		shell("pcgr --pcgr_dir . --output_dir {output_dir}/pcgr/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal} --sample_id {wildcards.subject} --debug --tumor_dp_tag TDP --tumor_af_tag TVAF --genome_assembly grch38 --input_vcf {SNVFILE}.gz --tumor_site {tumor_site} --tumor_purity 0.9 --tumor_ploidy 2.0 --include_trials --assay WGS --estimate_signatures --estimate_msi_status --estimate_tmb --input_cna {CNAFILE} --cpsr_report {input.cpsr_json} --show_noncoding")

# Used by PCGR for reporting known germline cancer susceptibility or related variants.
# The germline ("Normal") BAM file was generated by the somatic variant Dragen call.
rule dragen_germline_snv_vcf:
	input:
		germline_bam = get_normal_bam
	output:
		output_dir+'/{subject}_{normal}.dna.germline.hard-filtered.vcf.gz'
	shell:
		"dragen -b {input.germline_bam} -r {refgenome} --output-directory {output_dir} --output-file-prefix {wildcards.subject}_{wildcards.normal}.dna.germline --enable-variant-caller true --intermediate-results-dir {temp_dir} --enable-map-align false"

rule dragen_bcl_conversion:
	# If the run was generated using a library prep kit that includes UMIs, this must be noted in the samplesheet to allow proper grouping of the reads.
	# See config.yaml for more details.
	input:
		csv=get_samplesheet
	output:
		analysis_dir+'/primary/{sequencer}/{run}/Reports/fastq_list.csv'
	shell:
		# Nota bene: Disabled for the moment to avoid big dragen job dependency cascades
		"echo bcl-convert --force --sample-sheet {input.csv} --bcl-input-directory {bcl_dir}/{sequencer}/{wildcards.run} --output-directory {analysis_dir}/primary/{sequencer}/{wildcards.run}"

# One or more libraries correspond to a single sample defined by the wildcard match values for a sequencing run. 
# Generate a file listing all those FASTQ files so they can be processed together, e.g. for reference mapping.
def make_sample_fastq_list_csv(fastq_list_csv, wildcards, sample_libraries):
	# We could generate a combined FASTQ file, but this takes up a lot of extra space and tmp dir might fill for somatic genomes for example
	# Instead we will generate a fastq-list file that contains all the libraries for the sample (you may have more than one as in XP loading you can't name same sample ID across lanes for same sample name) 
	sample_fastq_list_csv = analysis_dir+'/secondary/'+wildcards.sequencer+'/'+wildcards.run+'/'+wildcards.sample+'_fastq_list.csv'
	with open(sample_fastq_list_csv, 'w') as f:
		with open(fastq_list_csv, 'r') as data_in:
			csv_file = csv.reader(data_in)
			# Print the header line
			f.write(",".join(next(csv_file)))
			f.write('\n')
       			for line in csv_file:
				if line[1] in sample_libraries:
					# We might have different sample names for the same sample like Li###-lane1
					# and this will cause problems downstream, e.g. when somatic mutation calling
					# and dragen enforces a simngle sample name in the source BAM. So replace any existing name.
					# Otherwise post hoc you'll need to do something like 
					# samtools reheader -c 'perl -ne "s/^(\@RG.*\s+SM:\S+)-lane\d/$1/"'
					line[2] = wildcards.sample
					f.write(",".join(line))
					f.write('\n')
	return sample_fastq_list_csv

# Return a tuple of two:
# Boolean as to whether this library should get the Unique Molecular Indices treatment, using info from the SampleSheet.csv
# The list of library names (that will be in the FASTQ list CSV) that correspond to the sample name identified by the wildcards 
def identify_dna_libraries(samplesheet, wildcards):
	# Figure out the library ID (which is in the FASTQ file names) from the sample name in the samplesheet
	# If XP loading was used, you can have multiple "libraries" called Li###-lane1, etc. for the same sample on a run as IEV requires unique sample IDs and normally we'd use the library ID
	has_UMIs = False
	sample_libraries = []
	sample_project_index = -1
	sample_name_index = -1
	sample_id_index = -1
	with open(samplesheet, 'r') as data_in:
		print("Gathering FASTQ file for "+wildcards.sample+" from "+samplesheet)
		csv_file = csv.reader(data_in)
		for line in csv_file:
			if "Sample_Name" in line and "Sample_ID" in line and "Sample_Project" in line: #it's the header for the sample list
				sample_name_index = line.index("Sample_Name")
				sample_id_index = line.index("Sample_ID")
				sample_project_index = line.index("Sample_Project")
			elif sample_project_index != -1 and "RNA" in line[sample_project_index]:
				continue
			elif sample_name_index != -1 and sample_id_index != -1 and line[sample_name_index] == wildcards.sample:
				sample_libraries.append(line[sample_id_index]) 
			elif len(line) > 1 and line[0] == "OverrideCycles" and "U" in line[1]:
				has_UMIs = True
		if sample_name_index == -1:
			raise NameError("Missing Sample_Name column in Illumina samplesheet "+samplesheet)
		if sample_id_index == -1:
			raise NameError("Missing Sample_ID column in Illumina samplesheet "+samplesheet)
	return has_UMIs, sample_libraries

rule dragen_dna_read_mapping:
	# This rule will check if the library was UMI-indexed or not, and change the mapping procedure accordingly.
	input:
		samplesheet_dir+'/{run}.csv',
		analysis_dir+'/primary/{sequencer}/{run}/Reports/fastq_list.csv',
		refgenome+'/reference.bin'
	output:
		analysis_dir+'/secondary/{sequencer}/{run}/{sample}.bam'
	run:
		library_info = identify_dna_libraries(input[0], wildcards)
		has_UMIs = library_info[0]
		#has_UMIs = False
		sample_libraries = library_info[1]
		print("Using UMIs: " + str(has_UMIs))
		print("Libraries: " + ", ".join(sample_libraries))
		this_sample_only_fastq_list_csv = make_sample_fastq_list_csv(input[1], wildcards, sample_libraries)
		if has_UMIs:
			shell("dragen -r {refgenome} --output-dir {analysis_dir}/secondary/{wildcards.sequencer}/{wildcards.run} --output-file-prefix {wildcards.sample} --enable-map-align true --enable-sort true --enable-map-align-output true --enable-bam-indexing true --umi-enable true --umi-correction-scheme=random --umi-min-supporting-reads 1 --umi-min-map-quality 1 --fastq-list {this_sample_only_fastq_list_csv} --fastq-list-all-samples true --intermediate-results-dir {temp_dir}")
		else:
			shell("dragen -r {refgenome} --output-dir {analysis_dir}/secondary/{wildcards.sequencer}/{wildcards.run} --output-file-prefix {wildcards.sample} --enable-map-align true --enable-sort true --enable-map-align-output true --enable-bam-indexing true --fastq-list {this_sample_only_fastq_list_csv} --fastq-list-all-samples true --intermediate-results-dir {temp_dir}")

#rule dragen_rna_read_mapping:
#        input:
#		samplesheet_dir+'/{run}.csv',
#		analysis_dir+'/primary/{sequencer}/{run}/Reports/fastq_list.csv',
#		refgenome+'/anchored_rna',
#		ref_exon_annotations
#	output:
#		analysis_dir+'/secondary/{sequencer}/{run}/{sample}.rna_spliced.bam'
#	shell:
#		"dragen -r {refgenome} --output-dir {analysis_dir}/secondary/{wildcards.sequencer}/{wildcards.run} --output-file-prefix {wildcards.sample}.rna_spliced --enable-sort true --enable-rna true --enable-map-align true --enable-map-align-output true --enable-bam-indexing true --annotation-file {ref_exon_annotations}"

rule dragen_exec_rna_quant_fusion_variant_calls:
	input:
		tumor_rna_bam = get_tumor_rna_bam
	output:
		"{output_dir}/{subject}_{tumor}_{normal}.rna.quant.sf",
		"{output_dir}/{subject}_{tumor}_{normal}.rna.fusion_candidates.features.csv"
	shell:
		"dragen -r {refgenome} --tumor-bam-input {input.tumor_rna_bam} --enable-map-align false --output-dir {output_dir} --output-file-prefix {wildcards.subject}-{wildcards.tumor}-{wildcards.sample}.rna.somatic --enable-variant-caller true --enable-rna-quantification true --enable-rna-gene-fusion true"

rule dragen_exec_somatic_snv_sv_and_cnv_calls:
	input:
		tumor_bam = get_tumor_bam,
		germline_bam = get_normal_bam
	output:
		"{output_dir}/{subject}_{tumor}_{normal}.dna.somatic.cnv.vcf.gz",
		"{output_dir}/{subject}_{tumor}_{normal}.dna.somatic.hard-filtered.vcf.gz",
		"{output_dir}/{subject}_{tumor}_{normal}.dna.somatic.sv.vcf.gz"
	shell:
		"dragen --enable-map-align false --bam-input {input.germline_bam} --tumor-bam-input {input.tumor_bam} -r {refgenome} --output-directory {output_dir} --output-file-prefix {wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic --enable-cnv true --intermediate-results-dir {temp_dir} --enable-variant-caller true --vc-enable-unequal-ntd-errors=true --vc-enable-trimer-context=true --enable-sv true --cnv-use-somatic-vc-baf true -f; mv {output_dir}/sv/results/variants/somaticSV.vcf.gz {output_dir}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.sv.vcf.gz; mv {output_dir}/sv/results/variants/somaticSV.vcf.gz.tbi {output_dir}/{wildcards.subject}_{wildcards.tumor}_{wildcards.normal}.dna.somatic.sv.vcf.gz.tbi"
