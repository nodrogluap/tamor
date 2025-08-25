import tempfile

configfile: "config/config.yaml"

# CAPTIV-8 is a study by the Terry Fox Research Institute's Marathon of Hope Cancer Centres Network that tests the effectiveness of immune 
# checkpoint inhibitor (atezolizumab, a.k.a. Tecentriq) in individuals with advanced cancer. 
# See https://www.marathonofhopecancercentres.ca/our-research/project/marathon-of-hope-cancer-centres-network-study-for-ontario-canadian-atezolizumab-precision-targeting-for-immunotherapy-intervention-(mohccn-o-captiv-8)
# 
# The study uses a new scoring system called IBV, which analyzes the genetic information of tumors to choose participants. 
# Currently, a score of 5 or greater is the qualifying criterion for review by a tumor molecular board.
#
# This rule reformats Dragen whole genome + transcriptome outputs into the metrics needed for creating the score.
rule generate_captiv8_data_matrix:
        resources:
                runtime=180,
                mem_mb=40000 
        input:
		# For oncogenic viral integration check.
		oncogenic_viral=config["output_dir"]+"/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.oncogenic_viral_integration_metrics.tsv",
		# For M1 vs M2 macrophage scoring. Also Consensus Molecular Subtype assignment (1-4) in the case of colorectal samples. 
		rsem=config["output_dir"]+'/{project}/{subject}/rna/',
                # For tumor mutational burden check.
                tmb_metrics_tsv=config["output_dir"]+"/pcgr/{project}/{subject}_{tumor}_{normal}/{subject}.pcgr.grch38.tmb.tsv",
		# For tumour variant allele frequency and SWItch/Sucrose Non-Fermentable gene mutation analysis (a subfamily of ATP-dependent chromatin remodeling complexes).
                snvs=config["output_dir"]+"/pcgr/{project}/{subject}_{tumor}_{normal}/{subject}.pcgr.grch38.snv_indel_ann.tsv.gz",
		# For copy number loss of SWI/SNF genes check.
		cnvs=config["output_dir"]+"/pcgr/{project}/{subject}_{tumor}_{normal}/{subject}.pcgr.grch38.cna_gene.tsv.gz"
                # For tumor immune microenvironment profiling (in particular, the CD8+ T cell proportion is used). 
                # CAPTIV-8's RICO workflow uses CIBERSORT, PCGR uses QuantiSeq.
                #immunedeconv=config["output_dir"]+"/pcgr/{project}/{subject}_{tumor}_{normal}/{subject}"
        output:
		# This output file is formatted for use in the captiv8.R script
                matrix_file=config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.captiv8.data_matrix.tsv'
        conda:
		# Reusing an existing env that has minimap2 in it.
                "../envs/djerba.yaml"
        run:
		germline_library_info = identify_libraries(False, False, wildcards)
                germline_sample_libraries = library_info[1]
		tumor_library_info = identify_libraries(False, True, wildcards)
                tumor_sample_libraries = library_info[1]
                rna_library_info = identify_libraries(True, True, wildcards)
                rna_sample_libraries = rna_library_info[1]

		# TODO: Reformat PCGR quantiseq info to CIBERSORT format
		cd8 = 0.5 
		cibersort_file = tempfile.NamedTemporaryFile()
		with open(cibersort_file.name, 'wt') as cibersort:
			print("\t","T cell CD8+",str(cd8), file=cibersort)

		# TODO: CAPTIV-8 wants only protein coding tumor mutational burden.
		tmb=20

		# TODO: See of there are mutations in the SWISNF genes
		swisnf = "yes"

                # Get the oncotree code from the sample metadata, and infer if it's haematological or colorectal (required as part of scoring).
		# TODO
		colorectal = "yes"

		# TODO
		lymphoid = "yes"

                # TODO Scan the reads in the BAMs for signs of oncogenic viral integration.
                oncogenic_viral_integration = "yes"
                
		# Collate all the gathered info in the correct order for feeding to captiv8.R
		with open(output.matrix_file, "wt") as matrix:
			print("\t".join("Patient",wildcards.subject), file=matrix)
			print("\t".join("Libraries",",".join(germline_sample_libraries,tumor_sample_libraries,rna_sample_libraries)), file=matrix)
			print("\t".join("Cibersort file",cibersort_file.name, file=matrix)
			print("\t".join("RSEM file", input.rsem, file=matrix)
			print("\t".join("TMB", str(tmb), file=matrix)
			print("\t".join("SWISNF", swisnf, file=matrix)
			print("\t".join("Colorectal", colorectal, file=matrix)
			print("\t".join("Lymphoid", lymphoid, file=matrix)
			print("\t".join("Oncogenic viral integration", oncogenic_viral_integration, file=matrix)
                
rule generate_captiv8_score:
	input:
		matrix=config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.captiv8.data_matrix.tsv'
	output:
		# The captiv8.R script has a fixed output file name, so create subdirs for each tumor-normal pair within a sample
		# so the score details can co-exists.
		log=config["output_dir"]+'/{project}/{subject}/captiv8/{subject}_{tumor}_{normal}/captiv8_score.txt',
		details=config["output_dir"]+'/{project}/{subject}/captiv8/{subject}_{tumor}_{normal}/captiv8_output.txt',
		metrics=config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.captiv8_metrics.csv'
	conda:
		"../envs/djerba.yaml"
	shell:
		"Rscript workflow/submodules/djerba/src/lib/djerba/plugins/captiv8/Rscripts/captiv8.R -m {input.matrix} -o "+config["output_dir"]+
                '/{project}/{subject}/captive8/{subject}_{tumor}_{normal} -b workflow/submodules/djerba/src/lib/djerba/plugins/captiv8/gencode.v31.ensg_annotation_w_entrez.bed'+
                '> {output.log}'
