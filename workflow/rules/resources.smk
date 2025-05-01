configfile: "config/config.yaml"

rule download_resource_files:
	output:
		config["ref_genome"]+'/anchored_rna',
		config["ref_exon_annotations"],
		config["ref_fasta"],
		"resources/oncotree-taxonomy-2021-11-02.rdf",
		"resources/ontology_mappings.txt",
	shell:
		"workflow/scripts/download_resources.py"
