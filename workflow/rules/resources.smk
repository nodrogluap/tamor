configfile: "config/config.yaml"

rule download_resource_files:
	output:
		config["ref_fasta"],
		"resources/oncotree-taxonomy-2021-11-02.rdf",
		"resources/ontology_mappings.txt"
	run:
		"workflow/scripts/download_resources.py"
