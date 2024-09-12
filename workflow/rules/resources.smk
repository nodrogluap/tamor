configfile: "config/config.yaml"

rule download_resource_files:
	output:
		config["ref_fasta"]
	run:
		"workflow/scripts/download_resources.py"
