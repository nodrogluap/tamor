configfile: "config/config.yaml"

rule get_metrics:
    output:
        "resources/metrics/dragen.tsv"
    conda:
        "../envs/metrics.yaml"
    log:
        notebook="logs/get_metrics.ipynb"
    notebook:
        "../notebooks/get_metrics.py.ipynb"

rule get_rna_stats:
    output:
        statsfile = "resources/metrics/rna_stats.tsv"
    log:
        "logs/get_rna_stats.csv"
    shell:
        "workflow/scripts/get_rna_stats.sh "+config["output_dir"]+" {output.statsfile}"
