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
