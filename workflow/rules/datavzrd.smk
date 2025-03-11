rule dragen_metrics_view_with_datavzrd:
    input:
        config="resources/datavzrd/dragen_plots_and_tables.yaml",
        table="resources/metrics/dragen.tsv"
    output:
        report(
            directory("results/datavzrd/dragen_metrics"),
            htmlindex="index.html",
            caption="../report/dragen_metrics.rst",
            category="Metrics",
            labels={"table": "Dragen metrics"},
        )
    log:
        "logs/datavzrd.log"
    wrapper:
        # This points to a standardize wrapper repo known to Snakemake
        "v4.7.2/utils/datavzrd"
