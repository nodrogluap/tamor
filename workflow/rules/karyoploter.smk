configfile: "config/config.yaml"

PERL_CMD = r"""perl -F\\t -ane 'print "$F[0]\t$F[1]\t",$F[1]+1,"\t$2\n" if $F[6] eq "PASS" and $F[$#F]=~/^(.+?:){3}([^:]+)/'"""

rule generate_karyoploter_cnv_plot:
        priority: 10
        input:
                somatic_snv_vcf=config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.hard-filtered.vcf.gz',
                somatic_cnv_vcf=config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.cnv.vcf.gz',
                # Evoke post germline and somatic cnv filtering (reduce artifacts in the plots).
                germline_filter_metrics_file=config["output_dir"]+'/{project}/{subject}/{subject}_{normal}.dna.germline.cnv_filter_metrics.csv',
                somatic_filter_metrics_file=config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.cnv_filter_metrics.csv'
        output:
                somatic_vaf=config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.hard-filtered.pass.vaf',
                jpeg=config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.cnv_plot.jpeg'
        conda:
                "../envs/plot_cnvs.yaml"
        resources:
                # Allow three hours for a plot
                runtime = 180,
                mem_mb = 50000
        shell:
                "{PERL_CMD} <(gzip -cd {input.somatic_snv_vcf}) > {output.somatic_vaf}; Rscript workflow/scripts/plot_cnvs.R {wildcards.project} {wildcards.subject} {wildcards.tumor} {wildcards.normal} {config[output_dir]}"
            
