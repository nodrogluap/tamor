configfile: "config/config.yaml"

rule generate_minimal_maf:
        input:
                somatic_snv_vcf = config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.hard-filtered.vcf.gz'
        output:
                maf = config["output_dir"]+'/{project}/{subject}/{subject}_{tumor}_{normal}.dna.somatic.hard-filtered.maf.txt'
        run:
                # From the VCF generate the minimal MAF file that can be subsequently run through the Genome Nexus Annotation Pipeline 
                # (https://github.com/genome-nexus/genome-nexus-annotation-pipeline) to generate a 32 column TCGA MAF formatted file suitable
                # for upload to cBioPortal.
                # Minimal MAF (per https://docs.cbioportal.org/file-formats/#create-the-cbioportal-mutation-data-file-with-genome-nexus-with-a-minimal-maf-file): 
                # Chromosome (Required): A chromosome number, e.g., "7".
                # Start_Position (Required): Start position of event.
                # End_Position (Required): End position of event.
                # Reference_Allele (Required): The plus strand reference allele at this position.
                # Tumor_Seq_Allele2 (Required): Primary data genotype.
                # Tumor_Sample_Barcode (Required): This is the sample ID. Either a TCGA barcode (patient identifier will be extracted), or for non-TCGA data, a literal SAMPLE_ID as listed in the clinical data file.
                # In addition to the above columns, it is recommended to have the read counts to calculate variant allele frequencies:
                # t_alt_count (Optional, but recommended): Variant allele count (tumor).
                # t_ref_count (Optional, but recommended): Reference allele count (tumor).
                shell("gzip -cd {input.somatic_snv_vcf} | perl -F\\t -ane '$F[0] =~ s/^chr//; print join(\"\\t\", $F[0], $F[1], $F[1]+length($F[3])-1, $F[3], $F[4], $2, $1, \"{wildcards.tumor}\"), \"\\n\" if not /^#/ and $F[$#F] =~ /^[^:]+:[^:]+:(\\d+),(\\d+)/' > {output.maf}")

