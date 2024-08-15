#!/usr/bin/env python
import subprocess
import os

# To be executed from Tamor root dir
with open("config/config.yaml", "r") as f:
    config = yaml.safe_load(f)

# Check which resources files specified in the config don't already exist on the filesystem where expected.

os.chdir("resources")
resource_dict = {
        "internal_pcgr": "http://insilico.hpc.uio.no/pcgr/pcgr_ref_data.grch38.20240621.tgz", 
        "internal_vep": "https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz",
        "ref_fasta": "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa", 
        "ref_exon_annotations": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz", 
        "ref_genome": "https://s3.amazonaws.com/webdata.illumina.com/downloads/software/dragen/resource-files/hg38-alt_masked.cnv.graph.hla.rna-9-r3.0-1.tar",
        "ref_ora": "https://s3.amazonaws.com/webdata.illumina.com/downloads/software/dragen/resource-files/misc/oradata_homo_sapiens.tar"
}
for resource_name in resource_names.keys():
    if(not os.path.isfile(config[resource_name])):
        # Adding -c so we can resume interrupted large downloads if possible
        os.subprocess(["wget", "-c", config[resource_name]])
        
