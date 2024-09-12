#!/usr/bin/env python
import subprocess
import os
from urllib.parse import urlparse
from collections import OrderedDict

# Check which resources files don't already exist on the filesystem where expected by default.
# Since Python 3.7 dictionaries are ordered, so checking if ref_fasta exists is a cheap check against completion of this script. 

os.chdir("resources")
resource_dict = {
        "internal_pcgr": ["http://insilico.hpc.uio.no/pcgr/pcgr_ref_data.grch38.20240621.tgz", "data", ["tar","xzf"]],
        "internal_vep": ["https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz","homo_sapiens", ["tar","xzf"]],
        "ref_exon_annotations": ["https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz", "gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz"],
        "ref_genome": ["https://s3.amazonaws.com/webdata.illumina.com/downloads/software/dragen/resource-files/hg38-alt_masked.cnv.graph.hla.rna-9-r3.0-1.tar", "kmer_cnv.bin", ["tar"]],
        "ref_ora": ["https://s3.amazonaws.com/webdata.illumina.com/downloads/software/dragen/resource-files/misc/oradata_homo_sapiens.tar", "oradata_homo_sapiens", ["tar"]],
        "ref_fasta": ["http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa", "GRCh38_full_analysis_set_plus_decoy_hla.fa"]
}
for resource_name, spec in resource_dict.items():
    if(not os.path.exists(spec[1])):
        print("Fetching missing resource " + resource_name +  " from " + spec[0])
        # Adding -c so we can resume interrupted large downloads if possible
        subprocess.run(["wget", "-c", spec[0], "-‐output-document", spec[0]+".part"])
        subprocess.run(["mv", spec[0]+".part", spec[0]])
        # Post-process as required
        if length(spec) == 3:
            a = urlparse(spec[0])
            spec[3].append(os.path.basename(a.path))
            subprocess.run(spec[3])
