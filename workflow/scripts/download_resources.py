#!/usr/bin/env python
import subprocess
import os
from urllib.parse import urlparse

# Check which resources files don't already exist on the filesystem where expected by default.

try:
    os.chdir("resources")
except OSError:
    raise SystemExit("FATAL: No 'resources' directory exists under the current working directory ("+os.getcwd()+
            ", please be sure to run this command from the top level of your Tamor install")
print("INFO: Changed working directory to 'resources'")

resource_dict = {
        "internal_pcgr": ["http://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20250221.grch38.tgz", "data/grch38", "tar", "xzf"],
        "internal_vep": ["https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz","homo_sapiens/113_GRCh38", "tar", "xzf"],
        "ref_exon_annotations": ["https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz", "gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz"],
        "ref_genome": ["https://s3.amazonaws.com/webdata.illumina.com/downloads/software/dragen/resource-files/hg38-alt_masked.cnv.graph.hla.rna-9-r3.0-1.tar", "kmer_cnv.bin", "tar", "xf"],
        "ref_ora": ["https://s3.amazonaws.com/webdata.illumina.com/downloads/software/dragen/resource-files/misc/oradata_homo_sapiens.tar", "oradata_homo_sapiens", "tar", "xf"],
        "ref_fasta": ["http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa", "GRCh38_full_analysis_set_plus_decoy_hla.fa"]
}
for resource_name, spec in resource_dict.items():
    if(not os.path.exists(spec[1])):
        a = urlparse(spec[0])
        local_path = os.path.basename(a.path)
        print("INFO: Fetching missing resource " + resource_name +  " to resources/" + local_path + " from " + spec[0])

        # Adding -c so we can resume interrupted large downloads if possible
        result = subprocess.run(["wget", "--progress=bar:force:noscroll", "-c", spec[0], "--output-document", local_path+".part"])
        if result.returncode != 0:
            raise SystemExit("FATAL: Aborting reference data download, received non-zero return code ("+str(result.returncode)+") from command: " +
            "wget --progress=bar:force:noscroll -c " + spec[0] + "--output-document "+ local_path+".part ")

        result = subprocess.run(["mv", local_path+".part", local_path])
        if result.returncode != 0:
            raise SystemExit("FATAL: Received non-zero return code ("+str(result.returncode)+") from command: mv " + local_path+".part" + local_path)

        # Post-process as required, and assume we can delete the source file after this.
        if len(spec) >= 3:
            post_process_cmd = spec[2:]
            post_process_cmd.append(local_path)

            print("INFO: Running post-processing command:", post_process_cmd, sep=" ")
            result = subprocess.run(post_process_cmd)
            if result.returncode != 0:
                raise SystemExit("FATAL: Aborting reference data download, received non-zero return code ("+str(result.returncode)+") from command: " + post_process_cmd)

            print("INFO: Removing source file after post-processing: "+ local_path)
            result = subprocess.run(["rm", local_path])
            if result.returncode != 0:
                print("WARNING: Received non-zero return code ("+str(result.returncode)+") from source cleanup command: rm " + local_path)

# Illumina Annotation Engine (used for Tumor Mutational Burden reporting) has its own downloading system, in a standard place on Dragen servers
nirvana_downloader_path = "/opt/edico/share/nirvana/Downloader"
if(not os.path.exists(nirvana_downloader_path)):
    raise SystemExit("FATAL: Aborting reference data download, did not find Nirvana downloader script in expected location for a Dragen server: " + nirvana_downloader_path)
if(not os.path.exists("nirvana")):
    result = subprocess.run(["mkdir", "nirvana"])
    if result.returncode != 0:
        raise SystemExit("FATAL: Cannot create missing output path: " + os.getcwd() + "/nirvana")
print("INFO: Fetching Nirvana annotation databases with " + nirvana_downloader_path)
result = subprocess.run([nirvana_downloader_path, "--ga", "GRCh38", "--out", "nirvana"])
if result.returncode != 0:
    raise SystemExit("FATAL: Received non-zero return code (" + str(result.returncode) + ") from command: " + nirvana_downloader_path + " --ga GRCh38 --out nirvana")
