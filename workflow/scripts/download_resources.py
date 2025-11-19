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
        "internal_pcgr": ["http://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20250314.grch38.tgz", "data/grch38", "tar", "xzf"],
        "internal_vep": ["https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz","homo_sapiens/113_GRCh38", "tar", "xzf"],
        "ref_exon_annotations": ["https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz", "gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz"],
        "ref_genome": ["https://s3.amazonaws.com/webdata.illumina.com/downloads/software/dragen/resource-files/hg38-alt_masked.cnv.graph.hla.rna-9-r3.0-1.tar", "kmer_cnv.bin", "tar", "xf"],
        "ref_ora": ["https://s3.amazonaws.com/webdata.illumina.com/downloads/software/dragen/resource-files/misc/oradata_homo_sapiens.tar", "oradata_homo_sapiens", "tar", "xf"],
        "ref_fasta": ["http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa", "GRCh38_full_analysis_set_plus_decoy_hla.fa"],
        "cnv_germline_b_allele_sites": ["https://webdata.illumina.com/downloads/software/dragen/resource-files/misc/hg38_1000G_phase1.snps.high_confidence.vcf.gz", "hg38_1000G_phase1.snps.high_confidence.vcf.gz"],
        "alu_vc_exclusions_bed": ["https://webdata.illumina.com/downloads/software/dragen/resource-files/bed-file-collection-1.0.0.tar", "bed-file-collection-1.0.0", "tar", "xvf"],
        "snv_systematic_noise": ["https://webdata.illumina.com/downloads/software/dragen/resource-files/systematic-noise-baseline-collection-1.1.0.tar", "systematic-noise-baseline-collection-1.1.0", "tar", "vxf"],
        "sv_systematic_noise": ["https://webdata.illumina.com/downloads/software/dragen/resource-files/sv-systematic-noise-baseline-collection-2.0.1.tar", "sv-systematic-noise-baseline-collection-2.0.1", "tar", "xvf"],
        "oncotree_mappings": ["https://raw.githubusercontent.com/cBioPortal/oncotree/refs/heads/master/scripts/ontology_to_ontology_mapping_tool/ontology_mappings.txt", "ontology_mappings.txt"],
        "oncotree_hierarchy": ["https://raw.githubusercontent.com/cBioPortal/oncotree/refs/heads/master/resources/rdf/oncotree-taxonomy-2021-11-02.rdf", "oncotree-taxonomy-2021-11-02.rdf"],
        "graphkb_canonical_transcripts": ["https://raw.githubusercontent.com/bcgsc/pori/feature/colab-notebooks/demo/mart_export.protein_coding.canonical.txt", "mart_export.protein_coding.canonical.txt"]
}

# Illumina Annotation Engine (used for Tumor Mutational Burden reporting) has its own downloading system, in a standard place on Dragen 4.2 servers, variable on v4.4+
is_dragen_v42 = True
nirvana_downloader_path = "/opt/edico/share/nirvana/Downloader"
install_path = ""
if(not os.path.exists(nirvana_downloader_path)):
    # In Dragen v4.4 it's changed to be from the PATH so that multiple Dragen versions can be supported on one server.
    result = subprocess.run(["/usr/bin/which", "dragen_info"], capture_output=True, text=True, check=True)
    if result.stderr:
        raise SystemExit("FATAL: Aborting reference data download, did not find Dragen install in path (assuming Dragen v4.4+)")   
    install_path = result.stdout.removesuffix("/bin/dragen_info\n")
    nirvana_downloader_path = install_path+"/share/nirvana/DataManager"
    if(not os.path.exists(nirvana_downloader_path)):
        raise SystemExit("FATAL: Aborting reference data download, did not find Nirvana downloader script in expected location for a Dragen server: " + nirvana_downloader_path)
    is_dragen_v42 = False
if(not os.path.exists("nirvana")):
    result = subprocess.run(["mkdir", "nirvana"])
    if result.returncode != 0:
        raise SystemExit("FATAL: Cannot create missing output path: " + os.getcwd() + "/nirvana")
print("INFO: Fetching Nirvana annotation databases with " + nirvana_downloader_path)
if is_dragen_v42:
    result = subprocess.run([nirvana_downloader_path, "--ga", "GRCh38", "--out", "nirvana"])
    if result.returncode != 0:
        raise SystemExit("FATAL: Received non-zero return code (" + str(result.returncode) + ") from command: " + nirvana_downloader_path + " --ga GRCh38 --out nirvana")
else:
    credentials_file = "credentials.json" # remember, we're already in the resources dir
    if(not os.path.exists(credentials_file)):
        raise SystemExit("FATAL: For Dragen v4.4+ an Illumina API key is required to download annotation resources, please generate (per https://help.dragen.illumina.com/product-guide/dragen-v4.4/nirvana#premium-sources) and place an Illumina API key JSON in resources/"+credentials_file)
    config_file = install_path+"/resources/annotation/all_annotations_GRCh38.json"
    if(not os.path.exists(config_file)):
        raise SystemExit("FATAL: For Dragen v4.4+ the expected default annotation config JSON is not where we expect it ("+config_file+")")
    subprocess.run([install_path+"/share/nirvana/DataManager", "download", "-r", "GRCh38", "--credentials-file",
                   credentials_file, "--dir", "nirvana", "--versions-config", config_file])

# Update the reference genome index as required for v4.4+
if(not is_dragen_v42):
        resource_dict["ref_genome"] = ["https://s3.us-east-1.amazonaws.com/webdata.illumina.com/downloads/software/dragen/references/genome-files/hg38-alt_masked.cnv.graph.hla.methyl_cg.rna-11-r5.0-1.tar.gz", "kmer_cnv.bin", "tar", "zxf"]
        resource_dict["snv_systematic_noise"] = ["https://webdata.illumina.com/downloads/software/dragen/resource-files/misc/systematic-noise-baseline-collection-2.0.0.tar", "systematic-noise-baseline-collection-2.0.0", "tar", "vxf"]
        resource_dict["sv_systematic_noise"]: ["https://webdata.illumina.com/downloads/software/dragen/resource-files/4.4/sv-systematic-noise-baseline-collection-v3.1.0-1.tar.gz", "WGS_hg38_v3.1.0_systematic_noise.sv.bedpe.gz", "tar", "zxvf"]

for resource_name, spec in resource_dict.items():
    if(not os.path.exists(spec[1])):
        a = urlparse(spec[0])
        local_path = os.path.basename(a.path)
        print("INFO: Fetching missing resource " + resource_name +  " to resources/" + local_path + " from " + spec[0])

        # Adding -c so we can resume interrupted large downloads if possible
        result = subprocess.run(["wget", "--progress=bar:force:noscroll", "-c", spec[0], "--output-document", local_path+".part"])
        if result.returncode != 0:
            raise SystemExit("FATAL: Aborting reference data download, received non-zero return code ("+str(result.returncode)+") from command: " +
            "wget --progress=bar:force:noscroll -c " + spec[0] + " --output-document "+ local_path+".part ")

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

