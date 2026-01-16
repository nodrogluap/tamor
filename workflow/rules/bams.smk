configfile: "config/config.yaml"

def get_file(project, subject, sample_name, is_rna, suffix):
        #print("Looking for " + project + " " + sample_name + " " + suffix)
        # Code that returns a list of bam files for based on *sample_name* 
        if is_rna:
                return config["output_dir"]+'/'+project+'/'+subject+'/rna/'+subject+'_'+sample_name+suffix
        else:
                return config["output_dir"]+'/'+project+'/'+subject+'/'+subject+'_'+sample_name+suffix

# Use a function to find the alignment file, note that maybe the BAM has been converted to CRAM at some point to save disk space, outside Snakemake
def get_normal_bam(wildcards):
        cram = get_file(wildcards.project, wildcards.subject, wildcards.normal, False, ".dna.germline.cram")
        if exists(cram):
                return cram
        return get_file(wildcards.project, wildcards.subject, wildcards.normal, False, ".dna.germline.bam")

def get_tumor_bam(wildcards):
        cram = get_file(wildcards.project, wildcards.subject, wildcards.tumor+"_"+wildcards.normal, False, ".dna.somatic_tumor.cram")
        if exists(cram):
                return cram
        return get_file(wildcards.project, wildcards.subject, wildcards.tumor+"_"+wildcards.normal, False, ".dna.somatic_tumor.bam")

def get_tumor_rna_bam(wildcards):
        cram = get_file(wildcards.project, wildcards.subject, wildcards.tumor, True, ".rna.cram")
        if exists(cram):
                return cram
        return get_file(wildcards.project, wildcards.subject, wildcards.tumor, True, ".rna.bam")

def get_aligned_input_param(aligned_reads_file):
        # Change the dragen parameter name depending on whether the input is a BAM or CRAM
        if aligned_reads_file.endswith(".cram"):
                return f"--cram-input {aligned_reads_file}"
        return f"--bam-input {aligned_reads_file}"
