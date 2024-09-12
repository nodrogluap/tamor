configfile: "config/config.yaml"

def get_file(project, subject, sample_name, is_rna, suffix):
        #print("Looking for " + project + " " + sample_name + " " + suffix)
        # Code that returns a list of bam files for based on *sample_name* 
        if is_rna:
                return config["output_dir"]+'/{project}/{subject}/rna/{subject}_'+sample_name+suffix
        else:
                return config["output_dir"]+'/{project}/{subject}/{subject}_'+sample_name+suffix

def get_normal_bam(wildcards):
        return get_file(wildcards.project, wildcards.subject, wildcards.normal, False, ".dna.germline.bam")

def get_tumor_bam(wildcards):
        return get_file(wildcards.project, wildcards.subject, wildcards.tumor, False, "_"+wildcards.normal+".dna.somatic.bam")

def get_tumor_rna_bam(wildcards):
        return get_file(wildcards.project, wildcards.subject, wildcards.tumor, True, ".rna.bam")

