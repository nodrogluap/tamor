configfile: "config/config.yaml"

def get_file(project, subject, sample_name, is_rna, suffix):
        #print("Looking for " + project + " " + sample_name + " " + suffix)
        # Code that returns a list of bam files for based on *sample_name* 
        if is_rna:
                return config["output_dir"]+'/{project}/{subject}/rna/'+sample_name+suffix
        else:
                return config["output_dir"]+'/{project}/{subject}/'+sample_name+suffix

def get_normal_bam(wildcards):
        return get_file(wildcards.project, wildcards.subject, wildcards.normal, False, ".bam")

def get_tumor_bam(wildcards):
        return get_file(wildcards.project, wildcards.subject, wildcards.tumor, False, ".bam")

def get_tumor_rna_bam(wildcards):
        return get_file(wildcards.project, wildcards.subject, wildcards.tumor, True, ".bam")

