import csv
import glob

dna_paired_samples = {}
rna_paired_samples = {}

# The format of the paired DNA samples file is tab-delimited, with subject ID in the zeroth column, tumor sample name in first column, then matched normal in the second. 
# Additional columns are metadata: Boolean indicating if some tumor is expected in the normal sample, and a number representing the tumor site 
# as per PCGR guidelines (see README.md for details), then finally project name in the last column.
with open(config["dna_paired_samples_tsv"], 'r') as data_in:
        tsv_file = csv.reader(decomment(data_in), delimiter="\t")
        for line in tsv_file:
                if len(line[0]) > 35 or len(line[0]) < 6:
                        raise NameError("Subject names are limited to between 6 and 35 characters by PCGR reporting requirements, please revise "+line[0])
                if '_' in line[0]:
                        raise NameError("Subject name cannot contain any underscores ('_'), please revise "+line[0])
                if '_' in line[1]:
                        raise NameError("Tumor sample name cannot contain any underscores ('_'), please revise "+line[1])
                if '_' in line[2]:
                        raise NameError("Germline sample name cannot contain any underscores ('_'), please revise "+line[2])
                dna_paired_samples_key = "_".join((line[0],line[1],line[3]))
                dna_paired_samples[dna_paired_samples_key] = [line[0],line[5],line[6],line[7],line[2],line[4]]

with open(config["rna_paired_samples_tsv"], 'r') as data_in:
        tsv_file = csv.reader(decomment(data_in), delimiter="\t")
        for line in tsv_file:
                rna_paired_samples_key = "_".join((line[0],line[1]))
                rna_paired_samples[rna_paired_samples_key] = [line[0],line[1],line[3]]

def dna_paired_sample_keys():
                dna_paired_samples.keys()

def dna_paired_sample_values():
        dna_paired_samples.values()

def rna_paired_sample_keys():
        rna_paired_samples.keys()

def rna_paired_sample_values():
        rna_paired_samples.values()

# Return a tuple of two:
# Boolean as to whether this library should get the Unique Molecular Indices treatment, using info from the SampleSheet.csv
# The list of library names (that will be in the FASTQ list CSV) that correspond to the sample name identified by the wildcards 
def identify_libraries(is_rna, is_tumor, wildcards):
        # Figure out the library ID (which is in the FASTQ file names) from the sample name in the samplesheet
        # If XP loading was used, you can have multiple "libraries" called Li###-lane1, etc. for the same sample on a run as IEV requires unique sample IDs and normally we'd use the library ID
        sample_has_UMIs = False
        sample_libraries = []
        if is_tumor and hasattr(wildcards, "tumor"):
                target_sample = wildcards.tumor
        elif hasattr(wildcards, "normal"):
                target_sample = wildcards.normal
        elif hasattr(wildcards, "sample"):
                target_sample = wildcards.sample
        else:
                raise NameError("Missing any useable wildcard attribute (tumor, normal, or sample) in call to identify_libraries")
        for samplesheet in glob.iglob(config["samplesheet_dir"]+'/*.csv'):
                run_has_UMIs = False
                with open(samplesheet, 'r') as data_in:
                        sample_project_index = -1
                        sample_name_index = -1
                        sample_id_index = -1
                        csv_file = csv.reader(data_in)
                        for line in csv_file:
                                if "Sample_Name" in line and "Sample_ID" in line and "Sample_Project" in line: #it's the header for the sample list
                                        sample_name_index = line.index("Sample_Name")
                                        sample_id_index = line.index("Sample_ID")
                                        sample_project_index = line.index("Sample_Project")
				# New NextSeq2000 format with sample ID and sample name squished into one field by convention (hyphen separated)
                                elif "Sample ID*" in line and "Project" in line:
                                        sample_id_index = line.index("Sample ID*")
                                        sample_project_index = line.index("Project")
                                elif sample_name_index != -1 and line[sample_name_index] == target_sample or sample_id_index != -1 and line[sample_id_index].endswith("-"+target_sample):
                                        if is_rna:
                                                if "RNA" in line[sample_project_index] and (not "test" in line[sample_project_index]):
                                                        sample_libraries.append(line[sample_id_index])
                                                        sample_has_UMIs = run_has_UMIs
                                        else:
                                                if not "RNA" in line[sample_project_index]:
                                                        sample_libraries.append(line[sample_id_index])
                                                        sample_has_UMIs = run_has_UMIs
                                elif len(line) > 1 and line[0] == "OverrideCycles" and "U" in line[1]:
                                        run_has_UMIs = True
                        # Uncomment below if you want to be pedantic.
                        #if sample_name_index == -1:
                        #       raise NameError("Missing Sample_Name column in Illumina samplesheet "+samplesheet)
                        #if sample_id_index == -1:
                        #       raise NameError("Missing Sample_ID column in Illumina samplesheet "+samplesheet)
        return sample_has_UMIs, list(set(sample_libraries)) # dedup

def get_rna_projects():
        return [tuple[2] for tuple in rna_paired_samples.values()]

def get_rna_subjects():
        return [tuple[0] for tuple in rna_paired_samples.values()]

def get_rna_sample_keys():
        return list(rna_paired_samples.keys())

def get_dna_projects():
        return [tuple[3] for tuple in dna_paired_samples.values()]

def get_dna_subjects():
        return [tuple[0] for tuple in dna_paired_samples.values()]

def get_dna_sample_keys():
        return list(dna_paired_samples.keys())

def get_tumor_site(wildcards):
        return dna_paired_samples["_".join((wildcards.subject, wildcards.tumor, wildcards.normal))][2]

def get_normal_contains_some_tumor(wildcards):
        return dna_paired_samples["_".join((wildcards.subject, wildcards.tumor, wildcards.normal))][1]

def get_tumor_has_pcr_duplicates(wildcards):
        return dna_paired_samples["_".join((wildcards.subject, wildcards.tumor, wildcards.normal))][4]

# In theory since we record the PCR status of the normal on multiple lines of the config file, we will use the first instance encountered as correct.
# TODO: enforce PCR status to be the same across all DNA config lines a normal appears on?
def get_normal_has_pcr_duplicates(wildcards):
	for key, tuplevalues in dna_paired_samples.items():
		if key.startswith(wildcards.subject+"_") and key.endswith("_"+wildcards.normal):
        		return tuplevalues[5]
	return False


