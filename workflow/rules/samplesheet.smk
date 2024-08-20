import csv
from pathlib import Path

configfile: "config/config.yaml"

def get_samplesheet(wildcards):
        return config["samplesheet_dir"]+'/'+wildcards.run+'.csv'

def load_samplesheet_info():
        if samplesheet_cache:
                return samplesheet_cache
        for path in Path(config["samplesheet_dir"]).iterdir():
                if path.is_file() and Path(path).suffix == '.csv':
                        #print("Reading file " + str(path))
                        samplesheet_lines = []
                        current_file = open(path, "r")
                        csv_file = csv.reader(current_file)
                        sample_name_index = -1
                        sample_id_index = -1
                        sample_project_index = -1
                        try:
                                for line in csv_file:
                                        # It's the header for the sample list.
                                        if 'Sample_Name' in line and 'Sample_ID' in line and 'Sample_Project' in line:
                                                sample_name_index = line.index('Sample_Name')
                                                sample_id_index = line.index('Sample_ID')
                                                sample_project_index = line.index('Sample_Project')
                                        # It's a line for which we might have some use.
                                        if sample_name_index != -1 and sample_id_index != -1 and sample_project_index != -1:
                                                if str(path) not in samplesheet_cache:
                                                        samplesheet_cache[str(path)] = samplesheet_lines
                                        samplesheet_lines.append(line)
                        # Usefully handle garbage text in the samplesheets rather than having a generic error/traceback.
                        except UnicodeDecodeError as ude:
                                lineGuess = ude.object[:ude.start].count(b'\n') + 1
                                sys.exit("Found a bad char (non-UTF-8) in file " + str(path) + " at byte " + str(ude.start) + " around line " + str(lineGuess))

                        if sample_name_index == -1:
                                print("Missing Sample_Name column in Illumina samplesheet "+str(path)+", skipping")
                        if sample_id_index == -1:
                                print("Missing Sample_ID column in Illumina samplesheet "+str(path)+", skipping")
                        if sample_project_index == -1:
                                print("Missing Sample_Project column in Illumina samplesheet "+str(path)+", skipping")
