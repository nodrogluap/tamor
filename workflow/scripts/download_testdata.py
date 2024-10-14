#!/usr/bin/env python
import subprocess
import yaml

# To be executed from Tamor root dir
with open("config/config.yaml", "r") as f:
    config = yaml.safe_load(f)

subprocess.run(["workflow/scripts/download_testdata.sh", config["ref_fasta"]])
