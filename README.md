# tamor
Illumina Dragen cancer genome and transcriptome analysis automation using Snakemake, integrating Personal Cancer Genome Report generation.

# Getting Started

Download this code base:

```bash
git clone https://github.com/nodrogluap/tamor
```

Install all the dependencies via conda or mamba (except the Dragen command itself, which as a commercial piece of software we're assuming is pre-installed on your server):

```bash
mamba env create -f conda_tamor.yml
```

Copy the ``config.yml.sample`` file to config.yml:

```bash
cp config.yml.sample config.yml
```

By default the config is set up to write result files under the current directory in ``data/output``, and is expecting the list of paired tumor-normal samples in a file called ``tumor_dna_paired_germline_dna_samples.tsv`` which has 5 columns:

```
subjectID<tab>cancer_sample_name<tab>germline_sample_name<tab>boolean_does_germline_data_contain_some_tumor<tab>PCGR_tissue_site_number
```

The cancer_sample_name and germline_sample_name must be exact names you used in your Illumina sequencing sample spreadsheet, as these sample sheets are the only metadata to which tamor has access. Place all the Illumina experiment sample sheets for your project into ``data/spreadsheets`` by default (see the ``samplesheets_dir`` setting in ``config.yml``).
Once the sample pairing file is ready, you can simply run Snakemake to generate the BAMs, VCFs, and CPSR/PCGR reports:
  
```bash
snakemake 
```
