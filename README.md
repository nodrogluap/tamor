# tamor
Illumina Dragen cancer genome and transcriptome analysis automation using Snakemake, integrating Personal Cancer Genome Report generation on hg38.

# Installation

0. Nota bene: These insteructions assume that you already have a Dragen server with software version 3.10 or higher and a working ``hg38`` genome index.

1. Download this code base:

```bash
git clone https://github.com/nodrogluap/tamor
```

2. Install all the dependencies via conda or [mamba](https://mamba.readthedocs.io/en/latest/installation.html) (my preference because it's much, much faster):

```bash
mamba env create -f conda_tamor.yml
```

3. Due to a quirk in the conda dependencies spec, you will need to install the latest version of the Perl ``zlib`` library module manually, and the R hg38 sequence module:

```bash
cpanm Compress::Raw::Zlib
R -e 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")'
```

4. Download and index the cancer databases that CPSR and PCGR rely on for annotation of your discovered sequence variants:

```bash
BUNDLE=pcgr.databundle.hg38.20220203.tgz
wget http://insilico.hpc.uio.no/pcgr/$BUNDLE
tar zxvf $BUNDLE
```

# Configuration

If you do not have the tamor directory leading in your shell's ``PATH`` variable, you will need to prepend it so tamor's ersatz ``bcftools`` command (place this in your .bashrc if you don't want to do this manually each etime):

```bash
export PATH=/where/you/have/put/tamor:$PATH
```

Copy the ``config.yml.sample`` file to config.yml:

```bash
cp config.yml.sample config.yml
```

This is the file that you can customize for your site-specific settings. By default the config is set up to write result files under the current directory in ``data/output``, and is expecting the list of paired tumor-normal samples in a file called ``tumor_dna_paired_germline_dna_samples.tsv`` which has 5 columns:

```
subjectID<tab>cancer_sample_name<tab>germline_sample_name<tab>boolean_does_germline_data_contain_some_tumor<tab>PCGR_tissue_site_number
```

The cancer_sample_name and germline_sample_name must be exact names you used in your Illumina sequencing sample spreadsheet, as these sample sheets are the only metadata to which tamor has access. Place all the Illumina experiment sample sheets for your project into ``data/spreadsheets`` by default (see the ``samplesheets_dir`` setting in ``config.yml``).
Once the sample pairing file is ready, you can simply run Snakemake to generate the BAMs, VCFs, and CPSR/PCGR reports:
  
```bash
snakemake 
```
