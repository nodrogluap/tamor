# Tamor

Rapid automated [Personal Cancer Genome Report](https://sigven.github.io/pcgr/) (PCGR) generation using 
[Illumina Dragen](https://www.illumina.com/products/by-type/informatics-products/dragen-secondary-analysis.html) + 
[Snakemake](https://snakemake.github.io/), handling both genomic and transcriptomic input data.

# tl;dr

Data for large scale tumor analysis projects can be spread over multiple DNA sequencing instrument runs, ``Tamor`` simplifies the process of analyzing them.

Tab-delimited files are configured by the user to associate tumor (DNA and/or RNA) and germline sequencing sample IDs with a study subject ID, along with a tissue-of-origin for the tumor. 
Somatic variants (including small nucleotide variants, structural variants and copy number variants) as well as gene expression reports are generated using 
1) these tab-delimited files
2) the Illumina sequencer output (BCL or FASTQ), and
3) the Illumina Experiment Manager samplesheets (CSV) for the sequencing runs

# Table of Contents
* [Prerequisites](#prerequisites)
* [Installation](#installation)
* [Testing](#testing)
* [Running Tamor](#running-tamor)
* [Configuration](#configuration)
* [Acknowledgements](#acknowledgements)

# Prerequisites

This workflow is intended for people with an Illumina Dragen hardware-accelerated (FPGA) system for high-throughput genomics analysis.  
If you don't have this hardware, this probably isn't for you.

Assuming you have a fresh Dragen server, you will need to download a 
(Dragen-formatted human reference genome)[https://support.illumina.com/sequencing/sequencing_software/dragen-bio-it-platform/product_files.html], 
preferably hg38. You can put this anywhere on the filesystem; by default Tamor will expect it in 
a subfolder of ``/usr/local/illumina/genomes`` (see [Configuration](#configuration)).

You will also need to set up a [slurm](https://slurm.schedmd.com/quickstart_admin.html#quick_start) queue so that jobs running on the Dragen FPGA don't collide with each other.
This is recommended by Illumina support, but not part of the Dragen documentation.

# Installation

0. [Install the mamba package manager](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) if you don't already have it on your system.

1. Create a mamba or conda environment for the latest Snakemake (8.something) and git:

```bash 
mamba create -c conda-forge -c bioconda -n snakemake snakemake git
mamba activate snakemake
```
2. Download the Tamor code:

```bash
git clone https://github.com/nodrogluap/tamor
cd tamor
```

3. Download (~22GB) the cancer annotation databases that [CPSR](https://github.com/sigven/cpsr) and PCGR rely on for annotating your discovered sequence variants:

```bash
workflow/scripts/download_resources.py
```

# Testing

Tamor follows the Snakemake [Distribution and Reproducibility](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) guidelines, so files are located in standardized locations.
The default config files are pre-configured for running a single test case from the [NCBI Short Read Archive](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA433607).  This case of apparent Chronic
Lymphocytic Leukemia (CLL) has both tumor DNA (30x coverage) and RNA data (27M) available, both 2x150bp paired-end Illumina.

If you would like to run the test case before reconfiguring Tamor to use your own data, you will need to download and format the CLL SRA records. 
This requires some additional specialty software not otherwise required by Tamor, so you will need to install a test mamba environment first.

```bash
mamba create env -f workflow/envs/test.yaml
mamba activate test
workflow/scripts/download_testdata.py
mamba deactivate test
```

This can take a few hours depending on your Internet connection speed, and requires at least 40GB of RAM to generate matched-pseudonormal FASTQ files from the cancer sample FASTQ files.

# Running Tamor

Once either the test data or your own (see [Configuration section below](#configuration)) is ready, you can run Snakemake to generate the 
[BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map) files, 
[VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) files, and 
[CPSR/PCGR reports](https://sigven.github.io/pcgr/index.html).

In a multi-user system, it is imperative to use a queuing system such as slurm to submit only one job at a time to Dragen v4.x. 
Once slurm is installed and configured on your Dragen system, Snakemake support for slurm is enabled by invoking like so:
  
```bash
snakemake --cluster sbatch --use-conda --cores=1
```

The default outputs are in a directory called ``results/pcgr/projectID/subjectID_tumorSampleID_germlineSampleID``. 
The most relevant document may be the self-contained Web page ``subjectID.pcgr.grch38.html``.

![Screenshot of a sample Personal Cancer Genome Report, Quarto version](docs/pcgr_screenshot.png)

# Configuration

## config/config.yaml
This is the file that you can customize for your site-specific settings. 
By default the config is set up to read input files from the ``resources`` folder, and write result files under the ``results`` folder. 
Tamor by default is expecting the input lists of paired tumor-normal samples in files called ``dna_samples.tsv`` and ``rna_samples.tsv``.


## config/dna_samples.tsv
Has 6 columns to be specified:

```
subjectID<tab>tumorSampleID<tab>germlineSampleID<tab>TrueOrFalse_germline_contains_some_tumor<tab>PCGRTissueSiteNumber<tab>ProjectID
```

The ``subjectID``, ``tumorSampleID`` and ``germlineSampleID`` must:

- *CONTAIN NO UNDERSCORES*
- The ``subjectID`` must be between 6 and 35 characters (due to a PCGR naming limitation)
- ``tumorSampleID`` and ``germlineSampleID`` must be the exact ``Sample_Name`` values you used in your Illumina sequencing sample spreadsheets (see samplesheet section below for details).

The *sixth* column is a unique project ID to which the subject belongs. For example if you have two cohorts of lung and breast cancer, 
assigning individuals to two projects would be logical. 
All project output files go into their own output folders, even if they were sequenced together on the same Illumina sequencing runs.

The *fourth* column of the paired input sample TSV file is usually ``False``, unless your germline sample is from a leukemia or perhaps 
a poor quality histology section from a tumor, in which case use ``True``. This instructs Dragen to consider low frequency variants 
in the germline sample to still show up as somatic variants in the tumor analysis output (see default of 0.05 under ``tumor_in_normal_tolerance_proportion`` in ``config.yaml``)

For the *fifth* column, the list of tissue site numbers for the version of PCGR included here is:

```
                        0  = Any
                        1  = Adrenal Gland
                        2  = Ampulla of Vater
                        3  = Biliary Tract
                        4  = Bladder/Urinary Tract
                        5  = Bone
                        6  = Breast
                        7  = Cervix
                        8  = CNS/Brain
                        9  = Colon/Rectum
                        10 = Esophagus/Stomach
                        11 = Eye
                        12 = Head and Neck
                        13 = Kidney
                        14 = Liver
                        15 = Lung
                        16 = Lymphoid
                        17 = Myeloid
                        18 = Ovary/Fallopian Tube
                        19 = Pancreas
                        20 = Peripheral Nervous System
                        21 = Peritoneum
                        22 = Pleura
                        23 = Prostate
                        24 = Skin
                        25 = Soft Tissue
                        26 = Testis
                        27 = Thymus
                        28 = Thyroid
                        29 = Uterus
                        30 = Vulva/Vagina
```

## config/rna_samples.tsv
Has 5 columns to be specified:

```
subjectID<tab>tumorRNASampleID<tab>matchedTumorDNASampleID<tab>ProjectID<tab>ImmuneDeconvCancerType
```

Where [ImmuneDeconv](https://omnideconv.org/immunedeconv/articles/immunedeconv.html)CancerType is 
one of the following (pick what seems closest if no exact match is available):

```
                        acc  = Adrenocortical carcinoma
                        blca = Bladder Urothelial Carcinoma
                        lgg  = Brain Lower Grade Glioma
                        brca = Breast invasive carcinoma
                        cesc = Cervical squamous cell carcinoma and endocervical adenocarcinoma
                        chol = Cholangiocarcinoma
                        coad = Colon adenocarcinoma
                        esca = Esophageal carcinoma
                        gbm  = Glioblastoma multiforme
                        hnsc = Head and Neck squamous cell carcinoma
                        kich = Kidney Chromophobe
                        kirc = Kidney renal clear cell carcinoma
                        kirp = Kidney renal papillary cell carcinoma
                        lihc = Liver hepatocellular carcinoma
                        luad = Lung adenocarcinoma
                        lusc = Lung squamous cell carcinoma
                        dlbc = Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
                        meso = Mesothelioma
                        ov   = Ovarian serous cystadenocarcinoma
                        paad = Pancreatic adenocarcinoma
                        pcpg = Pheochromocytoma and Paraganglioma
                        prad = Prostate adenocarcinoma
                        read = Rectum adenocarcinoma
                        sarc = Sarcoma
                        skcm = Skin Cutaneous Melanoma
                        stad = Stomach adenocarcinoma
                        tgct = Testicular Germ Cell Tumors
                        thym = Thymoma
                        thca = Thyroid carcinoma
                        ucec = Uterine Corpus Endometrial Carcinoma
                        uvm  = Uveal Melanoma
                        ucs  = Uterine Carcinosarcoma
```

## resources/samplesheets
These sample sheets are the only other metadata to which Tamor has access. Place all the 
[Illumina experiment sample sheets](https://support.illumina.com/downloads/sample-sheet-v2-template.html) for your project 
into ``resources/spreadsheets`` by default (see the ``samplesheets_dir`` setting in ``config/config.yaml``). They must be 
called ``runID.csv``, where runID is typically the Illumina folder name in the format ``YYMMDD_machineID_SideFlowCellID``.

Tamor can start with either BCL files or FASTQ. If you are starting with BCLs, the full Illumina experiment output folders (which contain the 
requisite ``Data/Intensities/Basecalls`` subfolder) are expected by in ``resources/bcls/runID`` (see ``bcl_dir`` setting in``config.yaml``). Tamor will perform BCL 
to FASTQ conversion, with the FASTQ output into ``results/analysis/primary/sequencer/runID`` (see ``analysis_dir`` setting in ``config.yaml``, and the 
default ``sequencer`` is ``HiSeq`` per the test data mentioned earlier). 

If instead you are providing the FASTQs directly as input to Tamor, they must also be in the ``resources/analysis/primary/sequencerName/runID`` directory, 
with a corresponding Illumina Experiment Manager samplesheet ``resources/spreadsheets/runID.csv``. *Why?* This is required because Tamor reads the sample 
sheet to find the correspondence between Sample_Name and Sample ID for each sequencing library, also analysis for DNA samples differs from that for RNA 
samples, so the sample sheet must also contain a ``Sample_Project`` column. Sample projects with names that contain "RNA" in them will be processed as such, 
all others are assumed to be DNA. The ``Sample_Project`` is not used for any other purpose than distinguishing RNA and DNA, and does not need to be the 
same as the ProjectIDs listed in the ``config`` folder files.

The samplesheet is also used to determine if Unique Molecular Indices were used to generate the sequencing libraries, which requires different handling in Dragen during genotyping downstream.

*If you provide FASTQ files directly, they must be timestamped later than the corresponding Illumina Experiment Manager spreadsheet, 
otherwise Snakemake will assume you've consequentially changed the spreadsheet and try to automatically regenerated all FASTQs 
for that run -- from potentially non-existent BCLs*.

# Acknowledgements

This project is being developed in support of the [Terry Fox Research Institute](https://www.tfri.ca/)'s 
[Marathon of Hope Cancer Care Network](https://www.marathonofhopecancercentres.ca/) activities within the 
[Prairie Cancer Research Consortium](https://www.marathonofhopecancercentres.ca/our-network/consortium/prairies-cancer-research-consortium).
