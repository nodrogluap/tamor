# tamor

Rapid automated [Personal Cancer Genome Report](https://sigven.github.io/pcgr/) generation using [Illumina Dragen](https://www.illumina.com/products/by-type/informatics-products/dragen-secondary-analysis.html) + [Snakemake](https://snakemake.github.io/).

# tl;dr

Data for large scale tumor analysis projects can be spread over multiple DNA sequencing instrument runs, ``tamor`` simplifies the process of analyzing them.

A tab-delimited file is provided by the user to associate tumor and germline sequencing sample IDs with a study subject ID, along with a tissue-of-origin for the tumor. PCGR somatic variant reports (including germline susceptibility sequence variants) are generated using 1) this tab-delimited file, 2) the Illumina sequencer output (BCL or FASTQ), and 3) the Illumina Experiment Manager samplesheets CSV for the sequencing runs.

Tumor RNA analysis is in development.

# Installation

0. Nota bene: These instructions assume that you already have a Dragen server with software version 3.10 or higher and a working ``hg38`` genome index.

1. Download this code base:

```bash
git clone https://github.com/nodrogluap/tamor
```

2. Install all the dependencies via conda or [mamba](https://mamba.readthedocs.io/en/latest/installation.html) (my preference because it's much, much faster):

```bash
mamba env create -f conda_tamor.yml
```

3. Due to quirks in the conda dependencies spec, you will need to install the latest version of the Perl ``zlib`` library module manually, and the R hg38 genome sequence module:

```bash
cpanm Compress::Raw::Zlib
R -e 'BiocManager::install("rtracklayer", force=TRUE);BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")'
```

4. Download (~22GB) the cancer databases that CPSR and PCGR rely on for annotation of your discovered sequence variants:

```bash
BUNDLE=pcgr.databundle.hg38.20220203.tgz
wget http://insilico.hpc.uio.no/pcgr/$BUNDLE
tar zxvf $BUNDLE
```

# Configuration

If you do not have the tamor directory leading in your shell's ``PATH`` variable, you will need to prepend it so tamor's ersatz ``bcftools`` command is used (place this in your .bashrc if you don't want to do this manually each time):

```bash
export PATH=/where/you/have/put/tamor:$PATH
```

Copy the ``config.yml.sample`` file to config.yml:

```bash
cp config.yml.sample config.yml
```

This is the file that you can customize for your site-specific settings. By default the config is set up to write result files under the current directory in ``output``, and is expecting the input list of paired tumor-normal samples in a file called ``tumor_dna_paired_germline_dna_samples.tsv`` which has 5 columns to be specified:

```
subjectID<tab>tumorSampleName<tab>germlineSampleName<tab>TrueOrFalse_germline_contains_some_tumor<tab>PCGRTissueSiteNumber
```

The ``subjectID``, ``tumorSampleName`` and ``germlineSampleName`` must:

- *CONTAIN NO UNDERSCORES*
- The ``subjectID`` must be between 6 and 35 characters (due to a PCGR naming limitation)
- ``tumorSampleName`` and ``germlineSampleName`` must be the exact ``Sample_Name`` values you used in your Illumina sequencing sample spreadsheets

These sample sheets are the only metadata to which tamor has access. Place all the Illumina experiment sample sheets for your project into ``data/spreadsheets`` by default (see the ``samplesheets_dir`` setting in ``config.yml``). They must be called ``runID.csv`` where runID is typically the Illumina folder name in the format ``YYMMDD_machineID_SideFlowCellID``.

Tamor can start with either BCL files or FASTQ. If you are starting with BCLs, the full Illumina experiment output folders (which contain the requisite ``Data/Intensities/Basecalls`` subfolder) are expected by in ``data/bcls/runID`` (see ``bcl_dir`` setting in``config.yaml``). Tamor will perform bcl to fastq conversion, with the FASTQ output into ``data/analysis/primary/sequencer/runID`` (see ``analysis_dir`` setting in ``config.yaml``, and the default ``sequencer`` is ``novaseq6000``). 

If instead you are providing the FASTQs directly as input to tamor, they must also be in the ``data/analysis/primary/sequencerName/runID`` directory, with a corresponding Illumina Experiment Manager samplesheet ``data/spreadsheets/runID.csv``. *Why?* This is required because tamor reads the sample sheet to find the correspondence between Sample_Name and Sample ID for each sequencing library (multiple Sample IDs can correspond to the same tumor sample e.g. if more than one prep was done, XP loading was used, or there are multiple barcodes for one sample for color balancing on small runs). The samplesheet is also used to determine if Unique Molecular Indices were used to generate the sequencing libraries, which requires different handling in Dragen during genotyping downstream.

*If you provide FASTQ files directly, they must be timestamped later than the corresponding Illumina Experiment Manager spreadsheet, otherwise Snakemake will assume you've consequentially changed the spreadsheet and try to automatically regenerated all FASTQs for that run -- from potentially non-existent BCLs*.

The list of tissue site numbers for the version of PCGR included here is:

```
                        0 = Any
                        1 = Adrenal Gland
                        2 = Ampulla of Vater
                        3 = Biliary Tract
                        4 = Bladder/Urinary Tract
                        5 = Bone
                        6 = Breast
                        7 = Cervix
                        8 = CNS/Brain
                        9 = Colon/Rectum
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

# Running a paired tumor-normal analysis

Any time you want to use tamor, you must be sure to have the conda/mamba environment loaded:

```bash
mamba activate pcgrr
```

Once the sample pairing file mentioned earlier is ready, you can simply run Snakemake to generate the FASTQs (optiuonally), BAMs, VCFs, and CPSR/PCGR reports:
  
```bash
snakemake --cores=2
```
The default outputs are in a directory called ``data/output/sampleID_tumorSampleName_germlineSampleName``. The most relevant document may be the self-contained Web page ``sampleID.pcgr_acmg.grch38.flexdb.html``.

![Screenshot of a sample Personal Cancer Genome Report, FlexDB version](docs/pcgr_screenshot.png)

# Acknowledgements

This project is being developed in support of the [Terry Fox Research Institute](https://www.tfri.ca/)'s [Marathon of Hope Cancer Care Network](https://www.marathonofhopecancercentres.ca/) activities within the [Prairie Cancer Research Consortium](https://www.marathonofhopecancercentres.ca/our-network/consortium/prairies-cancer-research-consortium).
