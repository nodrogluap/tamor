# Using Tamor without any local filesystem changes

In some instances, it may be advantageous to run Tamor without making any changes to the user's home directory (e.g. due to file quotas or lack of backups). 
The conda/mamba system used to manage software dependencies in Tamor will normally put all these files under the home directory, but this can be avoided and dependencies can be written instead
to any remote filesystem mount point with a few extra Linux shell setup steps.

## Initialization
To create a standalone Tamor on a mounted filesystem at /path/to/your/mount/point without a local mamba install, on the target machine that has /path/to/your/mount/point remote-mounted (e.g. via NFS or CIFS) do the following:

### 1. Install the standalone micromamba executable

```bash
cd /path/to/your/mount/point
export MAMBA_ROOT_PREFIX=`pwd`/mambaforge
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

### 2. Activate micromamba in the current shell
```bash
eval "$(micromamba shell hook --shell bash)"
micromamba activate
```

### 3. Install core Snakemake and git dependencies in a micromamba env called "tamor"
This should install ~187 packages from ~155MB of downloads on the mountpoint.

```bash
micromamba create -c conda-forge -c bioconda -n tamor mamba snakemake git git-lfs wget conda=24.7.1
```

### 4. Download the Tamor code
This will use git from the micromamba env called "tamor" that we just installed

```bash
micromamba activate tamor
git clone --recurse-submodules https://github.com/nodrogluap/tamor
```

## Running Tamor
Anytime you would like to run Tamor, you must activate micromamba from the mount point:

```
cd /path/to/your/mount/point
export MAMBA_ROOT_PREFIX=`pwd`/mambaforge
eval "$(micromamba shell hook --shell bash)"
micromamba activate
micromamba activate tamor
cd tamor
snakemake -j 1 --dry-run
```
