#!/usr/bin/env bash
cd resources
BUNDLE=pcgr_ref_data.grch38.20240621.tgz
wget http://insilico.hpc.uio.no/pcgr/$BUNDLE
tar zxvf $BUNDLE
CACHE=homo_sapiens_vep_112_GRCh38.tar.gz
wget https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/${CACHE}
tar zxvf $CACHE
