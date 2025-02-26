# script to visualize germline and somatic cnv calls and supporting BAFs/depth data for a moh sample 
# mamba env export -n karyoploter > workflow/envs/karyoploter.yaml

### TODO add .dna.germline.improper.pairs.bw and .dna.somatic.tumor.improper.pairs.bw counts
### TODO add .dna.germline.seg.called.merged; .dna.somatic.baf.seg; .dna.somatic.seg; .dna.somatic.seg.tumorDepths
### see https://help.dragen.illumina.com/product-guides/dragen-v4.3/dragen-dna-pipeline/cnv-calling/cnv-output

# ref https://bioconductor.org/packages/release/bioc/vignettes/CopyNumberPlots/inst/doc/CopyNumberPlots.html#plotcopynumbercalls

rm(list = ls())

library(rtracklayer)
library(karyoploteR)
library(CopyNumberPlots)
library(regioneR)
library(vcfR)
library(stringr)
library(zoo)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#all.genes <- genes(txdb) 

args <- commandArgs(trailingOnly = TRUE)
print(args)

project = args[1]
subject = args[2]
tumor = args[3]
normal = args[4]
output_dir = args[5]

dna_germline_path = paste0(output_dir,'/',project,'/',subject,'/',subject,'_',normal)
dna_somatic_path = paste0(output_dir,'/',project,'/',subject,'/',subject,'_',tumor,'_',normal)


# germline normalized depth file
bw_file <- paste0(dna_germline_path, ".dna.germline.tn.bw") 
if (file.exists(bw_file)) {
  germline.tn.bw <- import.bw(bw_file)
  names(mcols(germline.tn.bw)) <- c("TN_depth")
} else {
  germline.tn.bw <- data.frame("chr1",1,1,0)
  names(germline.tn.bw) <- c("chr","start","end","TN_depth")
  germline.tn.bw <- toGRanges(germline.tn.bw)
}

# germline bw file
bw_file <- paste0(dna_germline_path, ".dna.germline.hard-filtered.baf.bw")
if (file.exists(bw_file)) {
  germline.baf.bw <- import.bw(bw_file)
  names(mcols(germline.baf.bw)) <- c("BAF")
} else {
  germline.baf.bw <- data.frame("chr1",1,1,0)
  names(germline.baf.bw) <- c("chr","start","end","BAF")
  germline.baf.bw <- toGRanges(germline.baf.bw)
}

# germline PASS cnv calls
vcf_file <- paste0(dna_germline_path, ".dna.germline.cnv.vcf.gz")
if (file.exists(vcf_file)) {
  cnv_file <- read.vcfR(vcf_file)
  chr <- getCHROM(cnv_file)
  start <- as.numeric(getPOS(cnv_file))
  id <- getID(cnv_file)
  end <- as.numeric(sapply(strsplit(id, "-"), "[", 2))
  loh <- str_detect(sapply(strsplit(id, ":"), "[", 2), "LOH")
  filter <- getFILTER(cnv_file)
  copyNumber <- as.numeric(extract.gt(cnv_file, element = "CN", IDtoRowNames = FALSE))
  cnvs <- data.frame(chr,start,end,filter,copyNumber,loh)
  pass.cnvs <- cnvs[cnvs$filter == "PASS", ]
  germline.cnv.gr <- toGRanges(data.frame(chr=pass.cnvs["chr"], start=pass.cnvs["start"], end=pass.cnvs["end"], cn=pass.cnvs["copyNumber"], loh=pass.cnvs["loh"]))
} else {
  germline.cnv.gr <- data.frame("chr1",1,1,2,FALSE)
  names(germline.cnv.gr) <- c("chr","start","end","copyNumber","loh")
  germline.cnv.gr <- toGRanges(germline.cnv.gr)
}

#### somatic files

# somatic normalized depth file
bw_file <- paste0(dna_somatic_path, ".dna.somatic.tn.bw")
if (file.exists(bw_file)) {
  somatic.tn.bw <- import.bw(bw_file)
  names(mcols(somatic.tn.bw)) <- c("TN_depth")
} else {
  somatic.tn.bw <- data.frame("chr1",1,1,0)
  names(somatic.tn.bw) <- c("chr","start","end","TN_depth")
  somatic.tn.bw <- toGRanges(somatic.tn.bw)
}

# somatic bw file
bw_file <- paste0(dna_somatic_path, ".dna.somatic.hard-filtered.baf.bw")
somatic.baf.bw <- import.bw(bw_file)
names(mcols(somatic.baf.bw)) <- c("BAF")

# baf in tumor of germline het snv's
bedgraph <- paste0(dna_somatic_path, ".dna.somatic.tumor.baf.bedgraph.gz")
if (file.exists(bedgraph)) {
  tumor.baf.bed <- read.table(bedgraph, sep = '\t', header = F, col.names = c("chr","start","end","BAF"))
  tumor.baf.gr <- toGRanges(data.frame(chr=tumor.baf.bed["chr"], start=tumor.baf.bed["start"], end=tumor.baf.bed["end"], baf=tumor.baf.bed["BAF"]))
} else {
  tumor.baf.gr <- data.frame("chr1",1,1,0)
  names(tumor.baf.gr) <- c("chr","start","end","BAF")
  tumor.baf.gr <- toGRanges(tumor.baf.gr)
}

# vaf of passing somatic snv's
vafpass <- paste0(dna_somatic_path, ".dna.somatic.hard-filtered.pass.vaf")
if (file.exists(vafpass)) {
  tumor.vaf.bed <- read.table(vafpass, sep = '\t', header = F, col.names = c("chr","start","end","BAF"))
  tumor.vaf.gr <- toGRanges(data.frame(chr=tumor.vaf.bed["chr"], start=tumor.vaf.bed["start"], end=tumor.vaf.bed["end"], baf=tumor.vaf.bed["BAF"]))
} else {
  print("No VAF file found")
  tumor.vaf.gr <- data.frame("chr1",1,1,0)
  names(tumor.vaf.gr) <- c("chr","start","end","BAF")
  tumor.vaf.gr <- toGRanges(tumor.vaf.gr)
}

# somatic PASS cnv calls
vcf_file <- paste0(dna_somatic_path, ".dna.somatic.cnv.vcf.gz")
if (file.exists(vcf_file)) {
  cnv_file <- read.vcfR(vcf_file)
  chr <- getCHROM(cnv_file)
  start <- as.numeric(getPOS(cnv_file))
  id <- getID(cnv_file)
  end <- as.numeric(sapply(strsplit(id, "-"), "[", 2))
  loh <- str_detect(sapply(strsplit(id, ":"), "[", 2), "LOH")
  filter <- getFILTER(cnv_file)
  copyNumber <- as.numeric(extract.gt(cnv_file, element = "CN", IDtoRowNames = FALSE))
  cnvs <- data.frame(chr,start,end,filter,copyNumber,loh)
  pass.cnvs <- cnvs[cnvs$filter == "PASS", ]
  somatic.cnv.gr <- toGRanges(data.frame(chr=pass.cnvs["chr"], start=pass.cnvs["start"], end=pass.cnvs["end"], cn=pass.cnvs["copyNumber"], loh=pass.cnvs["loh"]))
} else {
  somatic.cnv.gr <- data.frame("chr1",1,1,2,FALSE)
  names(somatic.cnv.gr) <- c("chr","start","end","copyNumber","loh")
  germline.cnv.gr <- toGRanges(somatic.cnv.gr)
}


#### create plot

jpeg(filename=paste0(dna_somatic_path, ".cnv_plot.jpeg"), width=30, height=20, units="in", res=750)
kp <- plotKaryotype(plot.type = 4, genome= "hg38")
kpAddMainTitle(kp, main=paste0(subject,'_',tumor,'_',normal), cex=2)

# Plot germline normalized depth signal
tn.ymax <- 3
tn.ymin <- -3
kpPoints(kp, data=germline.tn.bw, y=germline.tn.bw$TN_depth, cex = 0.1, ymax=tn.ymax, ymin=tn.ymin, r0=0.86, r1=1, col = "darkslategray4")
kpAxis(kp, numticks = 5, r0=0.86, r1=1, tick.len = 5e6, side=2, ymax=tn.ymax, ymin=tn.ymin)
kpAddLabels(kp, labels="Germline Norm Depth", r0=0.86, r1=1, data.panel = 1, srt=90, pos=3, label.margin = 0.03)
kpAbline(kp, h=0.5, col="red", r0=0.86, r1=1)

# Plot germline BAF
kpPoints(kp, data=germline.baf.bw, y=germline.baf.bw$BAF, cex=0.1, r0=0.7, r1=0.84)
kpAxis(kp, numticks = 5, r0=0.7, r1=0.84, tick.len = 5e6, side=2)
kpAddLabels(kp, labels="Germline BAF", r0=0.7, r1=0.84, data.panel = 1, srt=90, pos=3, label.margin = 0.03)

# Plot germline cnvs
plotCopyNumberCalls(kp, germline.cnv.gr, cn.column = "copyNumber", cn.colors="red_blue", r0=0.62, r1=0.66, loh.height = 0, labels = NA)
kpAddLabels(kp, labels="Germline Copy Number", r0=0.6, r1=0.68, data.panel = 1, srt=90, pos=3, label.margin = 0.03)

# Plot somatic normalized depth
tn.ymax <- 3
tn.ymin <- -3
kpPoints(kp, data=somatic.tn.bw, y=somatic.tn.bw$TN_depth, cex = 0.1, ymax=tn.ymax, ymin=tn.ymin, r0=0.42, r1=0.56, col = "darkslategray4")
kpAxis(kp, numticks = 5, r0=0.42, r1=0.56, tick.len = 5e6, side=2, ymax=tn.ymax, ymin=tn.ymin)
kpAddLabels(kp, labels="Somatic Norm Depth", r0=0.42, r1=0.56, data.panel = 1, srt=90, pos=3, label.margin = 0.03)
kpAbline(kp, h=0.5, col="red", r0=0.42, r1=0.56)

# Plot somatic BAF
kpPoints(kp, data=somatic.baf.bw, y=somatic.baf.bw$BAF, cex=0.3, r0=0.26, r1=0.4)
kpAxis(kp, numticks = 5, r0=0.26, r1=0.4, tick.len = 5e6, side=2)
kpAddLabels(kp, labels="Somatic BAF", r0=0.26, r1=0.4, data.panel = 1, srt=90, pos=3, label.margin = 0.03)

# Add label for somatic VAF
kpAddLabels(kp, labels="Mean VAF", r0=0.26, r1=0.4, data.panel = 1, srt=90, pos=3, label.margin = 0.02, col="#FF0000B3")

# Calculate and plot rolling mean VAF of somatic SNV's (30 is a nice bin size but breaks if fewer in chr)  
for(chr in seqlevels(kp$genome)) {
  chr.dp <- sort(keepSeqlevels(x = tumor.vaf.gr, value = chr, pruning.mode = "coarse"))
  if (length(chr.dp) > 30) {
    rmean <- rollmean(chr.dp$BAF, k = 30, align = "center")  
    kpLines(kp, chr=chr, x=start(chr.dp)[1:(length(chr.dp))], y=rmean, col="#FF000099", r0=0.26, r1=0.4)
  } else {
    rmean <- rollmean(chr.dp$BAF, k = length(chr.dp), align = "center")
    kpLines(kp, chr=chr, x=start(chr.dp)[1:(length(chr.dp))], y=rmean, col="#FF000099", r0=0.26, r1=0.4)
  }
}

# Plot tumor BAF at germline Het sites
kpPoints(kp, data=tumor.baf.gr, y=tumor.baf.gr$BAF, cex=0.1, r0=0.1, r1=0.24, col="cornflowerblue")
kpAxis(kp, numticks = 5, r0=0.1, r1=0.24, tick.len = 5e6, side=2)
kpAddLabels(kp, labels="Somatic BAF of\nGermline Het SNVs", r0=0.1, r1=0.24, data.panel = 1, srt=90, pos=3, label.margin = 0.03)

# Plot somatic cnvs and loh events
plotCopyNumberCalls(kp, somatic.cnv.gr, cn.column = "copyNumber", cn.colors="red_blue", r0=0, r1=0.06, labels = NA, loh.color = "coral")
kpAddLabels(kp, labels="Somatic Copy Number", r0=0, r1=0.08, data.panel = 1, srt=90, pos=3, label.margin = 0.03)
kpAddLabels(kp, labels="LOH", r0=0, r1=0.02, data.panel = 1, pos=3, side = "right")

# Add cnv legend
cn.cols <- getCopyNumberColors(colors = "red_blue")
legend("topright", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols), title="Copy Number")

dev.off()

warnings()
