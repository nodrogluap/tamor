# Go to where Tamor expects the FASTQs in its default config.
cd results/analysis/primary/HiSeq

# Use the SRA tools (install in the test conda env) to fetch a publicly available sequencing dataset.
# Short Read Archive retrieval of sequencing data by accession - CLL DNA
prefetch SRR6702602

# Short Read Archive retrieval of sequencing data by accession - matched CLL RNA
prefetch SRR6702601

# Unpack the archive data into compressed FASTQ sequence/quality files (this also takes a long time as the raw sequences are reconstructed as a diff relative to the reference human genome):
fastq-dump -F --gzip --split-3 SRR6702602
fastq-dump -F --gzip --split-3 SRR6702601

# There is no "normal" sample here, only tumor. 
# So generate a normal-ish input, by removing most reads with somatic mutations from the FASTQ files, using a k-mer table of 30 million common nucleotide 
# variants in the human population, and the human reference genome k-mers (k=15):

# Get the common human population variant kmers (30mers) from the FastGT tool.
wget http://bioinfo.ut.ee/FastGT/downloads/kmer_db_WG30238282.db

# Calculate the reference genome k-mers with the Jellyfish program (installed in the test conda env). Needs ~40Gb memory.
jellyfish count -m 15 -s 100M -t 16 -C -o human.15mers.jf /work/chgi_common/reference_data/fasta/hg19.fa
# Export the internal Jellyfish file formatted kmer table to text, capturing just the kmers but not their counts.
jellyfish dump -c human.15mers.jf | cut -f 1 > human.15mers.txt

# Cut the variant kmers down to all 15mers we don't have already from the reference genome.
perl -F\\t -ane 'BEGIN{%h=split /(\s)/, `cat human.15mers.txt`}for my $kmer (@F[2..$#F]){for(my $i = 0; $i+15 < length($kmer);$i++){$k = substr($kmer, $i, 15); print "$k\n" unless exists $h{$k} or $p{$k}++}};' kmer_db_WG30238282.db > snp15

# Keep only reads without any uniq 15mers. The following need about 80GB RAM.
gzip -cd SRR6702602_2.fastq.gz | perl -ne 'BEGIN{open(H,"cat human.15mers.txt snp15|"); while(<H>){chomp; $kmer{$_}=1; tr/ACGT/TGCA/; $_ = reverse($_); $kmer{$_} = 1};$/="\n\@"}$uniq = 0; ($seq) = /^[^\n]+\n([^\n]+)/; while($seq=~/.{15}/g){if(not exists $kmer{$&}){$uniq = 1; last}}if(not $uniq){s/^\@// if $. == 1; chomp; print "\@$_\n"}' | gzip > SRR6702602_2.nonsomatic.fastq.gz
gzip -cd SRR6702602_1.fastq.gz | perl -ne 'BEGIN{open(H,"cat human.15mers.txt snp15|"); while(<H>){chomp; $kmer{$_}=1; tr/ACGT/TGCA/; $_ = reverse($_); $kmer{$_} = 1};$/="\n\@"}$uniq = 0; ($seq) = /^[^\n]+\n([^\n]+)/; while($seq=~/.{15}/g){if(not exists $kmer{$&}){$uniq = 1; last}}if(not $uniq){s/^\@// if $. == 1; chomp; print "\@$_\n"}' | gzip > SRR6702602_1.nonsomatic.fastq.gz 

# Only include read pairs where both reads appear to be not have any somatic mutations.
gzip -cd SRR6702602_1.nonsomatic.fastq.gz | perl -ne 'print "$1\n" if /^\@(\d+)/' > r1
gzip -cd SRR6702602_2.nonsomatic.fastq.gz | perl -ne 'print "$1\n" if /^\@(\d+)/' > r2
sort r1 r2 | uniq -d > keeper_read_ids
mkdir SRX3676780
gzip -cd SRR6702602_1.nonsomatic.fastq.gz | perl -ne 'BEGIN{%keep = split /(\n)/s, `cat keeper_read_ids`; $/="\n\@"} if(/^\@?(\d+)/ and $keep{$1}){s/^\@// if $. == 1; chomp; print "\@$_\n"} ' | gzip -c - > SRX3676780/SRR6702602_1.pseudonormal.fastq.gz
gzip -cd SRR6702602_2.nonsomatic.fastq.gz | perl -ne 'BEGIN{%keep = split /(\n)/s, `cat keeper_read_ids`; $/="\n\@"} if(/^\@?(\d+)/ and $keep{$1}){s/^\@// if $. == 1; chomp; print "\@$_\n"} ' | gzip -c - > SRX3676780/SRR6702602_2.pseudonormal.fastq.gz

# Clean up large intermediate files
rm SRR6702602_1.nonsomatic.fastq.gz SRR6702602_2.nonsomatic.fastq.gz r1 r2 keeper_read_ids
# Place the new generated sample FASTQs (SRR*) into their respective run folders (SRX*) as specified in the default config/dna_samples.tsv and config/rna_samples.tsv specs
# RNA
mkdir SRX3676781
mv SRR6702601_1.fastq.gz SRR6702601_2.fastq.gz SRX3676781
# DNA
mv SRR6702602_1.fastq.gz SRR6702602_2.fastq.gz SRX3676780

# Generate a fastq_list.csv file like bcl-convert normally would, as this file is used by Tamor to pass FASRQ file paths along to Dragen. Absolute paths are preferred by Dragen.
PWD=`pwd`
mkdir SRX3676781/Reports
echo "RGID,RGSM,RGLB,Lane,Read1File,Read2File" >> SRX3676781/Reports/fastq_list.csv
echo "GCTAGACTAT.ATGTCGTATT.1,SRR6702601,UnknownLibrary,1,$PWD/SRX3676781/SRR6702601_1.fastq.gz,$PWD/SRX3676781/SRR6702601_2.fastq.gz" >> SRX3676781/Reports/fastq_list.csv
mkdir SRX3676780/Reports
echo "RGID,RGSM,RGLB,Lane,Read1File,Read2File" >> SRX3676780/Reports/fastq_list.csv
echo "GCTAGACTAT.ATGTCGTATT.1,SRR6702602,UnknownLibrary,1,$PWD/SRX3676780/SRR6702602_1.fastq.gz,$PWD/SRX3676780/SRR6702602_2.fastq.gz" >> SRX3676780/Reports/fastq_list.csv
echo "ACTTCAAGCG.TTCATGGTTC.1,scrubbed-SRR6702602,UnknownLibrary,1,$PWD/SRX3676780/SRR6702602_1.pseudonormal.fastq.gz,$PWD/SRX3676780/SRR6702602_2.pseudonormal.fastq.gz" >> SRX3676780/Reports/fastq_list.csv
