#!/usr/bin/env perl

# This program rebuilds from a CRAM or BAM the source FASTQ.gz files used as input to the Dragen command that 
# did the mapping, according to the fastq_list.csv provided.
# In contrast to a generic BAM -> FASTQ conversion where all the read 1s and read2s are put into a single pair of files,
# this script regenerates the lane-split reads files, and if the BAM contained reads from multiple sequencing runs, the
# run-splitting is recapitulated by assigning the output FASTQ file for a read based on its read name (which has the run ID embedded in it).
# The upshot of this is that the provided fastq_list.csv can then be used to start another Dragen run, or the FASTQs can be used in other 
# downstream applications as if they just came off the sequencer.

use strict;
use warnings;
use File::Basename;
use IO::Compress::Gzip; # standard since 1998 (Perl 5.005)

@ARGV == 2 or @ARGV == 3 or die "Usage: $0 <input.cram or input.bam> <output_fastq_list.csv> [ref_genome.fa]\n",
                                "Reference genome is required if the input is in CRAM format.\n";

my $aln_src_file = $ARGV[0];
if($aln_src_file =~ /.cram$/){
	if(@ARGV == 2){
		die "FATAL: CRAM input detected based on file extension of $aln_src_file, ",
		    "but no third command-line argument (ref_genome.fa) was provided\n";
	}
	$aln_src_file = "-T $ARGV[2] " . $aln_src_file;
}
my $fastq_list_csv = $ARGV[1];

# First, read the FASTQ list CSV, sanity checking the specs, and then open file handles for all the expected outputs.
open(CSV, $fastq_list_csv)
  or die "Cannot open $fastq_list_csv for reading: $!\n";
# A CSV file could contain multiple lanes and runs, e.g. a primary run supplemented with a top-up run from the same library:
# RGID,RGSM,RGLB,Lane,Read1File,Read2File
# ATATGCATGT.CCAGGCACCA.1,PR-CY-PCA-0006-0002-N,PR-CY-PCA-0006-0002-N,1,novaseq6000/221209_A00906_0334_AHYJGCDSX3/Li37654_S5_L001_R1_001.fastq.gz,novaseq6000/221209_A00906_0334_AHYJGCDSX3/Li37654_S5_L001_R2_001.fastq.gz
# ATATGCATGT.CCAGGCACCA.2,PR-CY-PCA-0006-0002-N,PR-CY-PCA-0006-0002-N,2,novaseq6000/221209_A00906_0334_AHYJGCDSX3/Li37654_S5_L002_R1_001.fastq.gz,novaseq6000/221209_A00906_0334_AHYJGCDSX3/Li37654_S5_L002_R2_001.fastq.gz
# ATATGCATGT.CCAGGCACCA.1,PR-CY-PCA-0006-0002-N,PR-CY-PCA-0006-0002-N,1,novaseq6000/221214_A00906_0339_BHYJH5DSX3/Li37654_S10_L001_R1_001.fastq.gz,novaseq6000/221214_A00906_0339_BHYJH5DSX3/Li37654_S10_L001_R2_001.fastq.gz
# ATATGCATGT.CCAGGCACCA.2,PR-CY-PCA-0006-0002-N,PR-CY-PCA-0006-0002-N,2,novaseq6000/221214_A00906_0339_BHYJH5DSX3/Li37654_S10_L002_R1_001.fastq.gz,novaseq6000/221214_A00906_0339_BHYJH5DSX3/Li37654_S10_L002_R2_001.fastq.gz

my $header = <CSV>;
chomp $header;
my @header_columns = split /,/, $header;
my $num_columns = scalar(@header_columns);
my $Read1FileColumnIndex;
my $Read2FileColumnIndex;
my $LaneColumnIndex;
for my $index (0..$#header_columns){
	if($header_columns[$index] eq "Read1File"){
		$Read1FileColumnIndex = $index
	}
	elsif($header_columns[$index] eq "Read2File"){
		$Read2FileColumnIndex = $index
	}
	elsif($header_columns[$index] eq "Lane"){
		$LaneColumnIndex = $index
	}
}
if(not defined $Read1FileColumnIndex){
	die "FATAL: Could not find the minimum required file path header column (Read1File) in the FASTQ list CSV file $fastq_list_csv\n";
}
if(not defined $LaneColumnIndex){
	die "FATAL: Could not find the required 'Lane' header column in the FASTQ list CSV file $fastq_list_csv\n";
}

# Key is flowcell:lane:[12] -> filepath
my %flowcell_lane_read2fastq_path;
my %flowcell_lane_read2fastq_encoding;
while(<CSV>){
	chomp;
	my @F = split /,/, $_;
	if($#F != $num_columns-1){
		die "FATAL: Found FASTQ list CSV row with ", scalar(@F), " columns, but expected $num_columns\n";
	}
	my $lane = $F[$LaneColumnIndex];
	my $fastq_r1 = $F[$Read1FileColumnIndex];
	my ($run_dir, $suffix) = $fastq_r1 =~ m(([^/]+)/[^/]+(\.fastq|\.fastq\.gz|\.ora)$);
	if(not defined $run_dir){
		die "FATAL: Could not find a possible run name in the directory path part of $fastq_r1, please ensure the FASTQ is nested in directory\n";
	}
	if(not defined $suffix){
		die "FATAL: Could not find an acceptable suffix (.fastq, .fastq.gz, or .ora) on $fastq_r1, conversion cannot proceed\n";
	}
	# The flowcell serial (9 letters and numbers) is embedded at the end of the run name.
	my ($flowcell) = $run_dir =~ /\d\d\d\d_(?:[AB])?([A-Z0-9]{9})$/;
	if(not defined $flowcell){
		die "FATAL: Cannot find a nine character Illumina flowcell serial number in the FASTQ directory path element $run_dir (for FASTQ file $fastq_r1)\n"
	}
	$flowcell_lane_read2fastq_path{"$flowcell:$lane:1"} = $fastq_r1;
	$flowcell_lane_read2fastq_encoding{"$flowcell:$lane:1"} = $suffix;

	if(defined $Read2FileColumnIndex){
		my $fastq_r2 = $F[$Read2FileColumnIndex];
		# Sanity check that the second read in the pair comes from the same flowcell and is compressed the same way.
		my ($run_dir2, $suffix2) = $fastq_r2 =~ m(([^/]+)/[^/]+(\.fastq|\.fastq\.gz|\.ora)$);
       	 	if(not defined $run_dir2){
       	         	die "FATAL: Could not find a possible run name in the directory path part of $fastq_r2, please ensure the FASTQ is nested in directory ($fastq_list_csv line $.)\n";
        	}
		if($run_dir2 ne $run_dir){
			die "FATAL: The output directory for both files of a read pair must be the same, but found $run_dir2 != $run_dir ($fastq_list_csv line $.)\n";
		}
       	 	if(not defined $suffix2){
       	         	die "FATAL: Could not find an acceptable suffix (.fastq, .fastq.gz, or .ora) on $fastq_r2, conversion cannot proceed ($fastq_list_csv line $.)\n";
       	 	}
		if($suffix2 ne $suffix){
			die "FATAL: The suffix for both files of a read pair must be the same, but found $run_dir2 != $run_dir ($fastq_list_csv line $.)\n";
		}
		$flowcell_lane_read2fastq_path{"$flowcell:$lane:2"} = $fastq_r2;
		$flowcell_lane_read2fastq_encoding{"$flowcell:$lane:2"} = $suffix2;
	}
}
close(CSV);

# Make sure all the files do NOT already exist, and are writeable before we decide to proceed.
for my $outfile (values %flowcell_lane_read2fastq_path){
	if(-e $outfile){
		die "FATAL: Expected output file $outfile exists, will not overwrite. Manually move or delete it to reconvert.\n";
	}
	my $outdir = dirname($outfile);
	if(not -w $outdir){
		die "FATAL: Output directory $outdir is not writeable.\n";
	}
}

# Now we can finally open all the output file handles.
my %flowcell_lane_read2fastq_filehandle;
for my $key (keys %flowcell_lane_read2fastq_path){
	my $encoding = $flowcell_lane_read2fastq_encoding{$key};
	my $filepath = $flowcell_lane_read2fastq_path{$key};
	my $fh;
	open($fh, ">$filepath")
	  or die "Cannot open $filepath for writing: $!\n";
  	# Will output gzip'ed file.
	if($encoding eq ".fastq.gz" or $encoding eq ".ora"){
		my $gz_fh = IO::Compress::Gzip->new($fh)
                  or die "Cannot create IO::Compress::Gzip filehandle: $!\n";
	  	$flowcell_lane_read2fastq_filehandle{$key} = $gz_fh;
	}
	# Plain FASTQ output.
	else{
		$flowcell_lane_read2fastq_filehandle{$key} = $fh;
	}
}

# View everything except secondary and supplemental alignments (should be all primary and unmapped).
open(SAM, "samtools view -F 0x900 $aln_src_file |")
  or die "Cannot run samtools: $!\n";
while(<SAM>){
	my @F = split /\t/, $_;
	my @id_fields = split /:/, $F[0];
	if(scalar(@id_fields) != 7){
		die "FATAL: Could not parse read ID $F[0] in the BAM/CRAM, does not conform to Illumina read ID format ",
		    "suited to rebuilding FASTQs (expected Instrument:RunID:FlowCellID:Lane:Tile:X:Y)\n";
	}
	# Assume it's read1 if not explicitly read2 in the SAM flags.
	my $is_read2 = int($F[1]) & 0x80; 
	my $key = "$id_fields[2]:$id_fields[3]:".($is_read2 ? "2" : "1");
	if(not exists $flowcell_lane_read2fastq_filehandle{$key}){
		die "FATAL: Read ID $F[0] in the input BAM/CRAM references a flowcell+lane combination ($id_fields[2]+$id_fields[3]) that ",
		    "was not specified in the provided FASTQ list CSV. Aborting, as the conversion would be incomplete. Please amend $fastq_list_csv accordingly.\n";
	}
	$flowcell_lane_read2fastq_filehandle{$key}->print("\@$F[0]\n$F[9]\n+\n$F[10]\n");
}
close(SAM);

for my $fh (values %flowcell_lane_read2fastq_filehandle){
	close($fh);
}
