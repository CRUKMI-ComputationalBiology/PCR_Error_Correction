#!/usr/bin/env perl

###############################################################################
#   Copyright (C) 2017, Cancer Research UK Manchester Institute
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#############################################################################                       

use strict;
use warnings;
use threads;
no strict qw(subs refs);

#use Statistics::R;
use Class::Struct;

use FindBin;
use lib ("$FindBin::Bin/PerlLib", "$FindBin::Bin/PerlLibAdaptors", "/home/ckim/LIBs/PerlModule/share/perl5/");
use File::Basename;
use File::Path qw(remove_tree rmtree);
use Cwd;
use Carp;
use Getopt::Long qw(:config no_ignore_case pass_through);
use Math::Round qw/round/;

open (STDERR, ">&STDOUT"); 
my $VERSION = "PCR Error Correctiopn (PEC) Algorithm Version 2, 14 Dec. 2017";


############################################################################
#  Configure the MPI command & parameters on your system
#

my $MPI_Command_main = "mpi_command";
my $MPI_Command_misc = "Any parameters for mpi_command";

my $CPU = 1;
my $MR_PAGE_SIZE = 1024;
my $KMER_SIZE = 25;

# option list:
my ($seqType, @left_files, @right_files, @single_files, 
    $output_directory, $outdir_freqTable, $outdir_SNVCall,
    $PE_Single, $flist_reads, $WhichStep 
   );

my %allowed =
    ( seqType       => 'fa or fq',
      PE_single	    => 'PE',
    );

my %allowed_check;
foreach my $all (keys %allowed) {
    my %h = map { (my $s = $_) =~ s/^or //; $s => 1 } split ', ', $allowed{$all};
    $allowed_check{$all} = \%h;
}

# default
#
$output_directory = &create_full_path("pcrError_out_dir");

my $ref_file;
my $ftype;
my $MAQ = 10;
my $BQScore = 20;
my $Overlap_eKmer = 0;
my $run_step=0;
my $TargetRegion;
my $FreqTable;
my $ErrorAllow;
my $AF_allow;
my $Depth_allow;
my $Chromosome;
my $SampleName="SampleName";
my $Soft_Allow;

$PE_Single = "PE";

## Performance monitoring options 
my $pm_logfile = "pcrError.timing";
my $pm_pcrError_start=0;
my $pm_pcrError_end=0;
my $pm_left_fa_size=0;
my $pm_right_fa_size=0;
my $pm_single_fa_size=0;
my $pm_pcrError_fa_size=0;
my $pm_read_count=0;

my $minTumorAlleleFreq  = 0.05;
my $minTumorVariantRead = 3;

my $minHnvRefRead = 10;
my $minHnvVariantRead = 5;
my $maxNormalAlleleFreq = 0.1;
my $maxNormalAlleleFreq_SNV = 0.01;

my $cutoff_log10Pv = 2;
my $alpha_power    = 0.01;

my $Polishing_mode = 1;

#my $start_dir = cwd();

my $usage = <<_EOUSAGE_;

####################################################################
# PCR Error Correction algorithm input parameters   
###################################################################
#
#  --ftype <string>      :type of reads: ( $allowed{seqType} ) 
#  --Ref <string>        :reference sequence in fasta format 
#
#  --output <string>     :output bam filename without ".bam"
#
#---------------------------------------------------------------------------------------
#
#  Input files have to be pair-end datasets upto four lanes
#
#  --input_list <string>  :input filename listing pair-end reads upto 4 lanes, see the github wiki for more detail
#
#-------------------------------------------------------------------------------------
#
#  --CPU <int>           :number of processes to run the algorithm 
#  --page_size <int>     :size of memory (in Mbytes, default is 1024) allocation of MapReduce C++ object
#  		          See documentation of MapReduce-MPI C++ library for more info.
#
#  --KMER_SIZE <int>     :kmer size for sedquence assembly step in PEC algorithm (default is 30)
#
#  --MAQ <int>           :cutoff value (default is 10) of mapping quality from BWA aligner
#
#  --Soft_Allow <float>  :Allowed soft clip proportion (default is 0.1)
#
#  --SM <string>         :sample name for the input dataset
#

_EOUSAGE_

;

my $ROOTDIR                  = "$FindBin::RealBin";
my $FASTA_SPLITTER_DIR       = "$ROOTDIR/Fasta_Splitter_PE";
my $BWA_DIR		     = "$ROOTDIR/BWA";
my $READ_CLUSTERING_DIR      = "$ROOTDIR/PEC_MapReduce";
my $UTILDIR                  = "$ROOTDIR/util";

unless (@ARGV) {
    die "$usage\n";
}

my $FULL_CLEANUP = 0;
my $NO_FASTOOL = 0;

&GetOptions( 

    "output=s" 		=> \$output_directory,

    "input_list=s"       => \$flist_reads,

    "ftype=s"		=> \$ftype,
    "SM=s"		=> \$SampleName,

    "PE_Single=s"     	=> \$PE_Single,

    'CPU=i' 		=> \$CPU,

    'OVERLAP_eKmer=i'	=> \$Overlap_eKmer,

    'KMER_SIZE=i' 	=> \$KMER_SIZE,
    'page_size=i' 	=> \$MR_PAGE_SIZE,
    'MAQ=i'		=> \$MAQ, 
    'BQScore=i'		=> \$BQScore,
    'Target=s'		=> \$TargetRegion,
    'FreqTable=s'	=> \$FreqTable,
    'Ref=s'		=> \$ref_file,
    'ErrorAllow=i'	=> \$ErrorAllow,
    'afAllow=f'		=> \$AF_allow,
    'depthAllow=i'	=> \$Depth_allow,
    'chromosome=i'	=> \$Chromosome,
    'Soft_Allow=f'	=> \$Soft_Allow,
    'Polising_Mode=i'   => \$Polishing_mode,	

);


my $curr_limit_settings = `/bin/sh -c 'ulimit -a' `; 
unless ($curr_limit_settings && $curr_limit_settings =~ /\w/) {
    $curr_limit_settings = `/bin/csh -c limit`; # backup, probably not needed.
}

print "Current settings:\n$curr_limit_settings\n\n";

##################################################################################
#

my $MKDIR_OUTDIR_FLAG = 0;
my $MKDIR_OUTREF_FLAG = 0;

main: {

    $output_directory = &create_full_path($output_directory);
      unless (-d $output_directory) {
           &process_cmd("mkdir -p $output_directory");
           $MKDIR_OUTDIR_FLAG = 1;
    }
      chdir ($output_directory) or die "Error, cannot cd to $output_directory";

      &create_BAM_Freq($flist_reads, $MAQ, $BQScore, $TargetRegion, $ftype, $ErrorAllow, $ref_file, $AF_allow, $Depth_allow, $Soft_Allow, $Polishing_mode);

    exit(0);

}

#####################################################################################

sub create_BAM_Freq {
   my ($flist_Ref, $MAQ, $BQScore, $TargetRegion, $Ftype, $ErrorAllow, $ref_file, $AF_allow, $Depth_allow, $Soft_Allow, $Polishing_mode) = @_;

    my %flist_index;
    my %flist_dir;
    my %flist_R1;
    my %flist_R2;

    my $index = 0;

    open(FP,'<',$flist_Ref) or die "Can't open '$flist_Ref': $!";
    while(<FP>) {
        chomp $_;
        my @arr = split /\s+/, $_;

        if( !(exists $flist_index{$arr[0]}) ) {
                $flist_index{$arr[0]} = 1;
                $flist_dir{$arr[0]}{$arr[1]} = $arr[2];
                $flist_R1{$arr[0]}{$arr[1]}  = $arr[3];
                $flist_R2{$arr[0]}{$arr[1]}  = $arr[4];
                $index++;
        } else {
                $flist_index{$arr[0]} += 1;
                $flist_dir{$arr[0]}{$arr[1]} = $arr[2];
                $flist_R1{$arr[0]}{$arr[1]}  = $arr[3];
                $flist_R2{$arr[0]}{$arr[1]}  = $arr[4];
        }

    }
    close FP;

    foreach my $idx (keys %flist_index ) {

        my $scratch_directory = $output_directory . "/scratch" . "_". $idx;
        my $BamFile = "PEC_outBam" . "_" . $idx;
        my $FreqTable = "BaseFreq.txt" . "." . $idx;

        unless (-d $scratch_directory) {
             &process_cmd("mkdir -p $scratch_directory");
        }

        my @input_R1;
        my @input_R2;

        my $findex = 0;
        foreach my $id (keys %{$flist_dir{$idx}}) {

            if( ($findex >= 0) && ($findex <=3) ) {
                $input_R1[$findex] = $flist_dir{$idx}{$id} . "/" . $flist_R1{$idx}{$id};
                $input_R2[$findex] = $flist_dir{$idx}{$id} . "/" . $flist_R2{$idx}{$id};

                $findex++;
            }

        }

        &run_bwa(\@input_R1, \@input_R2, $scratch_directory, $MAQ, $ftype, $PE_Single, $ref_file);
        &run_readClustering($BamFile, $scratch_directory, $MAQ, $BQScore, $TargetRegion, $BamFile, $ftype, $ErrorAllow, $ref_file, $AF_allow, $Depth_allow, $Soft_Allow, $Overlap_eKmer, $findex, $Polishing_mode, $FreqTable);
        &process_cmd("rm -r $scratch_directory");

    }

    return;




}

####################################################################################

sub run_readClustering {
    my ($outfile, $tmpDIR, $MAQ, $BQScore, $TargetRegion, $SampleName, $Ftype, $ErrorAllow, $ref_file, $AF_allow, $Depth_allow, $Soft_Allow, $Overlap_eKmer, $num_lane, $Polishing_mode, $FreqTable) = @_;

    print "Assigned page size = ", $MR_PAGE_SIZE, "\n";

    my $cmd_readClustering = "$MPI_Command_main $MPI_Command_misc -n $CPU $READ_CLUSTERING_DIR/pec_algorithm --type $Ftype";

    if($MAQ)            { $cmd_readClustering .= " --MAQ $MAQ"; }
    if($BQScore)        { $cmd_readClustering .= " --BQScore $BQScore"; }
    if($MR_PAGE_SIZE)   { $cmd_readClustering .= " --PageSize $MR_PAGE_SIZE"; }
    if($KMER_SIZE)      { $cmd_readClustering .= " --K $KMER_SIZE"; }
    if($SampleName)     { $cmd_readClustering .= " --SM $SampleName"; }
    if($FreqTable)      { $cmd_readClustering .= " --FreqTable $FreqTable"; }
    if($TargetRegion)   { $cmd_readClustering .= " --Target $TargetRegion"; }
    if($ErrorAllow)     { $cmd_readClustering .= " --ErrorAllow $ErrorAllow"; }
    if($ref_file)       { $cmd_readClustering .= " --Ref $ref_file"; }
    if($AF_allow)       { $cmd_readClustering .= " --AF_allow $AF_allow"; }
    if($Depth_allow)    { $cmd_readClustering .= " --Depth_allow $Depth_allow"; }
    if($Soft_Allow)     { $cmd_readClustering .= " --SoftClip $Soft_Allow"; }
    if($Overlap_eKmer)  { $cmd_readClustering .= " --nKmer $Overlap_eKmer"; }  
    if($num_lane)       { $cmd_readClustering .= " --nLane $num_lane"; }
    if($Polishing_mode) { $cmd_readClustering .= " --OnPolishing $Polishing_mode"; }

    $cmd_readClustering .= " --Scratch $tmpDIR";

    &process_cmd($cmd_readClustering);

    if(defined($outfile) && (length $outfile) ) {
        &process_cmd("$UTILDIR/combine_samformat.sh cluster $CPU $outfile $ref_file");
    }

   return;
}


sub run_bwa {
  my ($tmp_input1, $tmp_input2, $tmpDIR, $MAQ, $Ftype, $PE_Single, $ref_file) = @_;

  my @input1 = @{$tmp_input1};
  my @input2 = @{$tmp_input2};

  my @output_read1; 
  my @output_read2;
  for(my $i=0; $i<=$#input1; $i++) {
    $output_read1[$i] = "reads_1_L". $i . ".input";
    $output_read2[$i] = "reads_2_L". $i . ".input";
  }

  my $cmd_splitter;
  my $cmd_alignMPI;

  if($PE_Single eq "PE") {

     for(my $i=0; $i<=$#input1; $i++) {
       if (defined $input1[$i] and length $input1[$i]) { 
	  my $lane_index = $i+1;
	
          my $source_read = $input1[$i];
          my $target_read = $output_read1[$i];

          if($input1[$i] =~ /\.gz/) {
	    &extract_cp_fastq($input1[$i], $output_read1[$i]);
          } else {
            my $cmd_merge_reads = "cp $source_read $target_read";
            &process_cmd($cmd_merge_reads);
          }

          $cmd_splitter = "$MPI_Command_main $MPI_Command_misc -n $CPU $FASTA_SPLITTER_DIR/Fasta_Splitter_PE -t $Ftype";
          $cmd_splitter = $cmd_splitter . " -r $target_read -l $lane_index";
          $cmd_splitter = $cmd_splitter . " -o $tmpDIR -pe 1";
          &process_cmd($cmd_splitter);

          $cmd_splitter = "rm $target_read";
          &process_cmd($cmd_splitter);

	  $source_read = $input2[$i];
          $target_read = $output_read2[$i];

          if($input2[$i] =~ /\.gz/) {
            &extract_cp_fastq($input2[$i], $output_read2[$i]);
          } else {
            my $cmd_merge_reads = "cp $source_read $target_read";
            &process_cmd($cmd_merge_reads);
          }

          $cmd_splitter = "$MPI_Command_main $MPI_Command_misc -n $CPU $FASTA_SPLITTER_DIR/Fasta_Splitter_PE -t $Ftype";
          $cmd_splitter = $cmd_splitter . " -r $target_read -l $lane_index";
          $cmd_splitter = $cmd_splitter . " -o $tmpDIR -pe 2";
          &process_cmd($cmd_splitter);

          $cmd_splitter = "rm $target_read";
          &process_cmd($cmd_splitter); 

      }
    }

  } elsif ($PE_Single eq "Single") {

     for(my $i=0; $i<=$#input1; $i++) {
       if (defined $input1[$i] and length $input1[$i]) {
          my $source_read = $input1[$i];
          my $target_read = $output_read1[$i];

          if($input1[$i] =~ /\.gz/) {
            &extract_cp_fastq($input1[$i], $output_read1[$i]);
          } else {
            my $cmd_merge_reads = "cp $source_read $target_read";
            &process_cmd($cmd_merge_reads);
          }

          $cmd_splitter = "$MPI_Command_main $MPI_Command_misc -n $CPU $FASTA_SPLITTER_DIR/Fasta_Splitter_PE -t $Ftype";
          $cmd_splitter = $cmd_splitter . " -r $target_read";
          $cmd_splitter = $cmd_splitter . " -o $tmpDIR -pe 0";
          &process_cmd($cmd_splitter);

          $cmd_splitter = "rm $target_read";
          &process_cmd($cmd_splitter);
      }
    }

  } else {
   die "Error, the process died with wrong PE file format!!!... Please check your input parameters....";
  }


  for(my $i=0; $i<=$#input1; $i++) {
       if (defined $input1[$i] and length $input1[$i]) {
	my $lane_index = $i+1;
 
        my $cmd_alignMPI = "$MPI_Command_main $MPI_Command_misc -n $CPU $BWA_DIR/BWA -t 1 -f $Ftype -q $MAQ -o $tmpDIR -pe";
  	if(    $PE_Single eq "PE" )    { $cmd_alignMPI = $cmd_alignMPI . " 1 -l $lane_index -r $ref_file";  &process_cmd($cmd_alignMPI); } 
  	elsif ($PE_Single eq "Single") { $cmd_alignMPI = $cmd_alignMPI . " 0 -l $lane_index -r $ref_file";  &process_cmd($cmd_alignMPI); }
  	else { die "Error, the process died with wrong PE file format!!!... Please check your input parameters...."; }

      }
  }

 return;
}


########
sub create_full_path {
    my ($file) = shift;
    if (ref($file) eq "ARRAY"){
       for (my $i=0;$i<scalar(@$file);$i++){
         $file->[$i] = &create_full_path($file->[$i]);
       }
       return @$file;
    }else{
      my $cwd = cwd();
      if ($file !~ m|^/|) { # must be a relative path
          $file = $cwd . "/$file";
      }
      return($file);
    }
}

###
sub process_cmd {
    my ($cmd) = @_;

    print "CMD: $cmd\n";

    my $start_time = time();
    my $ret = system($cmd);
    my $end_time = time();

    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }
    
    print "CMD finished (" . ($end_time - $start_time) . " seconds)\n";    

    return;
}

###
sub extract_cp_fastq {
	my @files = @_;

        my $cmd_merge_reads = "cp $files[0] .";
        &process_cmd($cmd_merge_reads);

        opendir(my $dh, ".") || die "Can't opendir : $!";
        my @fname = grep { /\.fastq.gz/ && -f "./$_" } readdir($dh);
        $cmd_merge_reads = "gunzip $fname[0]";
        &process_cmd($cmd_merge_reads);

        closedir $dh;

        opendir(my $fqdh, ".") || die "Can't opendir : $!";

        my @fqname = grep { /\.fastq/ && -f "./$_" } readdir($fqdh);
        $cmd_merge_reads = "mv $fqname[0] $files[1]";
        &process_cmd($cmd_merge_reads);

        closedir $fqdh;

}
