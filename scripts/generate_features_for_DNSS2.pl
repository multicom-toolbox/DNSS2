#!/usr/bin/perl
# Script: 	P1_generate_features.pl
# Author: 	Jie Hou
# Made:		11/03/17
#
#
# Input:
#
# Dependencies:

my $GLOBAL_PATH;
BEGIN { $GLOBAL_PATH = '/storage/htc/bdm/Collaboration/jh7x3/DeepCov_SS_SA_project/DNSS2.0/'; }
use lib "$GLOBAL_PATH/lib/";
use strict;
use Time qw(formatted_localtime);
use DN_SSpred2 qw(generate_feature_for_convolution_DNSS2_exp timeprint sort_scores score_dnss check_err); 
use Getopt::Long;

my ($help, $replace, $outdir, $wind, $iters, $bound, $reduced, $boost, $skipx, $arch, $pssm, $pssm_uniref50, $hmmEm, $hmmTr, $aln_sspro, $aln_sspro_unif90, $aln_sspro_unif50, $aln_hhblits, $atch, $seq, $tag); 
my $opt = GetOptions(	'-r!', \$replace,
			'-h!', \$help,
			'-help!', \$help,
			'-out:s', \$outdir,
			'-wind:i', \$wind,
			'-window:i', \$wind, 
			'-w:i', \$wind,
			'-bound!', \$bound,
			'-red!', \$reduced,
			'-x:i', \$skipx,
			'-boost:i', \$boost,
			'-iter:i', \$iters,
			'-pssm:i', \$pssm,
			'-pssm_uniref50:i', \$pssm_uniref50,
			'-hmmEm:i', \$hmmEm,
			'-hmmTr:i', \$hmmTr,
			'-aln_sspro:i', \$aln_sspro,
			'-aln_sspro_unif90:i', \$aln_sspro_unif90,
			'-aln_sspro_unif50:i', \$aln_sspro_unif50,
			'-aln_hhblits:i', \$aln_hhblits,
			'-atch:i', \$atch,
			'-seq:i', \$seq,
			'-tag:s', \$tag,
			'-arch:s', \$arch );

if ($help || !$outdir){
	print_help();
	exit;
}

$outdir .= "/" unless ($outdir =~ m/\/$/);

#####################################
$iters = 1 unless ($iters);
$pssm = 1 unless (defined($pssm));
$pssm_uniref50 = 1 unless (defined($pssm_uniref50));
$hmmEm = 1 unless (defined($hmmEm));
$hmmTr = 1 unless (defined($hmmTr));
$aln_sspro = 1 unless (defined($aln_sspro));
$aln_sspro_unif90 = 1 unless (defined($aln_sspro_unif90));
$aln_sspro_unif50 = 1 unless (defined($aln_sspro_unif50));
$aln_hhblits = 1 unless (defined($aln_hhblits));
$atch = 1 unless (defined($atch));
$seq = 1 unless (defined($seq)); ## DNSS exclude this, but I add here first
$wind = 7 unless ($wind);
$boost = 1 unless ($boost);
$skipx = 0 unless ($skipx);
$bound = 1 unless (defined($bound));
$reduced = 0 unless ($reduced);
my $target = 3*$boost;
$arch = "X,X,$target" unless ($arch);
my @archinfo = split(/,/, $arch);
#my $lwind = 17;
#my $larch = "600,600,600,250,3";
#####################################

my $trainfile = "/storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/datasets/NewTrainTest_20181027/data/lists/dnss2_train.lst";
my $testfile = "/storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/datasets/NewTrainTest_20181027/data/lists/dnss2_val.lst";
my $casp9file = "/storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/datasets/NewTrainTest_20181027/data/lists/casp9.lst";
my $casp10file = "/storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/datasets/NewTrainTest_20181027/data/lists/casp10.lst";
my $ssadir = "/storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/datasets/NewTrainTest_20181027/data/ssa/";
#my $pssmdir = "$GLOBAL_PATH/datasets/pssm/";
my $pssmdir = "/storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/datasets/feature_construct_cullpdb/uniref90_pssm/";
my $pssm_uniref50dir = "/storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/datasets/feature_construct_cullpdb/uniref90_pssm/";
my $hmmEmdir = "/storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/datasets/feature_construct_cullpdb/uniclust30_2017_10_hmm/";
my $hmmTrdir = "/storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/datasets/feature_construct_cullpdb/uniclust30_2017_10_hmm/";
my $msa2profile_ssprodir = "/storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/datasets/feature_construct_cullpdb/uniclust30_2017_10_hmm/";
my $msa2profile_sspro_uniref90_dir = "/storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/datasets/feature_construct_cullpdb/uniclust30_2017_10_hmm/";
my $msa2profile_sspro_uniref50_dir = "/storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/datasets/feature_construct_cullpdb/uniclust30_2017_10_hmm/";
my $msa2profile_hhblitsdir = "/storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/datasets/feature_construct_cullpdb/uniclust30_2017_10_hmm/";
my $overdir = "/storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/datasets/NewTrainTest_20181027/data/feature_cullpdb/"; ## this is jie's working directory 
$outdir = $overdir . $outdir;
my $logfile = $outdir . "log.txt";

my @features = ($pssm,$pssm_uniref50, $hmmEm, $hmmTr, $aln_sspro, $aln_sspro_unif50, $aln_sspro_unif90, $aln_hhblits, $atch, $seq, $bound, $skipx);
my @options = (\@features, $wind, $boost, $reduced);
my @dirs = ($ssadir, $pssmdir, $pssm_uniref50dir, $hmmEmdir, $hmmTrdir, $msa2profile_ssprodir, $msa2profile_sspro_uniref50_dir, $msa2profile_sspro_uniref90_dir, $msa2profile_hhblitsdir);

`mkdir $overdir` unless (-d $overdir);
`mkdir $outdir` unless (-d $outdir);

timeprint($logfile, "Performing SS train & test using the parameters:\nTag: $tag\nOutdir: $outdir\nIterations: $iters\nReduced: $reduced\nReplace: $replace\nSkip X: $skipx\n\nWindow: $wind\nPssm: $pssm\nPssm_uniref50: $pssm_uniref50\nHMMEM: $hmmEm\nHMMTR: $hmmTr\naln_sspro: $aln_sspro\naln_sspro_unif50: $aln_sspro_unif50\naln_sspro_unif90: $aln_sspro_unif90\naln_hhblits: $aln_hhblits\nAtch: $atch\nSeq: $seq\nArch: $arch\nBoost: $boost\nBound: $bound\n");

`rm -r $outdir/*` if ($replace);


my $protlist = `cat $trainfile`;
my @trainlist = split (/\n/, $protlist);
my $protlist = `cat $testfile`;
my @testlist = split (/\n/, $protlist);

generate_feature_for_convolution_DNSS2_exp(	$logfile, $outdir, \@trainlist, \@testlist,	\@dirs, \@options);


#my $protlist = `cat $casp9file`;
#my @trainlist = split (/\n/, $protlist);
#my $protlist = `cat $casp10file`;
#my @testlist = split (/\n/, $protlist);

#generate_feature_for_convolution_DNSS2_exp(	$logfile, $outdir, \@trainlist, \@testlist,	\@dirs, \@options);



################################################################
################################################################

sub print_help {
	print "\nHelp Summary for P1_generate_features.pl script\n";
	print "Written by Jie Hou\n";
	print "Made:		11/03/17\n";
	print "\nDescription:\n";
	print "This script generates training and testing features for secondary structure prediction.\n";
	print "\nRequired input:\n";
	print "\t-out	: Subdirectory of Deep1Dconv_ss/features_win1/ to save intermediate files to.\n";
	print "\t-wind	: Indicate the window size to use.\n";
	print "\nOptions:\n";
	print "\t-r	: Indicates that previously existing files will be replaced.\n";
	print "\t-x	: Deal with X res (0: keep, 1: skip lines, 2: skip windows)\n";
	print "\t-red	: Only predict non-gap residues.\n";
	print "\t-boost	: Indicate the size of boost window.\n";
	print "\t-pssm	: Include pssm features.\n";
	print "\t-hmmEm	: Include hmm emission features.\n";
	print "\t-hmmTr	: Include hmm transition.\n";
	print "\t-atch	: Include atchley factors.\n";
	print "\t-seq	: Include seq features.\n";
	print "\t-bound	: Include boundary features.\n";
	print "\t-arch	: Indicate the architecture to use.\n";
	print "\t-iter	: Indicate the number of iterations to attempt.\n";
	print "\t-help	: Print this help message.\n";
	print "\n\n";
}
