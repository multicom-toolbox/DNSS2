#!/usr/bin/perl
# Script: 	sstest.pl
# Author: 	Jie Hou
# Made:		9/25/13
# Last Mod:	11/4/17
#
# Input:
#
# Dependencies:
my $GLOBAL_PATH;
BEGIN {
$GLOBAL_PATH='/faculty/jhou4/tools/DNSS2/';
}

use strict;
use lib "$GLOBAL_PATH/lib";
use strict;
use Time qw(formatted_localtime);
use DN_SSpred2 qw(newDN timeprint sort_scores score_dnss check_err all_dnss_files score_server_prediction); 
use Getopt::Long;

my ($help,  $predir, $testfile, $tag,$server,$eva_ext,$seq_ext,$ssa);
my $opt = GetOptions('-h!', \$help,
			'-help!', \$help,
			'-indir:s', \$predir,
			'-list:s', \$testfile,
			'-ssa:s', \$ssa,
			'-server:s', \$server,
			'-eva_ext:s', \$eva_ext,
			'-seq_ext:s', \$seq_ext,
			'-tag:s', \$tag);

if ($help || !$tag || !$predir  || !$ssa || !$testfile){
	print_help();
	exit;
}

$predir .= "/" unless ($predir =~ m/\/$/);
$ssa .= "/" unless ($ssa =~ m/\/$/);

my $logfile = $predir . "log.txt";
#####################################

`mkdir $predir` unless (-d $predir);

my @Qscores;
my @Sovscores;

my $thistag = $tag;

my $protlist = `cat $testfile`;
my @testlist = split (/\n/, $protlist);


################################################################
##	Generating accuracy
################################################################

# Probabilities are reformatted into usable file aligning the sequence
# with the corresponding predicted most likely secondary structure
my @dirarray = ($predir);

timeprint ($logfile, "Evaluating predictions...");
my ($AvgQ3, $AvgSov) = score_server_prediction($predir, \@testlist, $thistag,$ssa, $server, $eva_ext,$seq_ext);
next if check_err($AvgQ3, $logfile);


################################################################
################################################################

sub print_help {
	print "\nHelp Summary for sstest.pl script\n";
	print "Written by Matt Spencer\n";
	print "\nDescription:\n";
	print "This script trains and tests a DN for secondary structure prediction.\n";
	print "\nRequired input:\n";
	print "\t-pred	: The directory contains ss predictions (*.pred, *.prob).\n";
	print "\t-out	  : The directory where dnss formatted predictions are saved.\n";
	print "\t-list	: The proteins list to summarize.\n";
	print "\t-tag	: Identifying tag - the score file will be named as such.\n";
	print "\nOptions:\n";
	print "\t-help	: Print this help message.\n";
	print "\n\n";
}
