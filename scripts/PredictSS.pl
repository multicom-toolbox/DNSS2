#!/usr/bin/perl
# Script:    PredictSS.pl
# Author:    Matt Spencer
# Date Made:    10/8/2013
# Last Mod:    10/29/2014
# Lastest Modified: 11/2/2016 by Jie Hou
my $GLOBAL_PATH;
BEGIN { 
$GLOBAL_PATH='/storage/htc/bdm/jh7x3/DNSS2_github/DNSS2/';
}

use strict;
use lib "$GLOBAL_PATH/lib";
use Time qw(formatted_localtime);
use DN_SSpred2 qw(generate_testfeature_for_convolution_DNSS2 make_testfeatures predict_SS2 timeprint make_dnss_file check_err); 
use SeqAlign qw(load_fasta);
use Getopt::Long;

my ($help, $file, @seq, @name, $indir, $outdir);
my $opt = GetOptions(    '-out:s', \$outdir,
            '-help!', \$help,
            '-seq:s@', \@seq,
            '-name:s@', \@name,
            '-indir:s', \$indir,
            '-file!', \$file );

# Required input is either sequence(s) or an input directory
print_help() if ($help || (!@seq && !$indir));


my $uniclust30db='/storage/htc/bdm/Collaboration/jh7x3/DeepCov_SS_SA_project/DNSS2.0/database/uniclust30_2017_10/uniclust30_2017_10';
my $uniref90db='/storage/htc/bdm/Collaboration/jh7x3/DeepCov_SS_SA_project/DNSS2.0/database/uniref90_20180920/uniref90.fasta';

##################################################
##  Setting up directory architecture
##################################################

$outdir = "." unless ($outdir);
$outdir .= "/" unless ($outdir =~ /\/$/);

# These intermediate directories already exist and are
# designed to remain as are. Files within these temporary
# dirs should be automatically deleted between every 
# prediction, though it shouldn't hurt if they aren't.

my $modeldir = "$GLOBAL_PATH/models/";
my $tempdir = "$outdir/temp/";
my $pssmdir = "$outdir/temp/pssm/";
my $hmmdir = "$outdir/temp/hmm/";
my $probdir = "$outdir/temp/prob/";
my $featdir = "$outdir/temp/feat/";
my $pfeatdir= "$outdir/temp/pfeat/";
my $pprobdir= "$outdir/temp/pprob/";

`mkdir $outdir` unless (-d $outdir);
`mkdir $tempdir` unless (-d $tempdir);
`mkdir $pssmdir` unless (-d $pssmdir);
`mkdir $hmmdir` unless (-d $hmmdir);
`mkdir $probdir` unless (-d $probdir);
`mkdir $featdir` unless (-d $featdir);
`mkdir $pfeatdir` unless (-d $pfeatdir);
`mkdir $pprobdir` unless (-d $pprobdir);

# A log file is generated for each prediction. These logs
# can be found at DNSS/logs/ and are labeled according to
# the time at which this program was called. There is 
# currently no mechanism for deleting old log files.

my $time = formatted_localtime();
$time =~ s/[\[\]]//g;
$time =~ s/\//_/g;
$time =~ s/:/_/g;
$time =~ s/\s/-/;
$time =~ s/ //;
my $logfile = "$GLOBAL_PATH/logs/$time.log";

my $ii=1;
while (-f $logfile){
    $logfile = "$GLOBAL_PATH/logs/$time.$ii.log";
    $ii++;
}

#my @DN1s = ($modeldir . "DN1.model.dat", $modeldir . "DN2.model.dat");
#my @DN1s_util = ($modeldir . "/util_models/DN1.model.util.dat", $modeldir . "/util_models/DN2.model.util.dat");
#my $DN2 = $modeldir . "DN3.model.dat";
#my $DN2_util = $modeldir . "/util_models/DN3.model.util.dat";
my @DN1s = ("deepss_1dconv","deepss_1dCRMN","deepss_1dFrac","deepss_1dInception","deepss_1dRCNN","deepss_1dResnet");

unless (-d $tempdir){
    print "Temp directory was unable to be created. Aborting.\n";
    timeprint($logfile, "Temp directory was unable to be created. Aborting.");
    exit;
}

timeprint($logfile, "Running PredictSS.pl", 2);

##################################################
##  Adjusting for different input options
##################################################

my @fasta;
if ($file) {            # Input specified as individual files
    @fasta = @seq;
    @name = ();
    foreach (@seq) {
        chomp($_);
        my @paths = split('/', $_);
        my @info = split(/\./, $paths[$#paths]);
        my $fullname = join(/./, @info[0..$#info-1]);
        push (@name, $fullname);
    }
}
elsif ($indir){            # Input specified as a directory of files
    $indir .= "/" unless ($indir =~ /\/$/);

    my $list = `ls $indir`;
    my @files = split(/\n/, $list);
    foreach (@files){
        push (@fasta, $indir . $_) if ($_ =~ /fasta$/);
        my @paths = split('/', $_);
        my @info = split(/\./, $paths[$#paths]);
        my $fullname = join(/./, @info[0..$#info-1]);
        push (@name, $fullname);
    }
}
else {                # Input simply a sequence
    for (my $ii=0; $ii<@seq; $ii++){
        my $thisname;
        if ($name[$ii]){
            $thisname = $name[$ii];    
        }
        else {
            $thisname = sprintf("%04s", $ii);
            push (@name, $thisname);
        }
        my $fastfile = $tempdir . "$thisname.fasta";
        push (@fasta, $fastfile);
        open (OUT, ">$fastfile");
        print OUT $seq[$ii];
        close OUT;
    }
}
my $device = 'CPU';
timeprint($logfile, "Calling CPU to run prediction!!!");

##################################################
##  SS Prediction on input files
##################################################

my $percent = 0;
my $probfile;
my $predfile;
for (my $ii=0; $ii<@fasta; $ii++){
    if ($ii/@fasta*100 >= $percent && $percent < 100){
        timeprint($logfile, "$percent % complete!");
        $percent += 10;
    }
    remove_intermediates();        # Temporary files removed between each prediction
    timeprint($logfile, "Predicting $name[$ii]...");

    unless (-f $fasta[$ii]){
        check_err("err: Fasta file not found: $fasta[$ii]", $logfile);
    }
    my $feat = $featdir . "$name[$ii].feat";
    ## Generate PSSM file
    my $pssm = $pssmdir . "$name[$ii].pssm";
    my $log = $pssmdir . "$name[$ii].log";
    if (-f $feat){
        goto SKIP;
    }
    
    timeprint($logfile, "Generating PSSM...\n$GLOBAL_PATH/scripts/generate-pssm_withDB.sh $fasta[$ii] $pssm $log $uniref90db\n");
    
    if(-f $pssm)
    {
      print "PSSM generation done!\n";
    }else{
      `$GLOBAL_PATH/scripts/generate-pssm_withDB.sh $fasta[$ii] $pssm $log $uniref90db`;
    }
    unless (-f $pssm){
        print "PSSM generation failed; using less stringent parameters.\n";
        `$GLOBAL_PATH/scripts/gen-pssm-less-stringent_withDB.sh $fasta[$ii] $pssm $log $uniref90db`;
        unless (-f $pssm){
            check_err("err: PSSM generation failed.", $logfile);
            next;
        }
    }
    
    ## Generate hhblit alignment
    timeprint($logfile, "Generating HHBLITS...\n$GLOBAL_PATH/scripts/generate-hmm_withDB.sh $fasta[$ii] $hmmdir$name[$ii] $uniclust30db\n");
    my $a3m = $hmmdir . "$name[$ii].a3m";
    my $hmm = $hmmdir . "$name[$ii].hmm";
    my $aln = $hmmdir . "$name[$ii].aln";
    if(-f $hmm)
    {
      print "HMM generation done!\n";
    }else{
      `$GLOBAL_PATH/scripts/generate-hmm_withDB.sh $fasta[$ii] $hmmdir$name[$ii] $uniclust30db`;
    }
    system("egrep -v \"^>\" $a3m | sed 's/[a-z]//g' > $aln");
    


    ## get alignment probability
    timeprint($logfile, "Generating Alignment Probability...\nperl $GLOBAL_PATH/scripts/msa2profile_sspro.pl $fasta[$ii] $aln  $hmmdir $GLOBAL_PATH/scripts/msa2profile_sspro $GLOBAL_PATH/scripts/sspro.model\n");
    `perl $GLOBAL_PATH/scripts/msa2profile_sspro.pl $fasta[$ii] $aln  $hmmdir $GLOBAL_PATH/scripts/msa2profile_sspro $GLOBAL_PATH/scripts/sspro.model`;


    
    ## get emission prob and transition prob
    timeprint($logfile, "Generating HMM Probability...\nperl $GLOBAL_PATH/scripts/gen_transition_prob_fromhmm.pl  $fasta[$ii] $hmm\n");
    `perl $GLOBAL_PATH/scripts/gen_transition_prob_fromhmm.pl  $fasta[$ii] $hmm`;


    SKIP:

    ## Generate Feature file
    $feat = $featdir . "$name[$ii].feat";


    #####################################
    my $para_iters = 1;
    my $para_pssm = 1;
    my $para_pssm_uniref50 = 0;
    my $para_hmmEm = 1;
    my $para_hmmTr = 1;
    my $para_aln_sspro = 0;
    my $para_aln_sspro_unif90 = 0;
    my $para_aln_sspro_unif50 = 0;
    my $para_aln_hhblits = 1;
    my $para_atch = 1;
    my $para_seq = 0;
    my $para_wind = 1;
    my $para_boost = 0;
    my $para_skipx = 0;
    my $para_bound = 0;
    my $para_reduced = 0;
    my $para_target = 3*$para_boost;
    my $para_arch = "X,X,$para_target";
    my @archinfo = split(/,/, $para_arch);
    #my $lwind = 17;
    #my $larch = "600,600,600,250,3";
    my $para_ssadir = "X";
    my $para_pssmdir = "$pssmdir";
    my $para_pssm_uniref50dir = "$pssmdir";
    my $para_hmmEmdir = "$hmmdir";
    my $para_hmmTrdir = "$hmmdir";
    my $para_msa2profile_ssprodir = "$hmmdir";
    my $para_msa2profile_sspro_uniref90_dir = "$hmmdir";
    my $para_msa2profile_sspro_uniref50_dir = "$hmmdir";
    my $para_msa2profile_hhblitsdir = "$hmmdir";
    my $overdir = $featdir;
    my $logfile = $overdir . "log.txt";

    my @para_features = ($para_pssm,$para_pssm_uniref50, $para_hmmEm, $para_hmmTr, $para_aln_sspro, $para_aln_sspro_unif50, $para_aln_sspro_unif90, $para_aln_hhblits, $para_atch, $para_seq, $para_bound, $para_skipx);
    my @para_options = (\@para_features, $para_wind, $para_boost, $para_reduced);
    my @para_dirs = ($para_ssadir, $para_pssmdir, $para_pssm_uniref50dir, $para_hmmEmdir, $para_hmmTrdir, $para_msa2profile_ssprodir, $para_msa2profile_sspro_uniref50_dir, $para_msa2profile_sspro_uniref90_dir, $para_msa2profile_hhblitsdir);
    
    `mkdir $overdir` unless (-d $overdir);
    
    timeprint($logfile, "Performing SS train & test using the parameters:\nOutdir: $outdir\nIterations: $para_iters\nReduced: $para_reduced\nSkip X: $para_skipx\n\nWindow: $para_wind\nPssm: $para_pssm\nPssm_uniref50: $para_pssm_uniref50\nHMMEM: $para_hmmEm\nHMMTR: $para_hmmTr\naln_sspro: $para_aln_sspro\naln_sspro_unif50: $para_aln_sspro_unif50\naln_sspro_unif90: $para_aln_sspro_unif90\naln_hhblits: $para_aln_hhblits\nAtch: $para_atch\nSeq: $para_seq\nArch: $para_arch\nBoost: $para_boost\nBound: $para_bound\n");
    
    my @testlist = ($name[$ii]);
    
    generate_testfeature_for_convolution_DNSS2(	$logfile, $overdir,  \@testlist,	\@para_dirs, \@para_options);
    

    ## Predict SS using 1st layer DNs
    my @probfiles;
    my @deletefiles;
    my $err=0;
    for (my $jj=0; $jj<@DN1s; $jj++){
            my $select_model = $DN1s[$jj];
            timeprint ($logfile, "Accessing model $select_model");
            
            
            my $model_in="$modeldir/model-train-$select_model.json"; 
            my $model_weight_in= "$modeldir/model-train-weight-$select_model-best-val.h5";
            
            my ($file1, $file2) = predict_SS2( $model_in, $model_weight_in, $feat, $probdir, $select_model);
            $err++ if (check_err($file1, $logfile));
            push (@probfiles, $file1);
            push (@deletefiles, $file2);
    }
    next if ($err);

    ## ensemble predictions
    timeprint ($logfile, "Ensemble predictions ...");
    
  
    my %prob_avg=();
    my $net_num = 0;
    for (my $jj=0; $jj<@DN1s; $jj++)
    {
      my $select_model = $DN1s[$jj];
      my $probfile = $probdir . $name[$ii].$DN1s[$jj].".prob";
      #print "$prob_file\n";
      if(!(-e $probfile))
      {
        print "Failed to find $probfile\n";
        next;
      }
      $net_num++;
      open(TMP,"$probfile") || die "Failed to find $probfile\n";
      my $c=0;
      while(<TMP>)
      {
        my $li=$_;
        chomp $li;
        $c++;
        my @tmp = split(/\s++/,$li);
        my $h_prob = $tmp[0];
        my $e_prob = $tmp[1];
        my $c_prob = $tmp[2];
        if(exists($prob_avg{$c}))
        {
          my @tmp2 =split(/\s/,$prob_avg{$c});
          $h_prob = $tmp2[0] + $h_prob;
          $e_prob = $tmp2[1] + $e_prob;
          $c_prob = $tmp2[2] + $c_prob;
          $prob_avg{$c} = "$h_prob $e_prob $c_prob";
        }else{
          $prob_avg{$c} = "$h_prob $e_prob $c_prob";
        }
      }
      close TMP;
      
    }
    my $out_file1="$probdir/$name[$ii].Ensemble.prob";
    my $out_file2="$probdir/$name[$ii].Ensemble.pred";
    open(OUT1,">$out_file1") || die "Failed to find $out_file1\n";
    open(OUT2,">$out_file2") || die "Failed to find $out_file2\n";
    foreach  my $indx (sort { $a <=> $b } keys %prob_avg) {
      my @tmp2 =split(/\s/,$prob_avg{$indx});
      my $max_index = maxindex(\@tmp2);
      print OUT1 sprintf("%.6f",$tmp2[0]/$net_num)." ".sprintf("%.6f",$tmp2[1]/$net_num)." ".sprintf("%.6f",$tmp2[2]/$net_num)."\n";
      if($max_index == 0)
      {
        print OUT2 "1 0 0\n";
      }elsif($max_index == 1)
      {
        print OUT2 "0 1 0\n";
      }elsif($max_index == 2)
      {
        print OUT2 "0 0 1\n";
      }else{
        die "wrong ss index $max_index\n";
      }
    }
    close OUT1;
    close OUT2;   
    
    
    my $probfile = "$probdir/$name[$ii].Ensemble.prob";
    ## Write output DNSS file
    my $dnss = $outdir . "$name[$ii].dnss";
    my $vdnss = $outdir . "$name[$ii].vdnss";
    my $header = ">$name[$ii]";
    my @dirarray = ($probfile);
    my $return = make_dnss_file (\@dirarray, $pssm, $dnss, $vdnss, $header);
    next if (check_err($return, $logfile));
    my $dnss2 = "$probdir/$name[$ii].Ensemble.dnss";
    my $vdnss2 = "$probdir/$name[$ii].Ensemble.vdnss";
    my $return = make_dnss_file (\@dirarray, $pssm, $dnss2, $vdnss2, $header);
    next if (check_err($return, $logfile));
    
    
    #### convert individual network to ss predictions 
    for (my $jj=0; $jj<@DN1s; $jj++)
    {
      my $select_model = $DN1s[$jj];
      my $probfile = $probdir . $name[$ii].$DN1s[$jj].".prob";
      #print "$prob_file\n";
      if(!(-e $probfile))
      {
        print "Failed to find $probfile\n";
        next;
      }
      my $dnss = $probdir . $name[$ii].$DN1s[$jj].".dnss";
      my $vdnss = $probdir . $name[$ii].$DN1s[$jj].".vdnss";
      my $header = ">$name[$ii]";
      my @dirarray = ($probfile);
      my $return = make_dnss_file (\@dirarray, $pssm, $dnss, $vdnss, $header);
      next if (check_err($return, $logfile));
      
    }
    
    timeprint ($logfile, "$name[$ii] prediction successful! Saved to $dnss and $vdnss");    
}
timeprint ($logfile, "100% complete!");


# This is used to remove temporary files between each prediction
sub remove_intermediates {
#    `rm -f temp/feat/* temp/pfeat/* temp/pprob/* temp/prob/* temp/pssm/*`;
}


sub print_help{


    print "#################################################################################################################\n";
    print "#                                                                                                               #\n";
    print "#   Software     :  DNSS (A Deep Learning Network Approach to ab initio Protein Secondary Structure Prediction) #\n";
    print "#   Release      :  2.0  (October 2018)                                                                         #\n";
    print "#                                                                                                               #\n";
    print "#   Author(s)    :  Jie Hou, Zhiye Guo, and Jianlin Cheng                                                       #\n";
    print "#   Copyright    :  Bioinformatics, Data Mining, and Machine Learning Lab (BDM)                                 #\n";
    print "#                    Department of Computer Science                                                             #\n";
    print "#                    University of Missouri, Columbia                                                           #\n";
    print "#                                                                                                               #\n";
    print "#################################################################################################################\n";

    print "\nRequired input:\n";
    print "\t-seq    : Sequence of interest (can have multiple inputs, e.g. -seq AAA -seq AAB -seq AAC\n";
    print "\nOutput:\n";
    print "\t.dnss    : File giving the sequence and SS prediction horizontally.\n";
    print "\t.vdnss    : Gives confidence levels for each prediction, additionally.\n";
    print "\nOptions:\n";
    print "\t-name    : Name of sequence (in same order). If files are given, file names will be used.\n";
    print "\t-file    : Indicates that fasta inputs are fasta files instead of sequences.\n";
    print "\t-indir : Indicate a directory full of fastas to predict.\n";
    print "\t-out    : Dictate the location of prediction files (default is . )\n";
    print "\t-help  : Print this help message.\n\n";
    print "* There are two ways to indicate the protein to predict:\n";
    print "----------------------------------------------------------------------------\n";
    print "(1) Predict from protein file:\n";
    print "   \n";
    print "   Usage:\n";
    print "   \$ perl run_DNSS2.pl -seq <file name>.fasta -file -out <output folder>\n";
    print "   \n";
    print "   Example:\n";
    print "   \$ source ~/python_virtualenv_DNSS2/bin/activate\n";
    print "   \$ perl run_DNSS2.pl -seq test/1GNY-A.fasta -file -out output/\n";
    print "          \n";
    print "----------------------------------------------------------------------------\n";
    print "(2) Predicting multiple proteins:\n";
    print "   \n";
    print "   \n";
    print "   Usage:\n";
    print "   \$ perl run_DNSS2.pl -indir <input directory> -out <output directory>\n";
    print "  \n";
    print "   Example:\n";
    print "   \$ source ~/python_virtualenv_DNSS2/bin/activate\n";
    print "   \$ perl run_DNSS2.pl -indir ./test/ -out ./output/\n";
    print "----------------------------------------------------------------------------\n";
    print "\n";
    print " If you have questions, please contact:   \n";     
    print " Jie Hou(jh7x3\@mail.missouri.edu)\n";
    print " Bioinformatics, Data Mining, and Machine Learning Lab (BDM)\n";
    print " Department of Computer Science\n";
    print " University of Missouri, Columbia\n";
    print " Email: chengji\@missouri.edu\n";

    print " ----------------------------------------------------------------------------\n";

    print " Citation: Spencer, Matt, Jesse Eickholt, and Jianlin Cheng. \"A deep learning network approach to ab initio protein secondary structure prediction.\" IEEE/ACM Transactions on Computational Biology and Bioinformatics (TCBB) 12.1 (2015): 103-112.\n";

    print "\n\n";
    exit;
}

sub maxindex {
  my( $aref, $idx_max ) = ( shift, 0 );
  $aref->[$idx_max] > $aref->[$_] or $idx_max = $_ for 1 .. $#{$aref};
  return $idx_max;
}
