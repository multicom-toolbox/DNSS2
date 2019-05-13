package DN_SSpred2;

# Module:	DN_SSpred.pm
# Author:	Matt Spencer,extended by Jie
# Made: 	8/30/13
# Last Mod:	11/08/13
#
# This module contains the functions necessary to use Deep Networks
# to predict Protein Secondary Structure. The functions are separated
# into three categories: the Training Procedure, Score Calculations, 
# and Other.
#
# The Training Procedure contains functions that will turn data into
# features, compile these features into training and testing files,
# train a DN with the data, and use the DN to predict structure.
#
# The Score Calculations allow the assessment of the prediction 
# process, using the Sov score and Q3 score to compare a predicted
# structure to a structure assumed to be native.
#
# Other Functions are more general functions that are used in 
# various places throughout the process.
#
# Functions:
#
#	TRAINING_PROCEDURE 
#	newDN : Makes Training File, test feat files, trains DN, tests DN
#	features : makes (test)feature arrays out of given input
#		make_features, make_testfeatures, all_testfeatures, print_features
#		linearize, linscale, normalize, normalize_array, get_avg
#	train_DN : calls the python training script with specs
#	test_DN : calls the python testing script with specs
#	predict_SS : same as test_DN for single prot predictions
#	all_dnss_files : makes dnss files out of prediction files
#		make_dnss_file
#
#	SCORE_CALCULATIONS 
#	score_probs : calls all_score with predefined file architecture
#	score_dnss : scores a list of prots from dnss files
#	all_score : gets scores of a list of prots. Returns avg and prints
#	sum_probs : combines multiple prob files into array of probs
#		add_probs, probs_to_ss
#	calc_sov : calculates Q and Sov scores from ref and eva seqs
#		len, seg_len, max, min, separate, get_first, get_last,
#		import_struct
#		import_ssa, import_pred, import_dnss, import_ss_sa, import_psi
#		calculate_RQ3, calculate_C, calculate_Q3, calculate_sov
#	calc_confusion_matrix : calculates the confusion matrix of a list
#		confusion_matrix, add_matrix, print_matrix
#	sort_scores : sorts avg scores according to sum of different scores
#		sort_scores_by_rank, bycol
#
#	OTHER_FUNCTIONS
#	load_results : gets scores from saved results.txt file
#	write_results : writes scores to saved results.txt file
#	shuffle : randomizes the order of elements in an array
#	timeprint : prints a message to a log file and to stderr
#	check_err : tests for an error message and probably prints it
#	get_targsize : Returns parameter size of target
my $GLOBAL_PATH;
BEGIN { $GLOBAL_PATH = '/storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/benchmark_tools/DNSS2/'; }
use lib "$GLOBAL_PATH/lib/";
use strict;
use PDBUtils;
use Time qw(formatted_localtime);
use ProteinUtils qw(scaled_log_odds_n_inf atchley_factor);
use SeqAlign qw(load_ssa load_pssm load_seq load_dnss align_and_print remove_gaps get_align_arrays apply_alignment print_comparison print_comparison_named print_separator print_double_separator align_seqs clean_seqs);
use Exporter;
use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(generate_feature_for_convolution generate_feature_for_convolution_DNSS2_exp  generate_testfeature_for_convolution_DNSS2 newDN make_features all_testfeatures make_testfeatures print_features
	train_DN test_DN all_dnss_files predict_SS  predict_SS2 make_dnss_file
	score_probs score_dnss sort_scores all_score
	calc_confusion_matrix print_matrix
	import_ssa import_dnss
	load_results write_results shuffle timeprint check_err);

#############################
#############################
##   TRAINING_PROCEDURE    ##
#############################
#############################

################################################################
# Name: newDN
# In:	logfile : file to write errors and progress to
# 	outdir  : home directory of the DN
# 	trlist_ref : ref of array containing training prot names
#	telist_ref : ref of array containing testing prot names
#	featdirs : ref of array containing directories used for features
#		ssadir, pssmdir, (optional: one to three predirs)
#	featopts : ref of array containing feature generation options
#		pssm, atch, seq
# Out:  predir : directory of predicted proteins from test list
################################################################

sub generate_feature_for_convolution { 
	my ($logfile, $outdir, $trlist_ref, $telist_ref, $featdirs, $featopts) = @_;

	################################################################
	##	Generating Training File
	################################################################

	$outdir .= "/" unless ($outdir =~ m/\/$/);
	my $TFfile;
	my ($ssadir, $pssmdir, @predirs) = @{ $featdirs };

	timeprint($logfile, "Generating Training File...");

	my $count = 0;
	my $percent = 0;

	#Get feature lines from each file, compile into array
	my @TFlines;
	my @trainlist = @{ $trlist_ref };
	foreach (@trainlist){
		if ($count/@trainlist*100 > $percent){
			#print STDERR "$percent % completed (train).\n";
			$percent += 50;
		}
		$count++;

		my $ssafile = $ssadir . "$_.ssa";
		my $pssmfile = $pssmdir . "$_.pssm";
		my @files = ($ssafile, $pssmfile);

		my $prot = $_;
		foreach my $dir (@predirs){
			next unless (-d $dir);
			my $probfile = $dir . "$prot.prob";
			push (@files, $probfile) if (-f $probfile);
		}

		my @features = make_features(\@files, $featopts);
		unless (check_err($features[0], $logfile)){
       $TFfile = $outdir . "$_.fea";
       open (OUT, ">$TFfile") or return "err: newDN: Couldn't open file $TFfile\n";
      	foreach (@features){
      		print OUT $_;
      	}
      	close OUT or return "err: newDN: Couldn't close file $TFfile\n";
		}
	}
	#print STDERR "100% training data completed.\n";

	################################################################
	##	Generating Testing File
	################################################################
	timeprint($logfile, "Generating Testing File...");

	$count = 0;
	$percent = 0;

	#Get feature lines from each file, compile into array
	@TFlines;
	my @testlist = @{ $telist_ref };
	foreach (@testlist){
		if ($count/@testlist*100 > $percent){
			#print STDERR "$percent % completed (test).\n";
			$percent += 50;
		}
		$count++;

		my $ssafile = $ssadir . "$_.ssa";
		my $pssmfile = $pssmdir . "$_.pssm";
		my @files = ($ssafile, $pssmfile);

		my $prot = $_;
		foreach my $dir (@predirs){
			next unless (-d $dir);
			my $probfile = $dir . "$prot.prob";
			push (@files, $probfile) if (-f $probfile);
		}

		my @features = make_features(\@files, $featopts);
		unless (check_err($features[0], $logfile)){
       $TFfile = $outdir . "$_.fea";
       open (OUT, ">$TFfile") or return "err: newDN: Couldn't open file $TFfile\n";
      	foreach (@features){
      		print OUT $_;
      	}
      	close OUT or return "err: newDN: Couldn't close file $TFfile\n";
		}
	}
	#print STDERR "100% testing data completed.\n";
}


sub generate_feature_for_convolution_DNSS2_exp { 
	my ($logfile, $outdir, $trlist_ref, $telist_ref, $featdirs, $featopts) = @_;

	################################################################
	##	Generating Training File
	################################################################

	$outdir .= "/" unless ($outdir =~ m/\/$/);
	my $TFfile;
	my ($ssadir, $pssmdir, $pssm_uniref50dir, $hmmEmdir, $hmmTrdir, $msa2profile_ssprodir, $msa2profile_sspro_uniref50_dir, $msa2profile_sspro_uniref90_dir, $msa2profile_hhblitsdir, @predirs) = @{ $featdirs };

	timeprint($logfile, "Generating Training File...");

	my $count = 0;
	my $percent = 0;

	#Get feature lines from each file, compile into array
	my @TFlines;
	my @trainlist = @{ $trlist_ref };
	foreach (@trainlist){
		if ($count/@trainlist*100 > $percent){
			#print STDERR "$percent % completed (train).\n";
			$percent += 50;
		}
		$count++;

		my $ssafile = $ssadir . "$_.ssa";
		my $pssmfile = $pssmdir . "$_.pssm";
		my $pssm_uniref50file = $pssm_uniref50dir . "$_.pssm";
		my $hmmEmfile = $hmmEmdir . "$_.hmm2prob";
		my $hmmTrfile = $hmmTrdir . "$_.hmm2prob";
		my $msa2profile_ssprofile = $msa2profile_ssprodir . "$_.align2profile";
		my $msa2profile_sspro_uniref50file = $msa2profile_sspro_uniref50_dir . "$_.align2profile";
		my $msa2profile_sspro_uniref90file = $msa2profile_sspro_uniref90_dir . "$_.align2profile";
		my $msa2profile_hhblitfile = $msa2profile_hhblitsdir . "$_.align2profile";
		my @files = ($ssafile, $pssmfile,$pssm_uniref50file,$hmmEmfile,$hmmTrfile,$msa2profile_ssprofile,$msa2profile_sspro_uniref50file,$msa2profile_sspro_uniref90file,$msa2profile_hhblitfile);

		my $prot = $_;
		foreach my $dir (@predirs){
			next unless (-d $dir);
			my $probfile = $dir . "$prot.prob";
			push (@files, $probfile) if (-f $probfile);
		}

		my @features = make_features(\@files, $featopts);
		unless (check_err($features[0], $logfile)){
       $TFfile = $outdir . "$_.fea";
       open (OUT, ">$TFfile") or return "err: newDN: Couldn't open file $TFfile\n";
      	foreach (@features){
      		print OUT $_;
      	}
      	close OUT or return "err: newDN: Couldn't close file $TFfile\n";
		}
	}
	#print STDERR "100% training data completed.\n";

	################################################################
	##	Generating Testing File
	################################################################
	timeprint($logfile, "Generating Testing File...");

	$count = 0;
	$percent = 0;

	#Get feature lines from each file, compile into array
	@TFlines;
	my @testlist = @{ $telist_ref };
	foreach (@testlist){
		if ($count/@testlist*100 > $percent){
			#print STDERR "$percent % completed (test).\n";
			$percent += 50;
		}
		$count++;

		my $ssafile = $ssadir . "$_.ssa";
		my $pssmfile = $pssmdir . "$_.pssm";
		my $pssm_uniref50file = $pssm_uniref50dir . "$_.pssm";
		my $hmmEmfile = $hmmEmdir . "$_.hmm2prob";
		my $hmmTrfile = $hmmTrdir . "$_.hmm2prob";
		my $msa2profile_ssprofile = $msa2profile_ssprodir . "$_.align2profile";
		my $msa2profile_sspro_uniref50file = $msa2profile_sspro_uniref50_dir . "$_.align2profile";
		my $msa2profile_sspro_uniref90file = $msa2profile_sspro_uniref90_dir . "$_.align2profile";
		my $msa2profile_hhblitfile = $msa2profile_hhblitsdir . "$_.align2profile";
		my @files = ($ssafile, $pssmfile,$pssm_uniref50file,$hmmEmfile,$hmmTrfile,$msa2profile_ssprofile,$msa2profile_sspro_uniref50file,$msa2profile_sspro_uniref90file,$msa2profile_hhblitfile);

		my $prot = $_;
		foreach my $dir (@predirs){
			next unless (-d $dir);
			my $probfile = $dir . "$prot.prob";
			push (@files, $probfile) if (-f $probfile);
		}

		my @features = make_features(\@files, $featopts);
		unless (check_err($features[0], $logfile)){
       $TFfile = $outdir . "$_.fea";
       open (OUT, ">$TFfile") or return "err: newDN: Couldn't open file $TFfile\n";
      	foreach (@features){
      		print OUT $_;
      	}
      	close OUT or return "err: newDN: Couldn't close file $TFfile\n";
		}
	}
	#print STDERR "100% testing data completed.\n";
}


sub generate_testfeature_for_convolution_DNSS2 { 
	my ($logfile, $outdir, $telist_ref, $featdirs, $featopts) = @_;

	################################################################
	##	Generating Training File
	################################################################

	$outdir .= "/" unless ($outdir =~ m/\/$/);
	my $TFfile;
	my ($ssadir, $pssmdir, $pssm_uniref50dir, $hmmEmdir, $hmmTrdir, $msa2profile_ssprodir, $msa2profile_sspro_uniref50_dir, $msa2profile_sspro_uniref90_dir, $msa2profile_hhblitsdir, @predirs) = @{ $featdirs };

	################################################################
	##	Generating Testing File
	################################################################
	timeprint($logfile, "Generating Testing File...");

	my $count = 0;
	my $percent = 0;

	#Get feature lines from each file, compile into array
	my @TFlines;
	my @testlist = @{ $telist_ref };
	foreach (@testlist){
		if ($count/@testlist*100 > $percent){
			#print STDERR "$percent % completed (test).\n";
			$percent += 50;
		}
		$count++;

		my $ssafile = $ssadir . "$_.ssa";
		my $pssmfile = $pssmdir . "$_.pssm";
		my $pssm_uniref50file = $pssm_uniref50dir . "$_.pssm";
		my $hmmEmfile = $hmmEmdir . "$_.hmm2prob";
		my $hmmTrfile = $hmmTrdir . "$_.hmm2prob";
		my $msa2profile_ssprofile = $msa2profile_ssprodir . "$_.aln2profile";
		my $msa2profile_sspro_uniref50file = $msa2profile_sspro_uniref50_dir . "$_.aln2profile";
		my $msa2profile_sspro_uniref90file = $msa2profile_sspro_uniref90_dir . "$_.aln2profile";
		my $msa2profile_hhblitfile = $msa2profile_hhblitsdir . "$_.aln2profile";
		my @files = ($ssafile, $pssmfile,$pssm_uniref50file,$hmmEmfile,$hmmTrfile,$msa2profile_ssprofile,$msa2profile_sspro_uniref50file,$msa2profile_sspro_uniref90file,$msa2profile_hhblitfile);

		my $prot = $_;
		foreach my $dir (@predirs){
			next unless (-d $dir);
			my $probfile = $dir . "$prot.prob";
			push (@files, $probfile) if (-f $probfile);
		}

		my @features = make_testfeatures(\@files, $featopts);
		unless (check_err($features[0], $logfile)){
       $TFfile = $outdir . "$_.feat";
       open (OUT, ">$TFfile") or return "err: newDN: Couldn't open file $TFfile\n";
      	foreach (@features){
      		print OUT $_;
      	}
      	close OUT or return "err: newDN: Couldn't close file $TFfile\n";
		}
	}
	#print STDERR "100% testing data completed.\n";
}




sub newDN {
	my ($logfile, $outdir, $trlist_ref, $telist_ref, $featdirs, $featopts, $arch_ref, $tag) = @_;

	################################################################
	##	Generating Training File
	################################################################

	$outdir .= "/" unless ($outdir =~ m/\/$/);
	my $modeldir = $outdir . "models/";
	my $TFfile = $modeldir . "ss_TF.txt";
	my $dnfile = $modeldir . "ss.model.dat";
	my ($ssadir, $pssmdir, @predirs) = @{ $featdirs };
	`mkdir $modeldir` unless (-d $modeldir);

	timeprint($logfile, "Generating Training File...");

	#Skip this section if training file is already made
	if (-f $TFfile || -f $dnfile){
		timeprint ($logfile, "Previously accomplished.");
		goto NEXT1;
	}
	my $count = 0;
	my $percent = 0;

	#Get feature lines from each file, compile into array
	my @TFlines;
	my @trainlist = @{ $trlist_ref };
	foreach (@trainlist){
		if ($count/@trainlist*100 > $percent){
			#print STDERR "$percent % completed.\n";
			$percent += 50;
		}
		$count++;

		my $ssafile = $ssadir . "$_.ssa";
		my $pssmfile = $pssmdir . "$_.pssm";
		my @files = ($ssafile, $pssmfile);

		my $prot = $_;
		foreach my $dir (@predirs){
			next unless (-d $dir);
			my $probfile = $dir . "$prot.prob";
			push (@files, $probfile) if (-f $probfile);
		}

		my @features = make_features(\@files, $featopts);
		unless (check_err($features[0], $logfile)){
			push (@TFlines, @features);
		}
	}

	# Lines are shuffled to prevent order bias. Ex: otherwise the first
	# hundred lines or so would all be from the same protein.
	shuffle (\@TFlines);
	open (OUT, ">$TFfile") or return "err: newDN: Couldn't open file $TFfile\n";
	foreach (@TFlines){
		print OUT $_;
	}
	close OUT or return "err: newDN: Couldn't close file $TFfile\n";

	#print STDERR "100% completed.\n";

	NEXT1:
	################################################################
	##	Making test features
	################################################################

	timeprint($logfile, "Generating Test Features...");

	my $testdir = $outdir . "testfeatures/";
	`mkdir $testdir` unless (-d $testdir);

	# Testfeature files are made from the test protein list
	# Impossible to reliably skip this step here, but fcn checks for 
	# existence of files before regenerating them
	my @errors = all_testfeatures($testdir, $featdirs, $telist_ref, $featopts);
	foreach (@errors){
		print STDERR "$_\n";
		`echo "$_\n" >> $logfile`;
	}

	################################################################
	##	Training DN
	################################################################

	timeprint($logfile, "Training DN...");

	my $dnfile = $modeldir . "ss.model.dat";
	my @arch = @{ $arch_ref };

	# This step is skipped if the model file already exists
	if (-f $dnfile){
		timeprint ($logfile, "Previously accomplished.");
		goto NEXT3;
	}

	my $archstring = "";
	foreach (@arch){
		$archstring .= "$_,";
	}
	$archstring = substr($archstring, 0, -1);

	# DN is trained with Training File data
	my @err = train_DN($modeldir, $TFfile, $archstring);
	if ($err[0] eq "err"){
		print STDERR "$err[1]\n";
		`echo "$err[1]\n" >> $logfile`;
	}

	NEXT3:
	################################################################
	##	Testing DNs
	################################################################

	timeprint($logfile, "Testing DN with test prots...");

	my $predir = $outdir . "pred/";
	`mkdir $predir` unless (-d $predir);
	my $dnssdir = $outdir . "dnss/";
	`mkdir $dnssdir` unless (-d $dnssdir);

	# DN is tested on test proteins.
	# Again, cannot skip here, but individual test can be skipped in fcn
	test_DN ($dnfile, $testdir, $predir, $telist_ref, $arch[$#arch]);

	################################################################
	##	Generating Predictions
	################################################################

	timeprint($logfile, "Generating DNSS prediction files...");

	# Probabilities are reformatted into usable file aligning the sequence
	# with the corresponding predicted most likely secondary structure
	my @dirarray = ($predir);
	all_dnss_files(\@dirarray, $pssmdir, $dnssdir, $telist_ref, $tag);

	return ($predir, $dnssdir);
}

################################################################
# Name: features
# In:	file_ref: array ref including:
#		ssafile: path to ssa file
#		pssmfile: path to pssm file
#		probfile1: (optional) path to prob file 1
#		probfile2: (optional) path to prob file 2
#	opt_ref: array ref including options:
#		feats_ref: array ref including binary values to indicate:
#			pssm : include pssm
#			atch : include atch
#			seq : include seq
#			bound : include boundary indicators
#			skipx : 0: include X
#				1: skip X lines
#				2: skip windows including X
#		window: window size
#		boostwind: boosting window size
#		reduced: make testfeats gapless
# Out:	array containing feature lines, or array of ("err", <err>)
# 
# Subs: make_features: use this public call to make normal features
#	make_testfeatures: use this public call to make test features
#	all_testfeatures: use to make a whole list of testfeat files
#	features: private algorithm to make features
#	print_features: prints given feature lines to given file
#	linearize: converts PSSM file to feature data
#	linscale: converts PSSM value to [0,1]
#	normalize: makes probs sum to 1 (used for layer features)
#	get_avg: finds average of probs predicting the same thing (used for LF)
################################################################

sub all_testfeatures {
	# Note that test feature generation does not require SSA dir input
	# featdir_ref : (($ssadir), $pssmdir, @probdirs)
	my ($outdir, $featdir_ref, $list_ref, $featopts) = @_;

	my $count = 0;
	my $percent = 0;
	my @errlines;
	my ($ssadir, $pssmdir, @probdirs) = @{ $featdir_ref };
	my @list = @{ $list_ref };

	# Generates a test file for each protein in the list
	foreach my $prot (@list){
		if ($count/@list*100 > $percent){
			#print STDERR "$percent % completed.\n";
			$percent += 50;
		}
		$count++;

		my $ssafile = 0; # If not given, input of 0 for this file is appropriate
		$ssafile = $ssadir . "$prot.ssa" if (-d $ssadir);
		my $testfeatfile = $outdir . "$prot.feat";
		my $pssmfile = $pssmdir . "$prot.pssm";
		my @files = ($ssafile, $pssmfile);

		# Skips this file generation if file already exists
		next if (-f $testfeatfile);

		foreach my $dir (@probdirs){
			next unless (-d $dir);
			my $probfile = $dir . "$prot.prob";
			push (@files, $probfile) if (-f $probfile);
		}

		# Generates feature lines for this protein
		my @features = make_testfeatures(\@files, $featopts);
		if (check_err($features[0])){ push (@errlines, $features[0]); next; }

		# Prints feature lines to designated file location
		my $err = print_features($testfeatfile, @features);
		push (@errlines, $err) if (check_err($err));
	}

	#print STDERR "100% completed.\n";
	# Returns error messages that were encountered during process
	return @errlines;
}

sub make_features {
	push (@_, 0);
	return features(@_);
}

sub make_testfeatures {
	push (@_, 1); # Safeguard in case want test but input ssa file anyway
	return features(@_);
}

sub features { 
	my ($file_ref, $opt_ref, $test) = @_;
	my ($ssafile, $pssmfile,$pssm_uniref50file,$hmmEmfile,$hmmTrfile,$msa2profile_ssprofile,$msa2profile_sspro_uniref50file,$msa2profile_sspro_uniref90file,$msa2profile_hhblitfile, @probfiles) = @{ $file_ref };
	my ($feats_ref, $window, $boost_wind, $reduced) = @{ $opt_ref  };
	my ($pssm,$pssm_uniref50, $hmmEm, $hmmTr, $aln_sspro, $aln_sspro_unif50, $aln_sspro_unif90, $aln_hhblits, $atch, $seq, $bound, $skipx) = @{ $feats_ref };
	# SSA file not necessary for making testfeatures, only features
	# Prob files never necessary
	return ("err: features: Couldn't find file $ssafile") unless ((-f $ssafile) || ($test && !$reduced));
	return ("err: features: Couldn't find file $pssmfile") unless (-f $pssmfile);
	return ("err: features: Couldn't find file $pssmfile") unless (-f $pssm_uniref50file);
	foreach (@probfiles){
		return ("err: features: Couldn't find file $_") unless (-f $_);
	}

	# Ensures that window size is odd. Default value is 15
	if ($window > 0){
		$window++ if ($window % 2 == 0);
	} else {
		$window = 15;
	}

	# Ensures that boost window size is odd. Default value is 1
	if ($boost_wind > 0){ 
		$boost_wind++ if ($boost_wind % 2 == 0); 
	}
	else { $boost_wind = 1; }

	#################################
	# Initialize Seq information
	#################################

	my @AAs = qw/A C D E F G H I K L M N P Q R S T V W Y/;
	my %SSHotCode;

	@{ $SSHotCode{"H"} }  = qw/0 0 1/;
	@{ $SSHotCode{"E"} }  = qw/0 1 0/;
	@{ $SSHotCode{"C"} }  = qw/1 0 0/;
	@{ $SSHotCode{"-"} }  = qw/0 0 0/; #only used if boosting

	# Makes binary code representing each residue
	my %AAHotCode;
	for (my $ii=0; $ii<20; $ii++){
		my @hotcode;
		for (my $jj=0; $jj<21; $jj++){
			if ($ii == $jj){
				push(@hotcode, 1);
			}
			else {
				push(@hotcode, 0);
			}
		}
		@{ $AAHotCode{$AAs[$ii]} } = @hotcode;
	}
	@{ $AAHotCode{"-"} } =  qw/0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1/;
	@{ $AAHotCode{"X"} } =  qw/0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1/; 


  #################################
	# Collect PSSM information
	#################################

	open (PSSM, $pssmfile) or return ("err: features: Couldn't open pssm file");
	my @lines = <PSSM>;
	my $pfilelen = @lines;
	chomp(@lines);
	close PSSM;

	my @pssmlines;
	@pssmlines = load_pssm($pssmfile);
  if($pssm_uniref50)
  {
  	#################################
  	# Collect PSSM information
  	#################################
  
  	open (PSSM, $pssm_uniref50file) or return ("err: features: Couldn't open pssm file");
  	@lines = <PSSM>;
  	$pfilelen = @lines;
  	chomp(@lines);
  	close PSSM;
  
  	@pssmlines;
  	@pssmlines = load_pssm($pssm_uniref50file);
  }
	#################################
	# Collect HMM information
	#################################
  #print "Loading $hmmEmfile\n";
  open (HMM, $hmmEmfile) or return ("err: features: Couldn't open hmm file $hmmEmfile");  
  my @hmmlines = <HMM>;
  close HMM;
  
  shift @hmmlines;
  
  my @hmmseq;
  my @hmmEm_prob;
  my @hmmTr_prob;
  my $i;
  my $line;
  foreach $line (@hmmlines)
  {
    chomp $line;
    my @tmp = split(/\t/,$line);
    my ($res) = $tmp[0];
    push (@hmmseq,$res);
    
    #### Emission prob
    my $hmmEm_prob_tmp = $tmp[2];
    for($i=3;$i<22;$i++)
    {
     #print $tmp[$i]." ";
      $hmmEm_prob_tmp .= " ".$tmp[$i];
    }
    #print "\n";
    push (@hmmEm_prob,$hmmEm_prob_tmp);
    
    #### Transition prob
    my $hmmTr_prob_tmp = $tmp[22];
    for($i=23;$i<29;$i++)
    {
      $hmmTr_prob_tmp .= " ".$tmp[$i];
    }
    push (@hmmTr_prob,$hmmTr_prob_tmp);
  }

	#################################
	# Collect MSA2PROFILE information
	#################################
  open (PROFILE, $msa2profile_ssprofile) or die ("err: features: Couldn't open $msa2profile_ssprofile file");
  my @msa2profilelines = <PROFILE>;
  close PROFILE;
  
  shift @msa2profilelines; 
  
  my @msa2profileseq_sspro;
  my @msa2profile_sspro_prob;
  my $i;
  foreach $line (@msa2profilelines)
  {
    chomp $line;
    my @tmp = split(/\t/,$line);
    my ($res) = $tmp[1];
    push (@msa2profileseq_sspro,$res);
    
    #### Emission prob
    my $msa2profile_prob_tmp = $tmp[2];
    for($i=3;$i<22;$i++)
    {
     #print $tmp[$i]." ";
      $msa2profile_prob_tmp .= " ".$tmp[$i];
    }
    #print "\n";
    push (@msa2profile_sspro_prob,$msa2profile_prob_tmp);
    
  }
	
	#################################
	# Collect MSA2PROFILE information
	#################################
  open (PROFILE, $msa2profile_sspro_uniref50file) or die ("err: features: Couldn't open $msa2profile_sspro_uniref50file file");
  my @msa2profilelines_uniref50 = <PROFILE>;
  close PROFILE;
  
  shift @msa2profilelines_uniref50;
  
  my @msa2profileseq_sspro_uniref50;
  my @msa2profile_sspro_prob_uniref50;
  my $i;
  foreach $line (@msa2profilelines_uniref50)
  {
    chomp $line;
    my @tmp = split(/\t/,$line);
    my ($res) = $tmp[1];
    push (@msa2profileseq_sspro_uniref50,$res);
    
    #### Emission prob
    my $msa2profile_prob_uniref50_tmp = $tmp[2];
    for($i=3;$i<22;$i++)
    {
     #print $tmp[$i]." ";
      $msa2profile_prob_uniref50_tmp .= " ".$tmp[$i];
    }
    #print "\n";
    push (@msa2profile_sspro_prob_uniref50,$msa2profile_prob_uniref50_tmp);
    
  }
		
	#################################
	# Collect MSA2PROFILE information
	#################################
  open (PROFILE, $msa2profile_sspro_uniref90file) or die ("err: features: Couldn't open $msa2profile_sspro_uniref90file file");
  my @msa2profilelines_uniref90 = <PROFILE>;
  close PROFILE;
  
  shift @msa2profilelines_uniref90;
  
  my @msa2profileseq_sspro_uniref90;
  my @msa2profile_sspro_prob_uniref90;
  my $i;
  foreach $line (@msa2profilelines_uniref90)
  {
    chomp $line;
    my @tmp = split(/\t/,$line);
    my ($res) = $tmp[1];
    push (@msa2profileseq_sspro_uniref90,$res);
    
    #### Emission prob
    my $msa2profile_prob_uniref90_tmp = $tmp[2];
    for($i=3;$i<22;$i++)
    {
     #print $tmp[$i]." ";
      $msa2profile_prob_uniref90_tmp .= " ".$tmp[$i];
    }
    #print "\n";
    push (@msa2profile_sspro_prob_uniref90,$msa2profile_prob_uniref90_tmp);
    
  }
	
	#################################
	# Collect MSA2PROFILE information
	#################################
  open (PROFILE, $msa2profile_hhblitfile) or die ("err: features: Couldn't open msa2profile_hhblitfile file");
  my @msa2profilelines = <PROFILE>;
  close PROFILE;
  
  shift @msa2profilelines;
  
  my @msa2profileseq_hhblit;
  my @msa2profile_hhblit_prob;
  my $i;
  foreach $line (@msa2profilelines)
  {
    chomp $line;
    my @tmp = split(/\t/,$line);
    my ($res) = $tmp[1];
    push (@msa2profileseq_hhblit,$res);
    
    #### Emission prob
    my $msa2profile_prob_tmp = $tmp[2];
    for($i=3;$i<22;$i++)
    {
     #print $tmp[$i]." ";
      $msa2profile_prob_tmp .= " ".$tmp[$i];
    }
    #print "\n";
    push (@msa2profile_hhblit_prob,$msa2profile_prob_tmp);
    
  }

	# The following decides what method will be used to transform the
	# information to the range [0, 1]. Other methods could be better.
#	my ($odds_ref, $inf_ref) = scaled_log_odds_n_inf(@lines);
	my ($odds_ref, $inf_ref) = linearize(@lines);

	my @odds = @{$odds_ref};
	my @inf = @{$inf_ref};

	for (my $ii=0; $ii<@inf; $ii++){
		$inf[$ii] = abs($inf[$ii]);
	}

	my @blankpssm;
	foreach (1..20){
		push (@blankpssm, "0.0");
	}
	push (@blankpssm, "0.0000");

	my @blankhmmEm;
	foreach (1..20){
		push (@blankhmmEm, "0.0");
	}
 
	my @blankhmmTr;
	foreach (1..7){
		push (@blankhmmTr, "0.0");
	}
 
	my @blankmsa2profile_sspro;
	foreach (1..20){
		push (@blankmsa2profile_sspro, "0.0");
	}
 
	my @blankmsa2profile_sspro_uniref50;
	foreach (1..20){
		push (@blankmsa2profile_sspro_uniref50, "0.0");
	}
 
	my @blankmsa2profile_sspro_uniref90;
	foreach (1..20){
		push (@blankmsa2profile_sspro_uniref90, "0.0");
	}
 
	my @blankmsa2profile_hhblit;
	foreach (1..20){
		push (@blankmsa2profile_hhblit, "0.0");
	}
 
	#################################
	# Calculate scalars and acthley factors
	#################################

	my $aas = "ACDEFGHIKLMNPQRSTVWY";
	my @scalar = ();
	my @offset = ();
	my @res = split(//, $aas);
		
	# In this loop, scalars and offsets are calculated for each of the 5
	# Atchley factors, which are later used to transform the factors to be
	# in the range [0, 1], required for DN training
	for(my $ii = 1; $ii <= 5; $ii++) {
		my $min = 15;
		my $max = - 15;
		for(my $jj = 0; $jj < @res; $jj++) {
		
			if(atchley_factor($res[$jj],$ii) < $min) {
				$min = atchley_factor($res[$jj], $ii);
			}
			if(atchley_factor($res[$jj],$ii) > $max) {
				$max = atchley_factor($res[$jj], $ii);
			}

		}

		$scalar[$ii] = (1/($max-$min));
		$offset[$ii] = ($min/($min-$max));
	}

	my %factors = ();

	# 5 factors for all residues are calculated and transformed by the above
	# scalar and offset before being saved to the hash.
	for (my $ii = 0; $ii < @res; $ii++) {
		my @thisatch;
		for(my $jj = 1; $jj <= 5; $jj++) {
			my $atch = atchley_factor($res[$ii], $jj) * $scalar[$jj] + $offset[$jj];
			$atch = sprintf("%.4f", $atch);
			push (@thisatch, $atch);
		}
		$factors{$res[$ii]} =  \@thisatch;
	}

	my @blankatch = (0, 0, 0, 0, 0);
	foreach (@blankatch){
		$_ = sprintf("%.4f", $_);
	}
	$factors{"-"} = \@blankatch;
	$factors{"X"} = \@blankatch;

	#################################
	# Get prob information, if given
	#################################

	# Prob files contain previously predicted probabilities for each secondary
	# structure for this protein. Due to the layered design of the prediction
	# process, the initial prediction will be performed, and the resulting prob
	# files will be used in a refinement prediction.
	my @probs;
	my @fileinfo;

	# Loop is deliberately implemented to handle an unknown number of prob files
	foreach my $file (@probfiles){
		open (PROB, $file) or return ("err: features: Coudln't open file $file.");

		my @lines;
		while (my $probline = <PROB>){
			# Normalization occurs to ensure probabilities sum to 1
			my @nums = normalize($probline);
			push (@lines, \@nums);
		}
		close PROB or return ("err: features: Coudln't close file $file.");
		push (@fileinfo, \@lines);		
	}

	# When multiple prob files are given, the average probability of each 
	# structure at each residue is calculated, giving equal weight to all
	# given prob files (a simple way to adjust weight would be to input
	# prob files for higher priority predictions twice)
	@probs = get_avg($pfilelen, \@fileinfo);
	check_probs("err", @probs);

	my @blankprob;
	for (my $ii=0; $ii<get_targsize(); $ii++){
		 push (@blankprob, "0.000000");
	}

	#################################
	# Collect SSA information and align
	#################################

	# SSA files contain an abbreviated version of the known data of a prot
	# Thus, ssa files are only used for training, since the testing procedure
	# predicts as if there was no known information, and prediction of new
	# proteins will not have the known information.
	my @ssalines;
	if ($ssafile && (!$test || $reduced)){
		open (SSA, $ssafile) or return ("err: features: Couldn't open ssa file");
		@ssalines = <SSA>;
		chomp(@ssalines);
		shift(@ssalines);
		close SSA;

		# Adds lines to SSA info when pssm seq is longer than ssa
		for (my $ii=@ssalines; $ii<@pssmlines; $ii++){
			my $index = $ii+1;
			my $line = "$index\t-\t-\t-";
			push (@ssalines, $line);
		}

		# SSA and PSSM seqs are aligned to be sure information is applied to
		# the correct feature lines
		my @lines = (\@pssmlines, \@ssalines);
		my @ssa_sequence = load_ssa($ssafile);
		my (@align_refs) = get_align_arrays(\@pssmlines, \@ssa_sequence);
		my ($new_pssm_ref, $new_ssa_ref) = apply_alignment(\@lines, \@align_refs);
		@pssmlines = @{ $new_pssm_ref };
		@ssalines = @{ $new_ssa_ref };

		for (my $ii=0; $ii<@ssalines; $ii++){
			$ssalines[$ii] = "-\t-\t-\t-" if ($ssalines[$ii] eq "-");
		}

		# PSSM data is adjusted according to above sequence alignment
		my @data = (\@odds, \@ssalines);
		my ($new_odds, $junk) = apply_alignment(\@data, \@align_refs);
		@odds = @{ $new_odds };

		my @data = (\@inf, \@ssalines);
		my ($new_inf, $junk) = apply_alignment(\@data, \@align_refs);
		@inf = @{ $new_inf };


    # hmmEm data is adjusted according to above sequence alignment
    # hmmTr data is adjusted according to above sequence alignment
    # aln_sspro data is adjusted according to above sequence alignment
    open (SSA, $ssafile) or return ("err: features: Couldn't open ssa file");
		my @ssalines2 = <SSA>;
		chomp(@ssalines2);
		shift(@ssalines2);
		close SSA;

		# Adds lines to SSA info when pssm seq is longer than ssa
		for (my $ii=@ssalines2; $ii<@hmmseq; $ii++){
			my $index = $ii+1;
			my $line = "$index\t-\t-\t-";
			push (@ssalines2, $line);
		}
    my @lines2 = (\@hmmseq, \@ssalines2);
    my @ssa_sequence2 = load_ssa($ssafile);
    my (@align_refs2) = get_align_arrays(\@hmmseq, \@ssa_sequence2);
    my ($new_hmmEm_ref2, $new_ssa_ref2) = apply_alignment(\@lines2, \@align_refs2);
    @hmmseq = @{ $new_hmmEm_ref2 };
    my @ssalines2 = @{ $new_ssa_ref2 };
    
    my @data2 = (\@hmmEm_prob, \@ssalines2);
    my ($new_hmmEm_prob2, $junk2) = apply_alignment(\@data2, \@align_refs2);
    @hmmEm_prob = @{ $new_hmmEm_prob2 };
    
    # hmmTr data is adjusted according to above sequence alignment
    my @data2 = (\@hmmTr_prob, \@ssalines2);
    my ($new_hmmTr_prob2, $junk2) = apply_alignment(\@data2, \@align_refs2);
    @hmmTr_prob = @{ $new_hmmTr_prob2 };

     
    # aln_sspro data is adjusted according to above sequence alignment
    open (SSA, $ssafile) or return ("err: features: Couldn't open ssa file");
		my @ssalines3 = <SSA>;
		chomp(@ssalines3);
		shift(@ssalines3);
		close SSA;

		# Adds lines to SSA info when pssm seq is longer than ssa
		for (my $ii=@ssalines3; $ii<@msa2profileseq_sspro; $ii++){
			my $index = $ii+1;
			my $line = "$index\t-\t-\t-";
			push (@ssalines3, $line);
		}
    my @lines3 = (\@msa2profileseq_sspro, \@ssalines3);
    my @ssa_sequence3 = load_ssa($ssafile);
    my (@align_refs3) = get_align_arrays(\@msa2profileseq_sspro, \@ssa_sequence3);
    my ($new_msa2profileseq_sspro_ref, $new_ssa_ref3) = apply_alignment(\@lines3, \@align_refs3);
    @msa2profileseq_sspro= @{ $new_msa2profileseq_sspro_ref };
    @ssalines3 = @{ $new_ssa_ref3 };
    
    my @data3 = (\@msa2profile_sspro_prob, \@ssalines3);
    my ($new_msa2profilesspro_prob, $junk3) = apply_alignment(\@data3, \@align_refs3);
    @msa2profile_sspro_prob = @{ $new_msa2profilesspro_prob };

    # aln_sspro data is adjusted according to above sequence alignment
    open (SSA, $ssafile) or return ("err: features: Couldn't open ssa file");
		my @ssalines4 = <SSA>;
		chomp(@ssalines4);
		shift(@ssalines4);
		close SSA;

		# Adds lines to SSA info when pssm seq is longer than ssa
		for (my $ii=@ssalines4; $ii<@msa2profileseq_sspro_uniref50; $ii++){
			my $index = $ii+1;
			my $line = "$index\t-\t-\t-";
			push (@ssalines4, $line);
		}
    my @lines4 = (\@msa2profileseq_sspro_uniref50, \@ssalines4);
    my @ssa_sequence4 = load_ssa($ssafile);
    my (@align_refs4) = get_align_arrays(\@msa2profileseq_sspro_uniref50, \@ssa_sequence4);
    my ($new_msa2profileseq_sspro_uniref50_ref, $new_ssa_ref4) = apply_alignment(\@lines4, \@align_refs4);
    @msa2profileseq_sspro_uniref50= @{ $new_msa2profileseq_sspro_uniref50_ref };
    @ssalines4 = @{ $new_ssa_ref4 };
    
    my @data4 = (\@msa2profile_sspro_prob_uniref50, \@ssalines4);
    my ($new_msa2profilesspro_prob_uniref50, $junk4) = apply_alignment(\@data4, \@align_refs4);
    @msa2profile_sspro_prob_uniref50 = @{ $new_msa2profilesspro_prob_uniref50 };
 
 
    # aln_sspro data is adjusted according to above sequence alignment
    open (SSA, $ssafile) or return ("err: features: Couldn't open ssa file");
		my @ssalines5 = <SSA>;
		chomp(@ssalines5);
		shift(@ssalines5);
		close SSA;

		# Adds lines to SSA info when pssm seq is longer than ssa
		for (my $ii=@ssalines5; $ii<@msa2profileseq_sspro_uniref90; $ii++){
			my $index = $ii+1;
			my $line = "$index\t-\t-\t-";
			push (@ssalines5, $line);
		}
    my @lines5 = (\@msa2profileseq_sspro_uniref90, \@ssalines5);
    my @ssa_sequence5 = load_ssa($ssafile);
    my (@align_refs5) = get_align_arrays(\@msa2profileseq_sspro_uniref90, \@ssa_sequence5);
    my ($new_msa2profileseq_sspro_uniref90_ref, $new_ssa_ref5) = apply_alignment(\@lines5, \@align_refs5);
    @msa2profileseq_sspro_uniref90= @{ $new_msa2profileseq_sspro_uniref90_ref };
    @ssalines5 = @{ $new_ssa_ref5 };
    
    my @data5 = (\@msa2profile_sspro_prob_uniref90, \@ssalines5);
    my ($new_msa2profilesspro_prob_uniref90, $junk5) = apply_alignment(\@data5, \@align_refs5);
    @msa2profile_sspro_prob_uniref90 = @{ $new_msa2profilesspro_prob_uniref90 };
 
 
       
    # aln_hhblit data is adjusted according to above sequence alignment
    open (SSA, $ssafile) or return ("err: features: Couldn't open ssa file");
		my @ssalines4 = <SSA>;
		chomp(@ssalines4);
		shift(@ssalines4);
		close SSA;

		# Adds lines to SSA info when pssm seq is longer than ssa
		for (my $ii=@ssalines4; $ii<@msa2profileseq_hhblit; $ii++){
			my $index = $ii+1;
			my $line = "$index\t-\t-\t-";
			push (@ssalines4, $line);
		}
    my @lines4 = (\@msa2profileseq_hhblit, \@ssalines4);
    my @ssa_sequence4 = load_ssa($ssafile);
    my (@align_refs4) = get_align_arrays(\@msa2profileseq_hhblit, \@ssa_sequence4);
    my ($new_msa4profileseq_hhblit_ref, $new_ssa_ref4) = apply_alignment(\@lines4, \@align_refs4);
    @msa2profileseq_hhblit= @{ $new_msa4profileseq_hhblit_ref };
    @ssalines4 = @{ $new_ssa_ref4 };
    
    my @data4 = (\@msa2profile_hhblit_prob, \@ssalines4);
    my ($new_msa2profilehhblit_prob, $junk4) = apply_alignment(\@data4, \@align_refs4);
    @msa2profile_hhblit_prob = @{ $new_msa2profilehhblit_prob };
    
         
		# Prior prediction data is adjusted according to seq alignment
		for (my $ii=@probs; $ii<@pssmlines; $ii++){
			$probs[$ii] = \@blankprob;
		}

		my @info = (\@probs, \@ssalines);
		my ($prob_ref, $junk) = apply_alignment(\@info, \@align_refs);
		@probs = @{ $prob_ref };
	}

	#################################
	# Loop to generate features
	#################################

	my @TFlines;
	my $halfwindow = ($window-1)/2;
	my $hBW = ($boost_wind-1)/2;
	my $min = 0;
	my $max = $#pssmlines;
  
  
	for (my $jj=0; $jj<@pssmlines; $jj++){

		# While making normal features, skip residues that have no native
		# structure, as these have no information to help train the DN.
		if (!$test || $reduced){
			my $ssainfo = $ssalines[$jj];
			my ($ssnum, $ssaa, $thisss, $thisrsa) = split(/\t/, $ssainfo);
			next if ($ssaa eq "-");
		}

		# Lines with X as target residue skipped if specified
		next if ($pssmlines[$jj] eq "X" && $skipx > 0);
		
		# Adds window residue features and builds PSSM feature array
		# Each of these arrays is built out of a single feature type
		# over the span of the whole window. Later, the arrays will be
		# concatenated together into a single string feature line
		my @printarray;
		my @targetfeat;
		my @boundfeat; # 1 if out-of-bounds
		my @seqfeat;
		my @pssmfeat;
		my @hmmEmfeat;
		my @hmmTrfeat;
		my @atchley;
		my @prob;
		my @hmmEmfeat;
		my @hmmTrfeat;
		my @msa2profile_ssprofeat;
		my @msa2profile_sspro_uniref50feat;
		my @msa2profile_sspro_uniref90feat;
		my @msa2profile_hhblitfeat;
	
		# For each line of the PSSM, loop indices spanning the specified window
		for (my $ii=$jj-$halfwindow; $ii<=$jj+$halfwindow; $ii++){

			# Generic information is printed when the window spans past 
			# the terminal residues on either side, or if the structure is
			# known to break the residue sequence of the protein in window
			if ($ii<$min || $ii>$max || $pssmlines[$ii] eq "-"){ 
				push (@targetfeat, @{$SSHotCode{"-"}}) if (abs($ii-$jj)<=$hBW);
				push (@boundfeat, 1);
				push (@seqfeat, @{$AAHotCode{"-"}});
				push (@pssmfeat, @blankpssm);
				push (@hmmEmfeat, @blankhmmEm);
				push (@hmmTrfeat, @blankhmmTr);
				push (@msa2profile_ssprofeat, @blankmsa2profile_sspro);
				push (@msa2profile_sspro_uniref50feat, @blankmsa2profile_sspro_uniref50);
				push (@msa2profile_sspro_uniref90feat, @blankmsa2profile_sspro_uniref90);
				push (@msa2profile_hhblitfeat, @blankmsa2profile_hhblit);
				push (@atchley, @blankatch);
				push (@prob, @blankprob);
				next;
			}

			# Residue value of X is given to all unexpected residue names
			# Windows including X are skipped if specified
			my $aa = $pssmlines[$ii];
			$aa = "X" unless ($aas =~ m/$aa/);
			goto SKIP if ($aa eq "X" && $skipx == 2);

#push (@targetfeat, $jj, $aa) if (abs($ii-$jj)<=$hBW);
		
			# Only check for SSA information if not testfeatures,
			# because SSA file is not required for test
			unless ($test) {
				my $ssainfo = $ssalines[$ii];
				my ($ssnum, $ssaa, $thisss, $thisrsa) = split(/\t/, $ssainfo);
				push (@targetfeat, @{$SSHotCode{$thisss}}) if (abs($ii-$jj)<=$hBW);
			}

			if (!$odds[$ii]){
				my $message = "err: features: Unexpected PSSM format at line $ii"; 
				return $message;
			}
      
      my @hmmEm_prob_arr = split(' ', $hmmEm_prob[$ii]);
      my @hmmTr_prob_arr = split(' ', $hmmTr_prob[$ii]);
      my @msa2profile_sspro_arr = split(' ', $msa2profile_sspro_prob[$ii]);
      my @msa2profile_sspro_uniref50_arr = split(' ', $msa2profile_sspro_prob_uniref50[$ii]);
      my @msa2profile_sspro_uniref90_arr = split(' ', $msa2profile_sspro_prob_uniref90[$ii]);
      my @msa2profile_hhblit_arr = split(' ', $msa2profile_hhblit_prob[$ii]);
      
            
			push (@seqfeat, @{$AAHotCode{$aa}});
			push (@boundfeat, 0);
			push (@pssmfeat, @{$odds[$ii]}, $inf[$ii]);
			push (@hmmEmfeat, @hmmEm_prob_arr);
			push (@hmmTrfeat, @hmmTr_prob_arr);
			push (@msa2profile_ssprofeat, @msa2profile_sspro_arr);
			push (@msa2profile_sspro_uniref50feat, @msa2profile_sspro_uniref50_arr);
			push (@msa2profile_sspro_uniref90feat, @msa2profile_sspro_uniref90_arr);
			push (@msa2profile_hhblitfeat, @msa2profile_hhblit_arr);
			push (@atchley, @{ $factors{$aa} });
			push (@prob, @{ $probs[$ii] });
		}

		# Adds debugging line. Do not activate for actual trial
		#push (@printarray, "Line $jj: $pssmlines[$jj]\n");

		# Adds target related features
		push (@printarray, @targetfeat) unless ($test);

		# Adds boundary and gap features
		push (@printarray, @boundfeat) if ($bound);

		# Adds prob related features
		push (@printarray, @prob) if (@probfiles);

		# Adds seq related features
		push (@printarray, @seqfeat) if ($seq);

		# Adds PSSM related features
    #print "Length of pssm: ",@pssmfeat."\n";
		push (@printarray, @pssmfeat) if ($pssm);
   
		push (@printarray, @pssmfeat) if ($pssm_uniref50);
   
		# Adds atchley features
		push (@printarray, @atchley) if ($atch);
   
		# Adds HMM Emission related features
    #print "Length of pssm: ",@pssmfeat."\n";
		push (@printarray, @hmmEmfeat) if ($hmmEm);
   
		# Adds HMM Transition related features
    #print "Length of pssm: ",@pssmfeat."\n";
		push (@printarray, @hmmTrfeat) if ($hmmTr);
   
		# Adds sspro msa profile related features
    #print "Length of pssm: ",@pssmfeat."\n";
		push (@printarray, @msa2profile_ssprofeat) if ($aln_sspro);
   
		# Adds sspro msa profile related features
    #print "Length of pssm: ",@pssmfeat."\n";
		push (@printarray, @msa2profile_sspro_uniref50feat) if ($aln_sspro_unif50);
   
		# Adds sspro msa profile related features
    #print "Length of pssm: ",@pssmfeat."\n";
		push (@printarray, @msa2profile_sspro_uniref90feat) if ($aln_sspro_unif90);
   
		# Adds hhblit msa profile related features
    #print "Length of pssm: ",@pssmfeat."\n";
		push (@printarray, @msa2profile_hhblitfeat) if ($aln_hhblits);


		# Concatenates features onto a feature line
		my $line;
		for (@printarray){
			$line .= "$_ ";
		}
		$line .= "\n";
		push (@TFlines, $line);
		SKIP:
	}

	return @TFlines;
}

sub check_probs {
	my ($str, @probs) = @_;

	for (my $ii=0; $ii<@probs; $ii++){
		next if ($probs[$ii] eq "-");
		my @vals = @{ $probs[$ii] };
		foreach (@vals){
			if ($_ > 1 || $_ < 0){
				print "$str: format wrong: @vals\n";
				return;
			}
		}
	}
}

sub print_features {
	my ($file, @features) = @_;

	open (OUT, ">$file") or return "err: print_features: Couldn't open file $file";
	foreach (@features){
		print OUT "$_";
	}
	close OUT or return "err: print_features: Couldn't close file $file";
}

################################################################
# Name:	linearize
# In:	lines : array of pssm file lines
# Out:	pssm, inf : arrays containing relevant PSSM information
################################################################

sub linearize {
	my @lines = @_;
	my @pssm =();
	my @inf = ();

	# Get rid of header
	shift @lines; shift @lines; shift @lines;

	# Loop over all lines until residues are depleted
	while(my $line_t = shift @lines) {
		if($line_t =~ m/^\s*$/) {
			last; 
		}

		my @fields = ();
		my $start = 9;
		my $len = 3;

		# Xs are treated differently because although the first part
		# of the PSSM matrix is all '-1', some additional data can
		# be derived from the second part of the matrix
		if (substr($line_t, 6, 1) eq "X"){
			$start = 71;
			$len = 4;
		}
		for(my $i =0; $i < 20; $i++) {
			my $val = substr($line_t, $start+$i*$len, $len);

			# Values are transformed to the range [0, 1]
			my $linvalue = linscale($val, $start); 
			push(@fields, sprintf("%.1f", $linvalue));
		}
		push(@pssm, [@fields]);

		# Gathers inf data and scales
		my $inf_t = substr($line_t, 151, 5);
		if($inf_t > 6) {
			push(@inf, "1.0000");
		} elsif( $inf_t <= 0) {
			push(@inf, "0.0000");
		} else {
			push(@inf, sprintf("%.4f", ($inf_t/6.0)));
		}
	}
	return(\@pssm, \@inf);
}

################################################################
# Name:	linscale
# In:	value: a number from PSSM
# 	source: number representing whether residue is "X"
# Out:	value transformed to [0,1]
################################################################

sub linscale {
	my ($value, $source) = @_;

	# Depending on what part of the matrix the value came from,
	# the domain of the scaling function must be adjusted
	return $value/100 if ($source == 71);
	return 0 if ($value < -5);
	return 1 if ($value > 5);
	return ($value*.1 + .5);
}

################################################################
# Name:	normalize
# In:	line: line from a prob file 
# Out:	array of nums from file, normalized to sum to 1
################################################################

sub normalize {
	my ($line) = @_;
	my @nums = split(/ /, $line);

	return normalize_array(@nums);
}

################################################################
# Name:	normalize_array
# In:	nums: array of numbers from a prob file 
# Out:	array of nums normalized to sum to 1
################################################################

sub normalize_array {
	my @nums = @_;
	my $sum = 0;
	foreach (@nums){
		$_ = 0 if ($_ =~ /nan/ || $_ < 0);
		$sum += $_;
	}

	my @newnums;
	foreach my $num (@nums){
		$num = $num/$sum if ($sum>1);
		$num = sprintf("%.6f", $num);
		push(@newnums, $num);
	}

	return @newnums;
}

################################################################
# Name:	get_avg
# In:	length: amount of lines expected
#	fileinfo_refsref: refs to arrays of prob file info
# Out:	string averages over each col over each row
################################################################

sub get_avg {
	my ($length, $fileinfo_refsref) = @_;
	my @fileinfo_refs = @{ $fileinfo_refsref };
	my $targsize = get_targsize();
	my @probs;
	my $count = 0;

	# This loop generates dummy prob information in absence of files
	if (@fileinfo_refs == 0){
		for (my $ii=0; $ii<$length; $ii++){
			my @array;
			for (my $jj=0; $jj<$targsize; $jj++){
				push (@array, "");
			}
			push (@probs, \@array);
		}
		return @probs;
	}

	# Averages are calculated regardless of how many files are input
	foreach (@fileinfo_refs){
		my @filelines = @{ $_ };

		for (my $ii=0; $ii<@filelines; $ii++){
			my @nums = @{ $filelines[$ii] };

			if (defined($probs[$ii])){ # Adds and averages 
				my @oldavg = @{ $probs[$ii] };
				my @newavg;
				for (my $jj=0; $jj<@nums; $jj++){
					$newavg[$jj] = ($count*$oldavg[$jj]+$nums[$jj])/($count+1);
				}
				$probs[$ii] = \@newavg;
			}
			else { # Initializes
				$probs[$ii] = \@nums;
			}
		}
		$count++;
	}

	# Concatenates probs into a string for each residue
	for (my $ii=0; $ii<@probs; $ii++){
		my @nums = @{ $probs[$ii] };

		my @newline;
		foreach (@nums){
			push (@newline, sprintf("%.6f", $_));
		}
		$probs[$ii] = \@newline;
	}

	return @probs;
}

################################################################
# Name: train_DN
# In:	dir : model directory
#	training_file : training file name (in dir)
#	arch : designated architecture
# Out:	Nothing or error. DN is trained if successful.
#
# Subs: ../feat-2-npy.py : python script converts text to numpy
#	../train-DN_arch.py : python script that trains DN
################################################################

sub train_DN {
	my ($dir, $training_file, $arch) = @_;

	my $name = substr($training_file, 0, -7);
	my @paths = split (/\//, $name);

	my $infile = $training_file;
	my $npyfile = "$name.npy";
	my $modelfile = $dir . "$paths[$#paths].model.dat";

	return("err: train_DN: Dir $dir not found") unless (-d $dir);
	return("err: train_DN: TF $training_file not found") unless (-f $training_file);
	
	# Training File converted to numpy first, then used as input for DN training
	`python2 ./feat-2-npy.py $infile $npyfile 1>&2` unless (-f $npyfile);
	`python2 ./trainDN.py $npyfile $modelfile $arch 1>&2`;

	#`rm $training_file`;
	#`rm $npyfile`;
}

################################################################
# Name: test_DN
# In:	model file : model file
#	feat dir : location of testfeature files
#	out dir : dir of desired prob files
#	list : array list of proteins to test
#	targ : number corresponding to the amount of targets
# Out:	Nothing or error. Files of probs are made if successful
# Subs:	../feat-2-npy.py : python script converts text to numpy
#	../test-DN.py : python script that makes prediction
################################################################

sub test_DN {
	my ($model, $featdir, $outdir, $list_ref, $targ) = @_;
	my @errlines;

	my @list = @{ $list_ref };

	$outdir .= "/" unless ($outdir =~ /\/$/);
	$featdir .= "/" unless ($featdir =~ /\/$/);
	my $npydir = $outdir . "npy/";
	`mkdir $npydir` unless (-d $npydir);

	my $count = 0;
	my $percent = 0;
	my $length = @list;

	# Makes prediction for each prot in input list
	foreach my $prot (@list){
		if ($count/$length*100 > $percent && $percent < 100){
			#print STDERR "$percent % complete.\n";
			$percent += 50;
		}
		$count++;

		my $featfile = $featdir . "$prot.feat";
		my $npy = $npydir . "$prot.npy";
		my $predfile = $outdir . "$prot.pred";
		my $probfile = $outdir . "$prot.prob";

		# Skips if feat file doesn't exist, or if prob file already does
		unless (-f $featfile){
			push (@errlines, "File $featfile not found");
			next;
		}
		next if (-f $probfile);

		# Feat files converted to numpy first, then used with DN
		`python2 ./feat-2-npy.py $featfile $npy`;
		my $output = `python2 ./testDN.py $npy $model $predfile $probfile $targ`;
	}
	#print STDERR "100% complete.\n";

	`rm -r $npydir`;

	return @errlines;
}

################################################################
# Name: predict_SS
# In:	model file : model file
#	feat file  : feature file
#	out dir	   : directory of output
# Out:	Name of the files made or error
# Subs:	./feat-2-npy.py : python script converts text to numpy
#	./test-DN.py    : python script that makes prediction
################################################################

sub predict_SS {
	my ($model, $featfile, $outdir, $tag) = @_;
	my @paths = split('/', $featfile);
	my @info = split(/\./, $paths[$#paths]);
	my $prot = join(/./, @info[$#info-1]);
	$prot .= $tag;

	my $npydir = $outdir . "npy/";
	`mkdir $outdir` unless (-d $outdir);
	`mkdir $npydir` unless (-d $npydir);

	my $npy = $npydir . "$prot.npy";
	my $probfile = $outdir . "$prot.prob";
	my $predfile = $outdir . "$prot.pred";

	unless (-f $featfile){
		return ("err: predict_SS: File $featfile not found");
	}
	
	`python2 ./feat-2-npy.py $featfile $npy`;
	`python2 ./test-DN.py $npy $model $predfile $probfile 3`;

	`rm -r $npydir`;
	return ($probfile, $predfile);
}

sub predict_SS2 {
	my ($model, $weight, $featfile, $outdir, $tag) = @_;
	my @paths = split('/', $featfile);
	my @info = split(/\./, $paths[$#paths]);
	my $prot = join(/./, @info[$#info-1]);
	$prot .= $tag;

	#my $npydir = $outdir . "npy/";
	`mkdir $outdir` unless (-d $outdir);
	#`mkdir $npydir` unless (-d $npydir);

	my $probfile = $outdir . "$prot.prob";
	my $predfile = $outdir . "$prot.pred";

	unless (-f $featfile){
		return ("err: predict_SS: File $featfile not found");
	} 
	`python2 $GLOBAL_PATH/scripts/test_dnss_single.py $model $weight $featfile $predfile $probfile &> $outdir/$tag.log`;
	return ($probfile, $predfile);
}

################################################################
# Name: all_dnss_files
# In:	predir_ref : array containing one or multiple prediction dirs
#	pssmdir : pssm dir
#	outdir : outdir path
#	list_ref : array of prots to make into files
#	tag : info to be printed in header of file
# Out:	None. Writes dnss files at designated outfile
################################################################

sub all_dnss_files {
	my ($predir_ref, $pssmdir, $outdir, $list_ref, $tag) = @_;
	my @predirs = @{ $predir_ref };
	my @list = @{ $list_ref };
	$pssmdir .= "/" unless $pssmdir =~ m/\/$/;
	$outdir .= "/" unless $outdir =~ m/\/$/;
	my $ssadir = "$GLOBAL_PATH/DNSS_dataset/ssa/";

	return "err: all_dnss_files: $pssmdir not found" unless (-d $pssmdir);
	return "err: all_dnss_files: $outdir not found" unless (-d $outdir);

	my $count = 0;
	my $percent = 0;
	foreach my $prot (@list){
		if ($count/@list*100 > $percent && $percent < 100){
			#print "$percent % complete.\n";
			$percent += 50;
		}
		$count++;

		my $pssmfile = $pssmdir . "$prot.pssm";
		my $header = "$prot : $tag";
		my $outfile = $outdir . "$prot.dnss";
		my $ssafile = $ssadir . "$prot.ssa";
		next if (-f $outfile);

		# Multiple pred files can be used to generate avg prediction
		my @predfiles;
		foreach my $dir (@predirs){
			$dir .= "/" unless $dir =~ /\/$/;
			my $tempfile = $dir . "$prot.prob";
			push (@predfiles, $tempfile) if (-f $tempfile);
		}
		if (@predfiles == 0) {
			print "all_dnss_files: no valid predfiles for $prot\n";
			next;
		}

		# Fcn checks that resulting file has one predicted ss per residue
		my $equal = make_dnss_file(\@predfiles, $pssmfile, $outfile, $header, $ssafile);
		unless ($equal){
			print "all_dnss_files: $prot has different size seq and ss!\n";
		}
	}
	#print "100% complete!\n";
}

################################################################
# Name: make_dnss_file
# In:	predfile_ref : array containing one or multiple prediction dirs
#	pssm : pssm file
#	out : outfile name (including path)
#	header : prot info, which will be printed as header of file
# Out:	None. Writes dnss file at designated outfile
################################################################

sub make_dnss_file {
	my ($predfile_ref, $pssmfile, $outfile, $vdnssfile, $header, $ssafile) = @_;
	my @predfiles = @{ $predfile_ref };
	my @seq;

	# PSSM residue sequence is loaded
	open (IN, $pssmfile) or return "err: make_dnss_file: Couldn't open file $pssmfile";
	<IN>; <IN>; <IN>;
	while ((my $line = <IN>) =~ m/^\s*(\d+)\s*(\w)/){
		my ($num, $aa) = $line =~ /^\s*(\d+)\s*(\w)/;
		push (@seq, $aa);
	}
	close IN or return "err: make_dnss_file: Couldn't close file $pssmfile";

	# If multiple prob files are used, their sum is calculated
	# Then this result is used to convert probabilities to structure
	my @probs = sum_probs(@predfiles);
	my @ss = probs_to_ss(@probs);
	my @best = best_probs(@probs);

	if (@ss != @seq){
#		print "\nBefore alignment:\n";
#		print_comparison(\@seq, \@ss);

		my @ssaseq = load_ssa($ssafile); 
		my @lines = (\@seq, \@ssaseq);
		my ($new_pssm_ref, $new_ssa_ref) = align_seqs(@lines);
		@seq = @{ $new_pssm_ref };
		@ssaseq = @{ $new_ssa_ref };

		my @newss;
		for (my $ii=0; $ii<@ssaseq; $ii++){
			if ($ssaseq[$ii] eq "-"){
				push (@newss, "-");
			}
			else {
				push (@newss, shift(@ss));
			}
		}
		@ss = @newss;

#		print "\nAfter alignment:\n";
#		print_comparison(\@ssaseq, \@seq, \@ss);
#		print "Continue? (y/n)\n";
#		my $response = <STDIN>;
#		unless ($response =~ /y/){
#			exit;
#		}
	}

	# Format of DNSS file is as follows:
	# >ProtName-A  	(given by input "header")
	# ACDFGHIKL	(sequence found in PSSM file)
	# CCCEEEHHH	(secondary structure prediction)
	open (OUT, ">$outfile") or return "err: make_dnss_file: Coudln't open file $outfile";
	print OUT "$header\n";
	print OUT @seq;
	print OUT "\n";
	print OUT @ss;
	print OUT "\n";
	close OUT or return "err: make_dnss_file: Coudln't close file $outfile";

	# Format of VDNSS file is as follows:
	# >ProtName-A  		(given by input "header")
	# A	C	.6969	(First col is residue)	
	# C	C	.8797	(Second col is ss pred)
	# D	H	.9900	(Third col is confidence)
	# F	E	.7989	(File is Tab separated)

	open (OUT, ">$vdnssfile") or return "err: make_dnss_file: Couldn't open file $vdnssfile";
	print OUT "$header\n";
	for (my $ii=0; $ii<@seq; $ii++){
		print OUT "$seq[$ii]\t$ss[$ii]\t$best[$ii]\n";
	}
	close OUT or return "err: make_dnss_file: Coudln't close file $vdnssfile";

	return 1;
}
#############################
#############################
##    SCORE_CALCULATIONS   ##
#############################
#############################

################################################################
# Name: score_probs
# In:	predir_ref : array ref containing one or multiple dirs of predictions
#		note that dirs need contain predictions of the same proteins (redundant)
#	list_ref : array ref of list of proteins to be scored
#	tag : the desired name of the output score file
# Out:	results: array containing (AvgQ3, AvgSov);
################################################################

sub score_probs {
	my ($predir_ref, $list_ref, $reduced, $tag) = @_;
	
	my $ssadir = "$GLOBAL_PATH/DNSS_dataset/ssa/";
	my $pssmdir = "$GLOBAL_PATH/DNSS_dataset/pssm/";
	my $dnssdir = "$GLOBAL_PATH/tmp/";

	# Temporary dnss files are made with predictions
	print "Making temp dnss files...\n";
	all_dnss_files($predir_ref, $pssmdir, $dnssdir, $list_ref);

	# No scorefile is made if no tag is specified
	my $scorefile = "";
	if ($tag){ $scorefile = "$GLOBAL_PATH/DNSS_dataset/ss_pred/scores/$tag.score"; }

	# Predictions are scored
	print "Scoring...\n";
	my (@results) = all_score ($ssadir, $dnssdir, $scorefile, $list_ref);

	# Temp files are removed
	`rm $dnssdir/*`;
	return @results;
}

################################################################
# Name: score_dnss
# In:	dnssdir : dir containing dnss files
#	list_ref : array ref of list of proteins to be scored
#	tag : the desired name of the output score file
# Out:	results: array containing (AvgQ3, AvgSov);
################################################################

sub score_dnss {
	my ($dnssdir, $list_ref, $tag) = @_;
	my $ssadir = "$GLOBAL_PATH/DNSS_dataset/ssa/";

	# No scorefile is made if no tag is specified
	my $scorefile = "";
	if ($tag){ $scorefile = "$dnssdir/$tag.score"; }
	#print "Saving to $scorefile\n";
	# Predictions are scored
	my @results = all_score ($ssadir, $dnssdir, $scorefile, $list_ref);
	return @results;	
}

################################################################
# Name: all_score
# In:	refdir : directory of ssa files
#	evadir : directory of eva files
#	sovfile : file to store scoring results
#	list : list of proteins
# Out:	Q3score, Sovscore
#	prints more information to sovfile
################################################################

sub all_score {
	my ($refdir, $evadir, $sovfile, $list_ref, $eva_ext) = @_;
	my @list = @{ $list_ref };
	$eva_ext = "dnss" unless ($eva_ext);

	my %Scorelines; 	# Holds ProtID => (scores)
	my %Qsegments;		# Holds percentile => # of occurrences
	my %Ssegments;		# Holds ^ for Sov scores instead of Q
	my $Qmatch = 0;		# Holds correctly scored residues (ignoring chains)
	my $Qtot = 0;		# Holds total residues compared (ignoring chains)
	my $Qsum = 0;		# Holds sum of Q scores by chain
	my $Sovsum = 0;		# Holds sum of Sov scores by chain
	my $num = 0;		# Holds number of chains scored

	my $count = 0;
	my $percent = 0;
	foreach my $prot (@list) {
		if ($count/@list*100 > $percent && $percent < 100){
			#print STDERR "$percent % completed.\n";
			$percent += 50;
		}
		$count ++;

		# Reference sequence and ss are loaded from SSA file
		my $reffile = $refdir . "$prot.ssa";
		unless (-f $reffile) {
			print "$prot: No SSA found. Skipping\n";
			next;
		}
		my @ref_seq = load_ssa($reffile);
		my @ref_ss = import_ssa($reffile);
		next if (check_err($ref_ss[0]));

		# Evaluee sequence and ss are loaded from eva file
		my $evafile = $evadir . "$prot.$eva_ext";
		unless (-f $evafile) {
			print "$prot: No eva file ($evafile) found. Skipping\n";
			next;
		}
		my @eva_seq = load_seq($evafile);
		my @eva_ss = import_struct($evafile);
		next if (check_err($eva_ss[0]));

		if (@eva_ss == 0){
			print "all_score: $prot eva file ss missing\n";
			next;
		}
		if (@ref_ss == 0){
			print "all_score: $prot ref file ss missing\n";
			next;
		}

		# Seqs are aligned to ensure that ss comparisons are relevant
		my (@align_arrays) = get_align_arrays(\@ref_seq, \@eva_seq);

		# Alignment is applied to seqs (for diagnostic purposes only)
		my @seqs = (\@ref_seq, \@eva_seq);
		my ($ref_ref, $eva_ref) = apply_alignment(\@seqs, \@align_arrays);
		@ref_seq = @{ $ref_ref };
		@eva_seq = @{ $eva_ref };

		# Alignment is applied to SSs
		my @seqs = (\@ref_ss, \@eva_ss);
		my ($ref_ref, $eva_ref) = apply_alignment(\@seqs, \@align_arrays);
		@ref_ss = @{ $ref_ref };
		@eva_ss = @{ $eva_ref };

		# Optional diagnostic print statements
#		print "\nComparing the following sequences/structures\n";
#		my @names = ("Ref Seq", "Ref SS", "Eva Seq", "Eva SS");
#		my @seqs = (\@ref_seq, \@ref_ss, \@eva_seq, \@eva_ss);
#		print_comparison_named(\@names, \@seqs);

		# Scores are calculated for this prot
		my ($thismat, $thistot, $thisQ, $thisSov) = calc_sov(\@ref_ss, \@eva_ss);
		next if (check_err($thismat));
		
		# Sums are updated
		$Qmatch += $thismat;
		$Qtot += $thistot;
		$Qsum += $thisQ;
		$Sovsum += $thisSov;
		$num++;

		# Score is placed into a percentile group for Q and Sov
		my $Qseg = int($thisQ/10);
		$Qsegments{$Qseg}++;

		my $Sseg = int($thisSov/10);
		$Ssegments{$Sseg}++;

		# Print file information is stored
		my $scoreout = "\n>$prot\n";
		$scoreout .= "!Q3  = $thisQ\n";
		$scoreout .= "!Sov = $thisSov\n";
		$Scorelines{$prot} = $scoreout;
	}
	#print "100% complete!\n";

	# Averages are calculated
	my $Qtotscore = $Qmatch/$Qtot*100;
	my $Qscore = $Qsum/$num;
	my $Sovscore = $Sovsum/$num;

	# Score file is written
	my @comments;
	push (@comments, "#Total Matches: $Qmatch\tTotal residues: $Qtot\n");
	push (@comments, "#Total Q3 score    : $Qtotscore\n");
	push (@comments, "#Average Q3 score  : $Qscore\n");
	push (@comments, "#Average Sov score : $Sovscore\n");
	push (@comments, "#\n#Number of proteins in each score range:\n");
	for (my $ii=0; $ii<10; $ii++){
		my $low = $ii. "0";
		my $high = $ii+1 . "0";
		my $range = "\t$low-$high |";
		push (@comments, $range);
	}
	push (@comments, "\n#Q3:");
	for (my $ii=0; $ii<10; $ii++){
		push (@comments, "\t$Qsegments{$ii}");
	}
	push (@comments, "\n#Sov:");
	for (my $ii=0; $ii<10; $ii++){
		push (@comments, "\t$Ssegments{$ii}");
	}
	push (@comments, "\n");

	my @keys = sort (keys %Scorelines);

	open (OUT, ">$sovfile") or return ($Qscore, $Sovscore);
	print OUT @comments;
	foreach my $prot (@keys) {
		print OUT $Scorelines{$prot};
	}
	close OUT or return ("err: all_score: Couldn't close file $sovfile");
	
	return ($Qscore, $Sovscore);
}

################################################################
# Name: sum_probs
# In:	probfiles : any number, but usually 3, prob files
# Out:	sumprobs : array containing the sum of probabilies 
# 
# Subs:	probs_to_ss : Converts probs to an ss char (C, E, H)
#	add_probs : adds a line of probs to an existing line
################################################################

sub sum_probs {
	my (@probfiles) = @_;
	my @sumprobs;
	my $targs = get_targsize();
	my $linecount = 0;

	# Sets template for initializing probs
	my @zeros;
	for (0..$targs-1){
		push (@zeros, 0);
	}

	# Loops over any number of prob files. Commonly 1 or 3 files are used
	foreach my $file (@probfiles) {
		open (IN, $file) or return ("err: sum_probs: Unable to open file $file");
		
		my @fileprobs;
		my $count = 0;
		while (my $line = <IN>){
			my @lineprobs = split(/ /, $line);
			my $boost = @lineprobs/$targs;
			my $half = ($boost-1)/2;

			# Adds probs based on how much boosting was done
			for (my $jj=0; $jj<$boost; $jj++){
				next if ($count+$jj-$half < 0);

				my $curr = \@zeros;
				if ($fileprobs[$jj+$count-$half]){
					$curr = $fileprobs[$jj+$count-$half];
				}

				# Separates one instance of a prediction and adds
				my $beg = $jj*$targs;
				my $end = ($jj+1)*$targs-1;
				my @add = @lineprobs[$beg..$end];
				my @probs = add_probs($curr, \@add);
				$fileprobs[$jj+$count-$half] = \@probs;
			}
			$count++;
		}

		close IN or return ("err: sum_probs: Unable to close file $file");

		# Normalizes probs to sum to 1 and adds this file's probs to totals
		for (my $ii=0; $ii<@fileprobs; $ii++){
			my @nums = @{ $fileprobs[$ii] };
			@nums = normalize_array(@nums);
			$fileprobs[$ii] = \@nums;

			my $curr = \@zeros;
			$curr = $sumprobs[$ii] if ($sumprobs[$ii]);
			my @sum = add_probs($curr, $fileprobs[$ii]);
			$sumprobs[$ii] = \@sum;
		}

		$linecount = $count-1;
	}

	# Normalizes probs to sum to 1 since they are sums
	for (my $ii=0; $ii<@sumprobs; $ii++){
		my @nums = @{ $sumprobs[$ii] };
		@nums = normalize_array(@nums);
		$sumprobs[$ii] = \@nums;
	}

	return @sumprobs[0..$linecount];
}

################################################################
# Name:	add_probs
# In:	2 array refs of probs that are to be added
# Out:	array containing added probs
################################################################

sub add_probs {
	my ($current_ref, $add_ref) = @_;
	my @current = @{ $current_ref };
	my @addition = @{ $add_ref };
	my @sum;

	for (my $ii=0; $ii<@addition; $ii++){
		$sum[$ii] = $current[$ii] + $addition[$ii];
	}
	return @sum;
}

################################################################
# Name:	probs_to_ss
# In:	array of probs
# Out:	array of predicted SS chars
################################################################

sub probs_to_ss {
	my (@probs) = @_;
	my @ss;
	my @options = ("C", "E", "H");

	foreach (@probs){
		my @nums = @{$_};
		my @sorted = sort {$b <=> $a} @nums;

		# Note that here the order of comparison matters in a rare
		# amount of cases. Sometimes probs are equal (e.g. "1 1 1")
		# in which case the current implementation will assign ties to
		# SS in this priority order: C, E, H
		for (my $ii=0; $ii<@nums; $ii++){
			if ($nums[$ii] == $sorted[0]){
				push (@ss, $options[$ii]);
				goto NEXT;
			}
		}
		NEXT:
	}
	return @ss;
}

################################################################
# Name:	best_probs
# In:	array of probs
# Out:	value of most likely prob
################################################################

sub best_probs {
	my (@probs) = @_;
	my @best;

	foreach (@probs){
		my @nums = @{$_};
		my @sorted = sort {$b <=> $a} @nums;
		push (@best, $sorted[0]);
	}
	return @best;
}

#########################################################
# Name: calc_sov
# In:	ref_ref : array ref containing reference ss sequence
#	eva_ref : array ref containing evaluee ss sequence
#		note: these refs should be fully aligned,
#			though need not be gapless
# Out:	array containing (Qmatches, Qtotal, Q3score, Sovscore) 
#
# subs: len, seg_len, max, min, separate, get_first, get_last,
#	import_ssa, import_ss_sa, import_pred, calculate_RQ3, 
#	calculate_C, calculate_Q3, calculate_sov
#########################################################

sub calc_sov {
	my ($ref_ref, $eva_ref) = @_;
	my @raw_ref_ss = @{ $ref_ref };
	my @raw_eva_ss = @{ $eva_ref };
	my ($Qmatch, $Qtotal, $Qscore, $Sovscore);

	if (@raw_ref_ss != @raw_eva_ss){
		return "err: calc_sov: input arrays not of equal size\n";
	}

	#Removes non-existent structure from comparison.
	#Thus, SS prediction of AAs that are removed in
	#known structure are not evaluated.
	my (@ref_ss, @eva_ss);
	for (my $ii=0; $ii<@raw_ref_ss; $ii++){
		if ($raw_ref_ss[$ii] ne "-" && $raw_eva_ss[$ii] ne "-"){
			push (@ref_ss, $raw_ref_ss[$ii]);
			push (@eva_ss, $raw_eva_ss[$ii]);
		}
	}

	# Optional diagnostic print statements
#	print "Reduced structures to be scored:\n";
#	my @names = ("Ref", "Eva");
#	my @seqs = (\@ref_ss, \@eva_ss);
#	print_comparison_named(\@names, \@seqs);

	# Q scores are calculated
	my @Qscores;
	($Qscore, $Qmatch, $Qtotal, @Qscores) = calculate_Q3(\@ref_ss, \@eva_ss);

	# Sov score is calculated
	$Sovscore = calculate_sov(\@ref_ss, \@eva_ss);

	return ($Qmatch, $Qtotal, $Qscore, $Sovscore);
}

sub calculate_sov {
	my ($ref_ref, $eva_ref) = @_;
	my @ref_ss = @{ $ref_ref };
	my @eva_ss = @{ $eva_ref };
	
	my @SSs = ("C", "E", "H");

	my $Ntotal = 0;
	my $Sumtotal = 0;
	my @Sovscores;

	# A score is calculated for the correctness of each kind of secondary structure
	foreach my $ss (@SSs){
		my $ref_segstart = -1;
		my $ref_segend = -1;
		my $eva_segstart = -1;
		my $eva_segend = -1;

		my @ref_segs;
		my @eva_segs;

		# Arrays are made contianing the segments of current secondary structure
		for (my $ii=0; $ii<@ref_ss; $ii++){
			if ($ref_ss[$ii] eq $ss){
				if ($ref_segstart == -1){
					$ref_segstart = $ii;
				}
				$ref_segend = $ii;
			}
			else {
				if ($ref_segstart != -1){
					my $segstring = "$ref_segstart,$ref_segend";
					push (@ref_segs, $segstring); 
					$ref_segstart = -1;
				}
			}

			if ($eva_ss[$ii] eq $ss){
				if ($eva_segstart == -1){
					$eva_segstart = $ii;
				}
				$eva_segend = $ii;
			}
			else {
				if ($eva_segstart != -1){
					my $segstring = "$eva_segstart,$eva_segend";
					push (@eva_segs, $segstring); 
					$eva_segstart = -1;
				}
			}
		}
		if ($ref_segstart != -1){
			my $segstring = "$ref_segstart,$ref_segend";
			push (@ref_segs, $segstring);
		}
		if ($eva_segstart != -1){
			my $segstring = "$eva_segstart,$eva_segend";
			push (@eva_segs, $segstring);
		}

		# Arrays are made containing the overlapping and non-overlapping regions
		my @overlaps;
		my @alones;
		my $rmark = 0;

		foreach my $evaseg (@eva_segs){
			my $done = 0;
			my $added = 0;
			while (!$done){
				if ($rmark == @ref_segs){
					push (@alones, $evaseg);
					$done = 1;
					next;
				}
				my $refseg = $ref_segs[$rmark];
				if (get_last($evaseg) < get_first($refseg)){
					if (!$added){
						push (@alones, $evaseg);
					}
					$done = 1;
				}
				elsif (get_first($evaseg) > get_last($refseg)){
					$rmark++;
					next;
				}
				else {
					my $segstring = "$evaseg;$refseg";
					push (@overlaps, $segstring);
					$added++;
					if (get_last($evaseg) > get_last($refseg)){
						$rmark++;
					}
					elsif (get_last($evaseg) < get_last($refseg)){
						$done = 1;
					}
					else {
						$rmark++;
						$done = 1;
					}
				}
			}
		}

		# Here the individual normalization is calculated:
		my $Ncount = 0;
		foreach my $segpair (@overlaps){
			my ($fl, $fr, $sl, $sr) = separate($segpair);

			$Ncount += len($fr, $fl);
		}
		foreach my $seg (@alones){
			my ($fl, $fr) = split(/,/, $seg);

			$Ncount += len($fr, $fl);
		}
		$Ntotal += $Ncount;

		# Here the individual summations are calculated.
		my $thisSum = 0;
		for (@overlaps){
			my ($fl, $fr, $sl, $sr) = separate($_);
			my $s1 = "$fl,$fr";
			my $s2 = "$sl,$sr";

			my $minov = len(max($fl, $sl), min($fr, $sr));
			my $maxov = len(min($fl, $sl), max($fr, $sr));

			my $diff = $maxov - $minov;
			my $int1 = int(seg_len($s1)/2);
			my $int2 = int(seg_len($s2)/2);
			my $delta = min($diff, $minov, $int1, $int2);

			my $numer = $minov + $delta;
			my $overlapNum = ($numer / $maxov) * seg_len($s1);
			$thisSum += $overlapNum;
		}
		$Sumtotal += $thisSum;
		if ($Ncount == 0){
			push (@Sovscores, 100);
		}
		else {
			push (@Sovscores, $thisSum/$Ncount*100);
		}
	}

	if ($Ntotal == 0){
		return ("err: calculate_sov: File insufficient");
	}
	my $Sovscore= 100*$Sumtotal/$Ntotal;
	return $Sovscore;

}

# Calculates the length between two indices
sub len {
	my ($first, $last) = @_;
	if ($first > $last){
		my $temp = $first;
		$first = $last;
		$last = $temp;
	}
	my $len = $last - $first + 1;
	return $len;
}

# Calculates the length of a segment
sub seg_len {
	my ($seg) = @_;
	my ($first, $last) = split (/,/, $seg);

	my $len = $last - $first + 1;
	return $len;
}

# Calculates the max of two variables
sub max {
	my ($x, $y) = @_;
	my $max = ($x, $y)[$x < $y];
	return $max;
}

# Calculates the min of any number of variables
sub min {
	my (@nums) = @_;
	my $min = $nums[0];
	for (@nums){
		$min = $_ if ($_ < $min);
	}
	return $min;
}

# Gives the four indices stored in a segment pair
sub separate {
	my ($seg_pair) = @_;

	my ($seg1, $seg2) = split (/;/, $seg_pair);
	my ($first_1, $last_1) = split (/,/, $seg1);
	my ($first_2, $last_2) = split (/,/, $seg2);

	my @indices = ($first_1, $last_1, $first_2, $last_2);
	return @indices;
}

# Returns the first index in a segment
sub get_first {
	my ($seg) = @_;

	my ($first, $last) = split (/,/, $seg);
	return $first;
}

# Returns the last index in a segment
sub get_last {
	my ($seg) = @_;

	my ($first, $last) = split (/,/, $seg);
	return $last;

}

# Determines file type and imports the SS sequence
sub import_struct {
	my ($file) = @_;
	
	my @info = split(/\./, $file);
	my $ext = $info[$#info];

	my @compatible = qw/ssa ss_sa dnss pred horiz/;
	unless (grep($ext, @compatible)) {
		return "err: import_struct: incompatible file type";
	}
	return import_ssa($file) if ($ext eq "ssa");
	return import_ss_sa($file) if ($ext eq "ss_sa");
	return import_dnss($file) if ($ext eq "dnss");
	return import_pred($file) if ($ext eq "pred");
	return import_psi($file) if ($ext eq "horiz");
}

# Imports the SS sequence from an SSA file
sub import_ssa {
	my ($ssafile) = @_;
	my @ssa_ss;

	open (IN, $ssafile) or return("err: import_ssa: Couldn't open file $ssafile");
	while (my $line = <IN>){
		next if ($line =~ m/#/);

		my @cols = split(/\t/, $line);
		push (@ssa_ss, $cols[2]);
	}
	close IN or return ("err: import_ssa: Couldn't close file $ssafile");

	return @ssa_ss;
}

# Imports the SS sequence from an SSPro file
sub import_ss_sa {
	my ($ss_safile) = @_;
	my @ss;

	open (IN, $ss_safile) or return "err: import_ss_sa: Couldn't open file $ss_safile";
	while (my $line = <IN>){
		if ($line =~ m/^[HEC]+$/){
			chomp($line);
			@ss = split(//, $line);
		}
	}
	close IN or return "err: import_ss_sa: Couldn't close file $ss_safile";
	return @ss;
}

# Imports the SS sequence from a DNSS file
sub import_dnss {
	my ($dnssfile) = @_;
	my @ss;

	open (IN, $dnssfile) or return "err: import_dnss: Couldn't open file $dnssfile";
	while (my $line = <IN>){
		if ($line =~ m/^[HEC-]+$/){ #### Note: Doesn't work with other chars
			chomp($line);
			@ss = split(//, $line);
		}
	}
	if (@ss == 0){
		print "WARNING: import_dnss did not find valid sec struct line!\n";
	}
	close IN or return "err: import_dnss: Couldn't close file $dnssfile";
	return @ss;
}

# Imports the SS sequence from a prob file
sub import_pred {
	my ($predfile) = @_;
	my @ss;
	my $amb;
	
	open (IN, $predfile) or return "err: import_pred: Couldn't open file $predfile";
	while (my $line = <IN>){
		next unless ($line =~ m/\d\s\d\s\d/);
		my @cols = split(' ', $line);

		# Note that predictions are made in the following priority,
		# which rarely matters, but sometimes multiple SSs are scored
		# with a prediction of 1
		if ($cols[0] == 1){
			push (@ss, "C");
		}
		elsif ($cols[1] == 1){
			push (@ss, "E");
		}
		elsif ($cols[2] == 1) {
			push (@ss, "H");
		}
		else {
			# This is here to look for predictions that were made in error
			$amb++;
			push (@ss, "!");
		}
	}
	if ($amb) {
#		print STDERR "$predfile contains $amb ambiguous predictions.\n";
	}

	close IN or return "err: import_pred: Couldn't close file $predfile";
	return @ss;
}

# Imports the SS sequence from a psi (.horiz) file
sub import_psi {
	my ($psifile) = @_;
	my @ss;

	open (IN, $psifile) or die "Couldn't open file $psifile\n";
	while (my $line = <IN>){
		next unless ($line =~ m/^Pred:/);
		my ($struct) = $line =~ /([CEH]+)/;
		my @moress = split(//, $struct);
		push (@ss, @moress);
	}
	return @ss;
}

# Calculates the Q3 of reduced sequence. This is obsolete since sequences are
# reduced before calculating the scores at all
sub calculate_RQ3 {
	my ($ref_ref, $eva_ref) = @_;
	my @ref_ss = @{ $ref_ref };
	my @eva_ss = @{ $eva_ref };

	my ($rfirst, $rlast, $efirst, $elast, $refcount, $evacount);
	$rfirst = 0;
	$efirst = 0;
	$refcount = 0;
	$evacount = 0;
	for (my $ii=0; $ii< @ref_ss; $ii++){
		if ($ref_ss[$ii] ne "-") {
			$rfirst = $ii unless ($rfirst);
			$rlast = $ii;
			$refcount ++;
		}
		if ($eva_ss[$ii] ne "-") {
			$efirst = $ii unless ($efirst);
			$elast = $ii;
			$evacount ++;
		}
	}

	my $first = ($rfirst, $efirst)[$rfirst < $efirst];
	my $last = ($rlast, $elast)[$rlast > $elast];
	my @new_ref = @ref_ss[$first..$last];
	my @new_eva = @eva_ss[$first..$last];

	my $num = ($refcount, $evacount)[$refcount > $evacount];
	my $den = ($refcount, $evacount)[$refcount < $evacount];
	my $size = $num/$den*100;
	my ($score, @others) = calculate_Q3(\@new_ref, \@new_eva);

	return ($size, $score);
}

# Calculates the C scores for each secondary structure
sub calculate_C {
	my ($ref_ref, $eva_ref) = @_;
	my @ref_ss = @{ $ref_ref };
	my @eva_ss = @{ $eva_ref };
	my @Cscores;
	my @SSs = qw/C E H/;

	foreach my $ss (@SSs){
		my $pp = 0;
		my $rr = 0;
		my $uu = 0;
		my $oo = 0;
		
		for (my $ii=0; $ii<@ref_ss; $ii++){
			if ($ref_ss[$ii] eq $ss){
				$pp++ if $eva_ss[$ii] eq $ss;
				$oo++ if $eva_ss[$ii] ne $ss;
			} else {
				$rr++ if $eva_ss[$ii] ne $ss;
				$uu++ if $eva_ss[$ii] eq $ss;
			}
		}
		my $score = ($pp*$rr - $uu*$oo)/(sqrt(($pp+$uu)*($pp+$oo)*($rr+$uu)*($rr+$oo)))*100;
		push (@Cscores, $score);
		print "!C$ss   = $score\n";
	}


	return @Cscores;
}

# Calculates the Q3 score of a pair of SS sequences
sub calculate_Q3 {
	my ($ref_ref, $eva_ref) = @_;
	my @ref_ss = @{ $ref_ref };
	my @eva_ss = @{ $eva_ref };

	my $Q3total = 0;
	my $Q3match = 0;

	my @SSs = ("C", "E", "H");
	my @Qscores;
	for (my $jj=0; $jj<@SSs; $jj++){
		my $thistotal = 0;
		my $thismatch = 0;

		for (my $ii=0; $ii<@ref_ss; $ii++){
			next if ($ref_ss[$ii] ne $SSs[$jj]);
			$thistotal++;
			if ($ref_ss[$ii] eq $eva_ss[$ii]){
				$thismatch++;
			}
		}

		if ($thistotal == 0) {
			push (@Qscores, 100);
		}
		else{
			push (@Qscores, $thismatch/$thistotal*100);
		}
		$Q3total += $thistotal;
		$Q3match += $thismatch;
	}
	
	my $Q3score = $Q3match/$Q3total*100;
	return ($Q3score, $Q3match, $Q3total, @Qscores);
}

################################################################
# Name: calc_confusion_matrix
# In:	list: array list of proteins
#	ref_dir: directory of reference files
#	eva_dir: directory of evaluee files
# Out:	Confusion matrix over all residues
################################################################

sub calc_confusion_matrix {
	my ($list_ref, $ref_dir, $ref_ext, $eva_dir, $eva_ext) = @_;
	my @list = @{ $list_ref };
	$ref_dir .= "/" unless ($ref_dir =~ /\/$/);
	$eva_dir .= "/" unless ($eva_dir =~ /\/$/);

	my %temp;
	my $matrix = \%temp;

	foreach my $prot (@list){
		# Reference sequence and ss are loaded from SSA file
		my $reffile = $ref_dir . "$prot.$ref_ext";
		unless (-f $reffile) {
			print "$prot: No SSA found. Skipping\n";
			next;
		}
		my @ref_seq = load_seq($reffile);
		my @ref_ss = import_struct($reffile);
		next if (check_err($ref_ss[0]));

		# Evaluee sequence and ss are loaded from eva file
		my $evafile = $eva_dir . "$prot.$eva_ext";
		unless (-f $evafile) {
			print "$prot: No eva file ($evafile) found. Skipping\n";
			next;
		}
		my @eva_seq = load_seq($evafile);
		my @eva_ss = import_struct($evafile);
		next if (check_err($eva_ss[0]));

		if (@eva_ss == 0){
			print "calc_confusin_matrix: $prot eva file ss missing\n";
			next;
		}
		if (@ref_ss == 0){
			print "calc_confusion_matrix: $prot ref file ss missing\n";
			next;
		}

		# Seqs are aligned to ensure that ss comparisons are relevant
		my (@align_arrays) = get_align_arrays(\@ref_seq, \@eva_seq);

		# Alignment is applied to seqs (for diagnostic purposes only)
		my @seqs = (\@ref_seq, \@eva_seq);
		my ($ref_ref, $eva_ref) = apply_alignment(\@seqs, \@align_arrays);
		@ref_seq = @{ $ref_ref };
		@eva_seq = @{ $eva_ref };

		# Alignment is applied to SSs
		my @seqs = (\@ref_ss, \@eva_ss);
		my ($ref_ref, $eva_ref) = apply_alignment(\@seqs, \@align_arrays);
		@ref_ss = @{ $ref_ref };
		@eva_ss = @{ $eva_ref };

		#my ($ref_ref, $eva_ref) = clean_seqs(\@ref_ss, \@eva_ss);
		#my ($ref_seq_ref, $eva_seq_ref) = clean_seqs(\@ref_seq, \@eva_seq);
		#@ref_ss = @{ $ref_ref };
		#@eva_ss = @{ $eva_ref };
		#@ref_seq = @{$ref_seq_ref };
		#@eva_seq = @{$eva_seq_ref };

		# Optional diagnostic print statements
		#print "\nComparing the following sequences/structures\n";
		#my @names = ("Ref Seq", "Ref SS", "Eva Seq", "Eva SS");
		#my @seqs = (\@ref_seq, \@ref_ss, \@eva_seq, \@eva_ss);
		#print_comparison_named(\@names, \@seqs);

		my $new_matrix = confusion_matrix($ref_ref, $eva_ref);
		$matrix = add_matrix($matrix, $new_matrix);
	}
	return $matrix;
}

################################################################
# Name: confusion_matrix
# In:	ref: reference SS sequence
# 	eva: evaluee SS sequence
# Out:	matrix: hash ref containing the confusion matrix.
################################################################

sub confusion_matrix {
	my ($ref_ref, $eva_ref) = @_;
	my @ref = @{ $ref_ref };
	my @eva = @{ $eva_ref };

	if (@ref != @eva){
		print "Ref and Eva not the same length!\n";
		exit;
	}

	my %matrix;

	for (my $ii=0; $ii<@ref; $ii++){
		if ($matrix{$ref[$ii]}{$eva[$ii]}){
			$matrix{$ref[$ii]}{$eva[$ii]}++;
		}
		else {
			$matrix{$ref[$ii]}{$eva[$ii]} = 1;
		}
	}

	return \%matrix;
}

################################################################
# Name: add_matrix
# In:	2 hash refs for matrices of the same shape
# Out:	hash ref contaiing the summed matrix.
################################################################

sub add_matrix {
	my ($mat_ref1, $mat_ref2) = @_;
	my %mat1 = %{ $mat_ref1 };
	my %mat2 = %{ $mat_ref2 };

	my @refs = keys (%mat2);

	foreach my $ref (@refs){
		my @evas = keys (%{ $mat2{$ref} });
		foreach my $eva (@evas){
			if ($mat1{$ref}{$eva}){
				$mat1{$ref}{$eva} += $mat2{$ref}{$eva};
			}
			else {
				$mat1{$ref}{$eva} = $mat2{$ref}{$eva};
			}
		}
	}
	return \%mat1;
}

################################################################
# Name: print_matrix
# In:	hash reference matrix
# Out:	none. Prints the matrix to stdout
################################################################

sub print_matrix {
	my ($ref) = @_;
	my %matrix = %{ $ref };

	my @keys = qw/C E H/;
	
	print "\tC\tE\tH\n";
	foreach my $ref (@keys){
		print "$ref\t";
		foreach my $eva (@keys){
			print "$matrix{$ref}{$eva}\t";
		}
		print "\n";
	}
}

################################################################
# Name: sort_scores
# In:	input: Array of array refs to scores (e.g. Q3 scores array, Sov array)
# Out:	bestindex: The position in the originally array of the best trial
#	data: An array containing (original index, Q3, Sov) sorted with best first
################################################################

sub sort_scores {
	my (@input) = @_;
	my $numscores = @input;

	my @data = ();
	my @temp = @{ $input[0] };
	my $numtrials = @temp;
	for (my $ii=0; $ii<$numtrials; $ii++){
		my @temp = ($ii);
		push (@data, \@temp);
	}

	for (my $ii=0; $ii<@input; $ii++){
		my @scores = @{ $input[$ii] };
		for (my $jj=0; $jj<@scores; $jj++){
			push ($data[$jj], $scores[$jj]);
		}
	}
	for (my $ii=0; $ii<@data; $ii++){
		my @scores = @{$data[$ii]};
		my $sum = 0;

		for (my $jj=1; $jj<@scores; $jj++){
			$sum += $scores[$jj];
		}

		push (@scores, $sum);
		$data[$ii] = \@scores;
	}

	@data = sort { bycol ($b, $a, $numscores+1) } @data;
	my @temp = @{ $data[0] };
	my $bestindex = $temp[0];
	return ($bestindex, @data);
}

################################################################
# Name: sort_scores_by_rank
# In:	Qref : array ref containing multiple Q scores
#	Sovref: array ref containing multiple corresponding Sov scores
# Out:	bestindex: The position in the originally array of the best trial
#	data: An array containing (original index, Q3, Sov) sorted with best first
################################################################

sub sort_scores_by_rank {
	my ($Qref, $Sovref) = @_;
	my @Qs = @{ $Qref };
	my @Sovs = @{ $Sovref };

	my @data;
	for (my $ii=0; $ii<@Qs; $ii++){
		my @temp = ($ii, $Qs[$ii], $Sovs[$ii]);
		push (@data, \@temp);
	}
	
	@data = sort { bycol ($a, $b, 1) } @data; 

	for (my $ii=0; $ii<@data; $ii++){
		my @temp = @{ $data[$ii] };
		push (@temp, $ii);
		$data[$ii] = \@temp;
	}

	@data = sort { bycol ($a, $b, 2) } @data;

	for (my $ii=0; $ii<@data; $ii++){
		my @temp = @{ $data[$ii] };
		$temp[$#temp] += $ii;
		$data[$ii] = \@temp;
	}

	@data = reverse sort { bycol ($a, $b, 3) } @data;

	my @temp = @{ $data[0] };
	my $bestindex = $temp[0];
	return ($bestindex, @data);
}

################################################################
# Name: bycol - sorting function
# In:	2 array refs, plus integer column
# Out:	Sorts according to value in $array[$column]
# Ex:	@sorted = sort { bycol ($a, $b, 3) } @unsorted;
################################################################

sub bycol {
	my ($a, $b, $col) = @_;
	my @ascores = @{ $a };
	my @bscores = @{ $b };

	my $scorea = $ascores[$col];
	my $scoreb = $bscores[$col];

	return ($scorea <=> $scoreb);
}

#############################
#############################
##    OTHER_FUNCTIONS      ##
#############################
#############################

################################################################
# Name:	load_results
# In:	results.txt file
# Out: 	AvgQ3, AvgSov as loaded from file
################################################################

sub load_results {
	my ($file) = @_;

	open (IN, $file) or return ("err: load_results: couldn't open file $file");
	my $line = <IN>;
	my (@scores) = split(/\t/, $line);
	close IN or return ("err: load_results: couldn't close file $file");

	return @scores;
}

################################################################
# Name:	write_results
# In:	file name, array of results
# Out: 	none or err
################################################################

sub write_results {
	my ($file, @results) = @_;

	open (OUT, ">$file") or return ("err, write_results: couldn't open file $file");
	foreach (@results){
		print OUT "$_\t";
	}	
	close OUT or return ("err: write_results: couldn't close file $file");
}

################################################################
# Name: shuffle
# In:	ref of array to be shuffled
# Out:	shuffles array (no reassignment necessary)
# Ex:	shuffle(\@myarray);
################################################################

sub shuffle {
	my $array = shift;
	for (my $ii=@$array-1; $ii>0; $ii--){
		my $jj = int rand($ii+1);
		next if ($ii == $jj);
		@$array[$ii,$jj] = @$array[$jj,$ii];
	}
}

################################################################
# Name: timeprint
# In:	logfile : File name to which messages will be printed
#	message : String message that will be printed
#	gaps : Will insert this number of blank lines before printing message
#	quiet : Boolean: True prevents messages from being printed to STDERR
# Out:	None. Prints message to STDERR (unless quiet) and to logfile
################################################################

sub timeprint {
	my ($logfile, $message, $gaps, $quiet) = @_;

	my $time = formatted_localtime();
	for (my $ii=0; $ii<$gaps; $ii++){
		print STDERR "\n" unless ($quiet);
		`echo "" >> $logfile`;
	}
	print STDERR "$time $message\n" unless ($quiet);
	`echo "$time $message" >> $logfile`;
}

################################################################
# Name:	check_err
# In:	message: possible error message
# 	logfile: relevant logfile (optional) - also -1 indicates quiet
#	print: Prints error message regardless of presence of "err"
# Out:	1 if error occurred, 0 otherwise
# Prnt:	Prints error message to log file and to std out
################################################################

sub check_err{
	my ($message, $logfile, $print) = @_;

	if ($message =~ m/err/){
		print STDERR "$message\n" unless ($logfile == -1);
		`echo "$message\n" >> $logfile` if (-f $logfile);
		return 1;
	}
	if ($print){
		print STDERR "$message\n";
		`echo "$message\n" >> $logfile` if (-f $logfile);
		return 1;
	}
}

################################################################
# Name: get_targsize
# In:	none
# Out:	the size of a single target. This is done so that the difference
#	between SSPRED, SAPRED and TORPRED can be easily converted
################################################################

sub get_targsize {
	return 3;
}

1;
