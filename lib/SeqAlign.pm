package SeqAlign;

# Module: 	SeqAlign.pm
# Author: 	Matt Spencer
# Made:		~3/1/13
# Last Mod:	11/07/13
#
# This module contains functions helpful for loading amino acid sequences
# from one of a variety of files, and is useful for comparing these 
# sequences, as well as providing a simple alignment that assumes that
# the sequences should be highly similar (less than 20% error).
#
# Functions:
#       get_type : indicates if file is compatible
#	load_seq : loads any valid file
#       load_ssa : gets AA sequence in ssa file
#       load_fasta : gets AA sequence in fasta file
#       load_sspro : gets AA sequence in sspro file
#	load_psi : gets AA sequence in psipred (horiz) file
#       load_pdb : gets AA sequence in pdb file
#       load_pssm : gets AA sequence in pssm file
#       load_dssp : gets AA sequence in dssp file
#	load_dnss : gets AA sequence in dnss file
#	clean_seqs : removes all gaps from a pair of sequences
#	align_seqs : makes two sequences align
#	align_and_print : aligns sequences and prints
#	align_and_print_named : aligns sequences and prints with labels
#       get_align_arrays : makes alignment between two sequences
#	remove_gaps : removes all gaps from a sequence
#	MAT_quick_mapping : maps seq1 to seq2 using recursive greediness
#       MAT_make_mapping : maps seq1 to seq2
#	greedy_align : greedy algorithm to find longest identical subseq
#       ALIGN_make_arrays : makes array alignment
#       ALIGN_remove_gaps : removes unnecessary gaps from alignment
#       align_score : gives percent of aligned residues
#	apply_alignment: makes new sequence based on known alignment
#       matched_score : gives score based on shorter segment
#       print_comparison : Prints AA sequences nicely
#       print_comparison_named : Prints AA sequences nicely with labels
#       print_array : Prints alignment arrays
#       print_separator : Prints dashed separator
#       print_double_separator : Prints two dashed separators
#       print_help : Prints help text
#
# Dependencies:
#	./PDBUtils.pm

use strict;
use PDBUtils;
use Exporter;
use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK =  qw(remove_gaps get_type load_seq load_ssa load_fasta load_sspro load_pdb load_pssm load_dssp load_dnss get_align_arrays get_align_arrays_absolute align_seqs align_and_print align_and_print_named apply_alignment align_score matched_score print_comparison print_comparison_named print_array print_separator print_double_separator print_help MAT_align_seqs clean_seqs);

################################################################
# Name : load_seq
# Takes: string: any valid file path
# Returns: array: AA sequence found in the file or 0
################################################################

sub load_seq {
	my ($file) = @_;

	my $type = get_type($file);
	return 0 if ($type eq "0");

	my @seq;

	@seq = load_pdb($file) if ($type eq "pdb");
	@seq = load_ssa($file) if ($type eq "ssa");
	@seq = load_sspro($file) if ($type eq "ss_sa");
	@seq = load_fasta($file) if ($type eq "fasta");
	@seq = load_pssm($file) if ($type eq "pssm");
	@seq = load_dnss($file) if ($type eq "dnss");
	@seq = load_psi($file) if ($type eq "horiz");

	return @seq;
}


################################################################
# Name : get_type
# Takes: Path name of a file
# Returns: int: 0 if there is no load function available
#	   string: suffix of file if it is compatible.
################################################################

sub get_type {
	my ($file) = @_;

	my @pathinfo = split (/\//, $file);
	my @fileinfo = split (/\./, $pathinfo[$#pathinfo]);
	my $suffix = $fileinfo[$#fileinfo];

	return 0 unless ($suffix eq "pdb" || $suffix eq "ssa" || $suffix eq "ss_sa" || $suffix eq "fasta" || $suffix eq "pssm" || $suffix eq "dnss" || $suffix eq "horiz");
	return $suffix;
}


################################################################
# Name : load_ssa
# Takes: string: ssa file path
# Returns: array: AA sequence found in ssa
################################################################

sub load_ssa {
	my ($ssa) = @_;

	open (SSA, $ssa) or die "Coulnd't open file $ssa\n";
	my @ssaseq;
	while (my $line = <SSA>){
		next if ($line =~ m/^#/);
		my @data = split(/\t/, $line);
		push (@ssaseq, $data[1]);
	}
	close SSA or die "Couldn't close file $ssa\n";

	return @ssaseq;
}

################################################################
# Name : load_fasta
# Takes: string: fasta file path
# Returns: array: AA sequence found in fasta
################################################################

sub load_fasta {
	my ($fasta) = @_;

	open (FASTA, $fasta) or die "Couldn't open file $fasta\n";
	my @fastaseq;
	my $line = <FASTA>;
	while ( $line = <FASTA>){
		last if $line =~ /^>/;
		chomp($line);
		push(@fastaseq, split(//, $line));
	}
	close FASTA or die "Couldn't close file $fasta\n";

	return @fastaseq;
}

################################################################
# Name : load_sspro
# Takes: string: sspro file path
# Returns: array: AA sequence found in sspro
################################################################

sub load_sspro {
	my ($ss_sa) = @_;

	open (SS_SA, $ss_sa) or die "Couldn't open file $ss_sa\n";
	my @ss_saseq;
	my $line = <SS_SA>;
	$line = <SS_SA>;
	chomp($line);
	@ss_saseq = split(//, $line);
	close SS_SA or die "Couldn't close file $ss_sa\n";

	return @ss_saseq;
}

################################################################
# Name : load_psi
# Takes: string: psipred file path
# Returns: array: AA sequence found in psipred
################################################################

sub load_psi {
	my ($psi) = @_;

	open (PSI, $psi) or die "Couldn't open file $psi\n";
	my @psiseq;
	while (my $line = <PSI>){
		next unless ($line =~ m/AA:/);
		my ($seq) = $line =~ /AA: (\w+)/;
		push(@psiseq, split(//, $seq));
	}
	close PSI or die "Couldn't close file $psi\n";

	return @psiseq;
}

################################################################
# Name : load_pdb
# Takes: string: pdb file path
# Returns: array: AA sequence found in pdb
################################################################

sub load_pdb {
	my ($pdb) = @_;

	open (PDB, $pdb) or die "Coudln't open file $pdb\n";

	my @pdbseq;
	my $resnum = 0;
	while (my $line = <PDB>){
		unless ($line =~ m/^ATOM/){
			next;
		}
		my $thisres = PDBUtils::get_res_num($line);
		if ($resnum == $thisres){
			next;
		}
		my $res = PDBUtils::get_res($line);
		my $aa = PDBUtils::get_res_code($res);
		push (@pdbseq, $aa);
		$resnum = $thisres;
	}

	close PDB or die "Coudlnt' close file $pdb\n";
	return @pdbseq;
}

################################################################
# Name : load_pssm
# Takes: string: pssm file path
# Returns: array: AA sequence found in pssm
################################################################

sub load_pssm {
	my ($pssm) = @_;
	open (PSSM, $pssm) or die "Couldn't open pssm file $pssm\n";

	my @pssmseq;
	<PSSM>; <PSSM>; <PSSM>;
	while (my $line = <PSSM>){
		last if ($line =~ m/^\s*$/);
		my ($res) = $line =~ /^\s+\d+\s+(\w)/;
		push (@pssmseq, $res);
	}

	close PSSM or die "Couldn't close file $pssm\n";
	return @pssmseq;
}

################################################################
# Name : load_dssp
# Takes: string: dssp file path
# Returns: array: AA sequence found in dssp
################################################################

sub load_dssp {
	my ($dssp) = @_;
	open (DSSP, $dssp) or die "Couldn't open dssp file $dssp\n";
	my @dsspseq;

	my $lastID = 0;
	my $skip=0;
	my $line;
	until ($line =~ m/#/){
		$line = <DSSP>;
	}
	while ($line = <DSSP>){
		my ($id, $aa, $ss, $sa) = $line =~
			/^\s*\d+\s+(\d+)...(.)..(.).{9}\s*\w+\s+\w+\s+(\d+)/;

		unless ($aa){
			$skip=1;
			next;
		}
		elsif ($id != $lastID+1){
			$skip=1;
		}

		if ($skip) {
			for (my $ii=$lastID+1; $ii<$id; $ii++){
				push (@dsspseq, "-");
			}
			$skip = 0;
		}

		push (@dsspseq, $aa);
		$lastID = $id;
	}

	close DSSP or die "Couldn't close file $dssp\n";
	return @dsspseq;
}

################################################################
# Name : load_dnss
# Takes: string: dnss file path
# Returns: array: AA sequence found in dnss
################################################################

sub load_dnss {
	my ($dnss) = @_;

	open (DNSS, $dnss) or die "Couldn't open file $dnss\n";
	my @dnsseq;
	my $line = <DNSS>;
	$line = <DNSS>;
	chomp($line);
	@dnsseq = split(//, $line);
	close DNSS or die "Couldn't close file $dnss\n";

	return @dnsseq;
}

sub clean_seqs {
	my (@seq_refs) = @_;

	for (my $ii=0; $ii<@seq_refs; $ii++){
		my @seq = @{ $seq_refs[$ii] };
		@seq = remove_gaps(@seq);
		$seq_refs[$ii] = \@seq;
	}
	return @seq_refs;
}

################################################################
# Name : align_and_print
# Takes: 2 array refs: names and AA sequence refs
# Returns: prints aligned sequences
################################################################

sub align_and_print {
	my (@seqs) = @_;
	
	my @aligned = align_seqs(@seqs);
	print_comparison(@aligned);
}

################################################################
# Name : align_and_print_named
# Takes: 2 array refs: names and AA sequence refs
# Returns: prints aligned sequences with names
################################################################

sub align_and_print_named {
	my ($names_ref, $seqs_ref) = @_;
	my @seqs = @{ $seqs_ref };
	
	my @aligned = align_seqs(@seqs);
	print_comparison_named($names_ref, \@aligned);
}

################################################################
# Name : align_seqs_absolute
# Takes: 2 array refs: AA sequences
# Returns: 2 array refs: new aligned sequences using slow but sure method
################################################################

sub align_seqs {
        my (@seq_refs) = @_;
	return 0 if (@seq_refs != 2);	

	my @align_refs = get_align_arrays(@seq_refs);
	my @newseqs = apply_alignment(\@seq_refs, \@align_refs);
	
	return @newseqs;
}

################################################################
# Name : get_align_arrays_absolute
# Takes: 2 array refs: AA sequences
# Returns: 2 array refs: alignment of sequences using slow but sure method
################################################################

sub get_align_arrays_absolute {
        my (@seq_refs) = @_;

	@seq_refs = clean_seqs(@seq_refs);
        my ($map_ref) = MAT_make_mapping(@seq_refs);
        my (@align_refs) = ALIGN_make_arrays($map_ref, @seq_refs);

	return @align_refs;
}

################################################################
# Name : align_seqs
# Takes: 2 array refs: AA sequences
# Returns: 2 array refs: new aligned sequences
################################################################

sub align_seqs {
        my (@seq_refs) = @_;
	return 0 if (@seq_refs != 2);	

	my @align_refs = get_align_arrays(@seq_refs);
	my @newseqs = apply_alignment(\@seq_refs, \@align_refs);
	
	return @newseqs;
}

################################################################
# Name : get_align_arrays
# Takes: 2 array refs: AA sequences
# Returns: 2 array refs: alignment of sequences
################################################################

sub get_align_arrays {
        my (@seq_refs) = @_;
	my @trans_refs;

	foreach (@seq_refs){
		my @seq = @{ $_ };
		@seq = transform(@seq);
		push (@trans_refs, \@seq);
	}

        my ($map_ref) = MAT_quick_mapping(@trans_refs);
        my (@align_refs) = ALIGN_make_arrays($map_ref, @trans_refs);

	return @align_refs;
}

sub quick_print {
	my ($map_ref, $seq_refref) = @_;
	my ($seq_ref1, $seq_ref2) = @{ $seq_refref};
	my @seq1 = @{ $seq_ref1 };
	my @seq2 = @{ $seq_ref2 };
	my %map = %{ $map_ref };

	my @newseq1;
	my @newseq2;

	my @keys = sort {$a <=> $b} keys (%map);
	foreach (@keys) {
		push (@newseq1, $seq1[$_]);
		push (@newseq2, $seq2[$map{$_}]);
	}

	print_comparison (\@newseq1, \@newseq2);
}

################################################################
# Name : MAT_quick_mapping
# Takes: 2 array refs: AA sequences
# Returns: hash containing mapping from seq1 to seq2
# Note : Does recursive greedy alignments until seqs are < 10 residues
################################################################

sub MAT_quick_mapping {
	my ($seq1_ref, $seq2_ref, $beg, $end) = @_;
	my %mapping;
	my $beg_const = .05;
	my $end_const = .05;

	my @seq1 = @{ $seq1_ref };
	my @seq2 = @{ $seq2_ref };
	my $len1 = @seq1;
	my $len2 = @seq2;

	if ( ($len1,$len2)[$len1>$len2] < 11){
		my $map_ref = {};
		return $map_ref if ($len1==0 || $len2==0);
		$map_ref = MAT_make_mapping($seq1_ref, $seq2_ref, $beg, $end);
		return $map_ref;
	}
		
	#Here the iteration of greedy alignment occurs 
	my ($index1, $index2, $length) = greedy_align($seq1_ref, $seq2_ref);
	if ($length == 0){
		return MAT_make_mapping($seq1_ref, $seq2_ref, $beg, $end);
	}

	for (my $ii=$index1; $ii<$index1+$length; $ii++){
		$mapping{$ii} = $index2-$index1+$ii;
	}

	if ($index1>0 && $index2>0){
		my @substr1 = @seq1[0..$index1-1];
		my @substr2 = @seq2[0..$index2-1];
		my $map_ref = MAT_quick_mapping(\@substr1, \@substr2, 0, $end_const);

		my %beg_map = %{ $map_ref };

		foreach (keys (%beg_map)){
			$mapping{$_} = $beg_map{$_};
		}
	}

	if ($index1<$len1 && $index2<$len2){
		my $start1 = $index1+$length;
		my $start2 = $index2+$length;
		my @substr1 = @seq1[$start1..$len1-1];
		my @substr2 = @seq2[$start2..$len2-1];
		my $map_ref = MAT_quick_mapping(\@substr1, \@substr2, $beg_const);

		my %end_map = %{ $map_ref };
		
		foreach (keys (%end_map)){
			$mapping{$start1+$_} = $start2 + $end_map{$_};
			$mapping{$start1+$_} = -1 if ($end_map{$_} == -1);
		}
	}

	return \%mapping;
}

################################################################
# Name : greedy_align
# Takes: 2 array refs: AA sequences
# Returns: start indices of alignment for each seq and the length
################################################################

sub greedy_align {
	my ($seq1_ref, $seq2_ref) = @_;

	my @seq1 = @{ $seq1_ref };
	my @seq2 = @{ $seq2_ref };
	my $len1 = @seq1;
	my $len2 = @seq2;
	my $minlen = ($len1, $len2)[$len1 > $len2];

	my $bestscore = 0;
	my ($bestibot, $bestjbot);
	
	my $done=0;
	my $imax = $len1;
	my $imin = $len1-1;
	my $jmin = 0;
	my $jmax = 1;
	until ($done){
		my @subseq1 = @seq1[$imin..$imax];
		my @subseq2 = @seq2[$jmin..$jmax];
		my $len = $imax - $imin;

		my ($starti, $startj, $highi, $highj);
		my $highscore = 0;
		my $starti = -1;
		my $startj = -1;
		my $score = 0;
		for (my $kk=0; $kk<$len; $kk++){
			if ($subseq1[$kk] eq $subseq2[$kk]){
				$starti = $imin + $kk if ($starti == -1);
				$startj = $jmin + $kk if ($startj == -1);
				$score++;
				if ($score > $bestscore){
					$bestscore = $score;
					$bestibot = $starti;
					$bestjbot = $startj;
				}
			}
			else {
				$starti = -1;
				$startj = -1;
				$score = 0;
			}
		}
		last if ($imax == 1);		

		if ($imin > 0){
			$imin--;
			$imax-- if ($imax-$imin > $minlen);
		}
		else {
			$imax-- if ($jmax == $len2);
		}

		if ($jmax < $len2){
			$jmax++;
			$jmin++ if ($jmax-$jmin > $minlen);
		}
		else {
			$jmin++ if ($imax-$imin < $minlen);
		}
	}

	return ($bestibot, $bestjbot, $bestscore);
}


################################################################
# Name : MAT_make_mapping
# Takes: 2 array refs: AA sequences
#	 beg and end: Scores to weight matches starting at beginning or
#		ending at end, respectively
# Returns: hash containing mapping from seq1 to seq2
# Note : Uses Needleman-Wunsh Algorithm (I think)
################################################################

sub MAT_make_mapping {
	my ($seq1_ref, $seq2_ref, $beg, $end) = @_;

	my @seq1 = @{ $seq1_ref };
	my @seq2 = @{ $seq2_ref };
	my $len1 = @seq1;
	my $len2 = @seq2;

	my @Matrix;
	for (my $ii=0; $ii<$len1+1; $ii++){
		my @temp;
		for (my $jj=0; $jj<$len2+1; $jj++){
			$temp[$jj] = 0;
		}
		$Matrix[$ii] = \@temp;
	}

	my $max = 0;
	my $ipos = 0;
	my $jpos = 0;
	for (my $ii=1; $ii<$len1+1; $ii++){
		for (my $jj=1; $jj<$len2+1; $jj++){
			my $score = 0;
			$score += 1 if ($seq1[$ii-1] eq $seq2[$jj-1]);
			$score += $end if ($len1-$ii == $len2-$jj);
			$score += $beg if ($ii == $jj);

			my $diag = $Matrix[$ii-1][$jj-1] + $score;

			my ($delmax, $index) = max_of_col(\@Matrix, $ii-2, $jj-1);
			my $del = $delmax + $score - (.25 + .01*($ii-1-$index));

			my ($insmax, $index) = max_of_row(\@Matrix, $ii-1, $jj-2);
			my $ins = $insmax + $score - (.25 + .01*($jj-1-$index));

			my ($value) = max_of_array($diag, $del, $ins);
			$Matrix[$ii][$jj] = $value;
			if ($value > $max){
				$max = $value;
				$ipos = $ii;
				$jpos = $jj;
			}
		}
	}
#####
##print_csv(\@Matrix, \@seq1, \@seq2);
#####
	my %used1;
	for (my $ii=0; $ii<$len1; $ii++){
		$used1{$ii} = 0;
	}
	my %used2;
	for (my $jj=0; $jj<$len2; $jj++){
		$used2{$jj} = 0;
	}

	my %mapping;
	my $max = $Matrix[$ipos][$jpos];
	my $firsti = $ipos-1;
	my $lasti = $ipos-1; 
	my $firstj = $jpos-1;
	my $lastj = $jpos-1;
	while($max > 0){
		$mapping{$ipos-1} = $jpos-1;
		$used1{$ipos-1}++;
		$used2{$jpos-1}++;
		$firsti = $ipos-1;
		$firstj = $jpos-1;

		($max, $ipos, $jpos) = get_next_max(\@Matrix, $ipos, $jpos);
	}

	while ($firsti > 0 && $firstj > 0){
		$firsti--;
		$firstj--;
		$mapping{$firsti} = $firstj;
		$used1{$firsti}++;
		$used2{$firstj}++;
	}
	
	foreach (keys %used1){
		if ($used1{$_} == 0){
			$mapping{$_} = -1;
			$used1{$_}++;
		}
	}
	return (\%mapping);
}

################################################################
# Name : print_csv
# Takes: 2D array of an alignment score matrix
# Returns: prints csv file (test/Matrix.csv)
################################################################

sub print_csv {
	my ($M, $seq1_ref, $seq2_ref) = @_;
	my @Matrix = @{ $M };
	my @seq1 = @{ $seq1_ref };
	my @seq2 = @{ $seq2_ref };
	my $csvline = ",";

	my $outfile = "test/Matrix.csv";
	open (OUT, ">$outfile") or die "Couldn't open file $outfile\n";
	print "Writing file $outfile\n";

	for (my $jj=0; $jj<@seq2; $jj++){
		$csvline .= "$seq2[$jj],";
	}
	print OUT "$csvline\n";

	for (my $ii=1; $ii<=@seq1; $ii++){
		$csvline = "$seq1[$ii-1],";
		for (my $jj=1; $jj<=@seq2; $jj++){
			$csvline .= "$Matrix[$ii][$jj],";
		}
		print OUT "$csvline\n";
	}

	close OUT or die "Couldn't close file $outfile\n";
}

################################################################
# Name : get_next_max
# Takes: 2D array of an alignment score matrix and i, j coordinates of current max
# Returns: max value and i, j coordinates of next maximum
################################################################

sub get_next_max{
	my ($mat_ref, $ipos, $jpos) = @_;
	my @Matrix = @{ $mat_ref };

	my $max = $Matrix[$ipos-1][$jpos-1];
	my $besti = $ipos-1;
	my $bestj = $jpos-1;
	if ($max == $Matrix[$ipos][$jpos]-1){
		return ($max, $besti, $bestj);
	}

	my ($delmax, $delind) = max_of_col(\@Matrix, $ipos-2, $jpos-1);
	my ($insmax, $insind) = max_of_row(\@Matrix, $ipos-1, $jpos-2);

	my ($max, $index) = max_of_array($max, $delmax, $insmax);
	return ($max, $besti, $bestj) if ($index == 0);
	return ($max, $delind, $jpos-1) if ($index == 1);
	return ($max, $ipos-1, $insind);
}

################################################################
# Name : max_of_col
# Takes: 2D array of score matrix, maximum row, column index
# Returns: maximum value found in that column in rows up to indicated maximum
################################################################

sub max_of_col {
	my ($mat_ref, $imax, $jj) = @_;
	my @Matrix = @{ $mat_ref };

	return 0 if ($imax < 0 || $jj<0);

	my @values;
	for (my $ii=0; $ii<=$imax; $ii++){
		push (@values, $Matrix[$ii][$jj]);
	}
	my ($maxval, $index) = max_of_array(@values);
	
	return ($maxval, $index);
}

################################################################
# Name : max_of_row
# Takes: 2D array of score matrix, row index, maximum column
# Returns: maximum value found in that row in columns up to indicated maximum
################################################################

sub max_of_row {
	my ($mat_ref, $ii, $jmax) = @_;
	my @Matrix = @{ $mat_ref };

	return 0 if ($jmax < 0 || $ii < 0);

	my @values;
	for (my $jj=0; $jj<=$jmax; $jj++){
		push (@values, $Matrix[$ii][$jj]);
	}
	my ($maxval, $index) = max_of_array(@values);
	
	return ($maxval, $index);
}

################################################################
# Name : max_of_array
# Takes: an array of values
# Returns: maximum value of the array as well as its index
################################################################

sub max_of_array {
	my @values = @_;

	my $max = 0;
	my $index = 0;
	for (my $ii=0; $ii<@values; $ii++){
		if ($values[$ii] > $max){
			$max = $values[$ii];
			$index = $ii;
		}
	}
	return ($max, $index);
}

################################################################
# Name : ALIGN_make_arrays
# Takes : hash_ref : mapping of sequence1 to sequence2
# 	: 2 array_refs : AA sequences
# Returns: 2 array refs containing alignment
################################################################

sub ALIGN_make_arrays {
	my ($mapping_ref, $seq1_ref, $seq2_ref) = @_;
	my %mapping = %{ $mapping_ref };
	my @seq1 = @{ $seq1_ref };
	my @seq2 = @{ $seq2_ref };

	my (@align1, @align2);

	#This bit is to make the alignment arrays
	my @keys = keys (%mapping);
	my @order1 = sort {$a <=> $b} @keys;
	my @used2;
	my $Aindex = 0;
	my $done = 0;

	#This loop assigns a value to each position of sequence1
	for (my $ii=0; $ii<@seq1; $ii++){
		my $pos = $order1[$Aindex];
		if ($pos != $ii){
			$align1[$ii] = $ii;
			$align2[$ii] = -1;
		}
		else {
			$align1[$ii] = $pos;
			$align2[$ii] = $mapping{$pos};
			$Aindex++;
			$used2[$mapping{$pos}] = 1 unless ($mapping{$pos} == -1);
		}
	}

# Old loop for above function. Was replaced when found to not assign align vals to
# seq1 right overhangs
#	foreach my $pos (@order1){
#		if ($pos != $Aindex) {
#			for ($Aindex..($pos-1)){
#				$align1[$Aindex] = $Aindex;
#				$align2[$Aindex] = -1;
#				$Aindex++;
#			}
#		}
#		$align1[$Aindex] = $pos;
#		$align2[$Aindex] = $mapping{$pos};
#		$Aindex++;	
#		$used2[$mapping{$pos}] = 1 unless ($mapping{$pos} == -1);
#	}

	#The following makes the beginning of seq2 be used if there is an overhang
	my $startpos = $align2[0];
	if ($startpos != 0){
		for (my $ii=$startpos-1; $ii>=0; $ii--){
			my $seq1length = @seq1;
			unshift(@align1, -1);
			unshift(@align2, $ii);
			$used2[$ii] = 1;
		}
	}
	
	#Determines which of seq2 are still unused
	my @unused2;
	for (my $ii=0; $ii<@seq2; $ii++){
		if ($used2[$ii] == 1){
			next;
		}
		else {
			push (@unused2, $ii);
		}
	}

	#Ensures that all  members of seq2 are used, inserting dashes into seq1
	my $Aindex = 0;
	foreach my $addition (@unused2){
		my $seq1length = @seq1;
		my $seq2length = @seq2;
		while (($align2[$Aindex] < $addition - 1)){ 
			$Aindex++;
			#last if ($Aindex >= @align2);
		}
		splice(@align1, $Aindex+1, 0, -1);
		splice(@align2, $Aindex+1, 0, $addition);
		$Aindex++;
	}

	my @seq_refs = (\@seq1, \@seq2);
	my @align_refs = (\@align1, \@align2);
	(@align_refs) = ALIGN_remove_gaps(\@seq_refs, \@align_refs);

	return (@align_refs);
}

################################################################
# Name : ALIGN_remove_gaps
# Takes: 2 array refs: array of sequences and array of alignment
# Returns: array of array refs of new alignments
################################################################

sub ALIGN_remove_gaps {
	my ($seq_ref_ref, $align_ref_ref) = @_;
	my @seq_refs = @{ $seq_ref_ref };
	my @align_refs = @{ $align_ref_ref };

	my @seq1 = @{ $seq_refs[0] };
	my @seq2 = @{ $seq_refs[1] };
	my @align1 = @{ $align_refs[0] };
	my @align2 = @{ $align_refs[1] };

	$seq1[@seq1] = "-" if ($seq1[$#seq1] ne "-");
	$seq2[@seq2] = "-" if ($seq2[$#seq2] ne "-");

	#Removes segments of mutual non-residues
	for (my $index=0; $index<@align1; $index++){
		while ($seq1[$align1[$index]] eq "-" && $seq2[$align2[$index]] eq "-"){
			splice(@align1, $index, 1);
			splice(@align2, $index, 1);
		}
	}	

	#Removes hidden mutual non-residues
	for (my $index=0; $index<@align1; $index++){
		if ($seq1[$align1[$index]] eq "-"){
			if ($seq2[$align2[abs($index-1)]] eq "-"){
				splice(@align1, $index, 1);
				splice(@align2, $index-1, 1);
			}
			elsif ($seq2[$align2[$index+1]] eq "-"){
				splice(@align1, $index, 1);
				splice(@align2, $index+1, 1);
			}
		}	
	}
	return (\@align1, \@align2);
}

################################################################
# Name : remove_gaps
# Takes: array of sequence
# Returns: array of sequence without gaps
################################################################

sub remove_gaps {
	my (@seq) = @_;
	my @newseq;

	foreach (@seq){
		push (@newseq, $_) unless ($_ eq "-");
	}
	return @newseq;
}


################################################################
# Name : apply_alignment
# Takes: 2 seq array refs and 2 alignment array refs
# Returns: 2 array refs containing aligned sequence
################################################################

sub apply_alignment {
	my ($seq_ref_ref, $align_ref_ref) = @_;
	my @seq_refs = @{ $seq_ref_ref };
	my @align_refs = @{ $align_ref_ref };

	my @seq1 = @{ $seq_refs[0] };
	my @seq2 = @{ $seq_refs[1] };
	my @align1 = @{ $align_refs[0] };
	my @align2 = @{ $align_refs[1] };
	
	my $a1len = @align1;
	my $a2len = @align2;

	my $seq1len = @seq1;
	my $seq2len = @seq2;

	$seq1[@seq1] = "-" if ($seq1[$#seq1] ne "-");
	$seq2[@seq2] = "-" if ($seq2[$#seq2] ne "-");

	my (@newseq1, @newseq2);
	foreach my $index (@align1){
		push (@newseq1, $seq1[$index]);
	}
	foreach my $index (@align2){	
		push (@newseq2, $seq2[$index]);
	}
	
	return (\@newseq1, \@newseq2);
}

################################################################
# Name : align_score
# Takes: 2 array refs: AA sequences
# Returns: float: percent of identically aligned residues
################################################################

sub align_score {
	my ($seq1_ref, $seq2_ref) = @_;
	my @seq1 = @{ $seq1_ref };
	my @seq2 = @{ $seq2_ref };

	my $match = 0;
	my $min = ($#seq1, $#seq2)[$#seq1 > $#seq2];

	return 0 if ($min <= 0);
	$min++;
	my $count;

	for (my $ii=0; $ii<$min; $ii++){
		if ($seq1[$ii] eq $seq2[$ii]){
			$match++;
			$count++;
		}
		elsif ($seq1[$ii] eq "X" || $seq2[$ii] eq "X") {
			if ($seq1[$ii] eq "-" || $seq2[$ii] eq "-") {
				$count++;
			}
		}
		else {
			$count++;
		}
	}
	my $score = $match/$count;
	return $score;
}

################################################################
# Name : matched_score
# Takes: 2 array refs: AA sequences
# Returns: float: align_score bounded by shorter AA segment
################################################################

sub matched_score {
	my ($seq1_ref, $seq2_ref) = @_;
	my @seq1 = @{ $seq1_ref };
	my @seq2 = @{ $seq2_ref };

	my $first1 = 0;
	my $last1 = 0;
	for (my $ii=0; $ii<@seq1; $ii++){
		if ($seq1[$ii] ne "-"){
			$last1 = $ii;
			if ($first1==0){
				$first1 = $ii;
			}
		}
	}

	my $first2 = 0;
	my $last2 = 0;
	for (my $ii=0; $ii<@seq2; $ii++){
		if ($seq2[$ii] ne "-"){
			$last2 = $ii;
			if ($first2==0){
				$first2 = $ii;
			}
		}
	}

	my $match = 0;
	my $first = ($first1, $first2)[$first1 < $first2];
	my $last = ($last1, $last2)[$last1 > $last2];

	my $len = $last - $first;
	return 0 if ($len <= 0);

	my @matched1 = @seq1[$first..$last];
	my @matched2 = @seq2[$first..$last];

	return align_score(\@matched1, \@matched2);
}

################################################################
# Name : transform
# Takes: array: sequence
# Returns: sequence after transforming "-" <=> "_"
################################################################

sub transform {
	my (@seq) = @_;

	for (my $ii=0; $ii<@seq; $ii++){
		if ($seq[$ii] eq "_"){ $seq[$ii] = "-"; }
		elsif ($seq[$ii] eq "-"){ $seq[$ii] = "_"; }
	}
	return @seq;
}	

################################################################
# Name : print_comparison
# Takes: 2D array ref: one or multiple AA sequences
# Returns: Prints sequences to stdout in refined format
################################################################

sub print_comparison {
	my (@sequences) = @_;

	my $past = 0;
	my $start = 0;
	while ($past < @sequences){
		print "         ";
		
		my $maxindex = 0;
		foreach (@sequences){
			my @seq = @{ $_ };
			$maxindex = @seq if (@seq > $maxindex);
		}

		my $count = 1;
		for (my $ii=$start; $ii<$start+50; $ii++){
			last if ($ii == $maxindex);
			$count-=10 while ($count >= 10);
			print $count;
			$count++;
		}

		print "\n";

		for (my $ii=0; $ii<@sequences; $ii++) {
			print "         ";
			my @thisseq = @{ $sequences[$ii] };

			for (my $jj=$start; $jj<$start+100; $jj++){
				print $thisseq[$jj];
				if ($jj >= $#thisseq){
					$past++;
					last;
				}
			}
			print "\n";
		}
		print "\n" unless ($past == @sequences);
		$start += 100;
	}
}

################################################################
# Name : print_comparison_named
# Takes: 1: array ref: names of the corresponding AA sequence
#	 2: 2D array ref: one or multiple AA sequences
# Returns: Prints sequences to stdout in refined format
################################################################

sub print_comparison_named {
	my ($name_ref, $sequence_ref) = @_;
	my @names = @{ $name_ref };
	my @sequences = @{ $sequence_ref };

	my $past = 0;
	my $start = 0;
	while ($past < @sequences){
		print "         ";
		my $maxindex = 0;
		foreach (@sequences){
			my @seq = @{ $_ };
			$maxindex = @seq if (@seq > $maxindex);
		}

		my $count = 1;
		for (my $ii=$start; $ii<$start+100; $ii++){
			last if ($ii == $maxindex);
			$count-=10 while ($count >= 10);
			print $count;
			$count++;
		}

		print "\n";

		for (my $ii=0; $ii<@sequences; $ii++) {

			print $names[$ii];
			for (my $jj=length($names[$ii]); $jj<9; $jj++){
				print " ";
			}
			my @thisseq = @{ $sequences[$ii] };

			for (my $jj=$start; $jj<$start+100; $jj++){
				print $thisseq[$jj];
				if ($jj >= $#thisseq){
					$past++;
					last;
				}
			}
			print "\n";
		}
		print "\n" unless ($past == @sequences);
		$start += 100;
	}
}

################################################################
# Name : print_array
# Takes: array of array refs: AA sequences
# Returns: Prints arrays
################################################################

sub print_array {
	my (@array_refs) = @_;

	print ">Alignment Arrays:\n";
	foreach my $ref (@array_refs){
		my @array = @{ $ref };
		print "@array\n";
	}
}

################################################################
# Name : print_separator
# Takes: -
# Returns: Prints a dashed line separator
################################################################

sub print_separator {
	print "\n-----------------------------------------------------------\n\n";
}

################################################################
# Name : print_double_separator
# Takes: -
# Returns: Prints two dashed line separators
################################################################

sub print_double_separator {
	print "\n-----------------------------------------------------------\n";
	print "-----------------------------------------------------------\n\n";
}


1;
