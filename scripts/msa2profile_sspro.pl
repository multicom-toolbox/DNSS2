#!/usr/bin/perl -w

use Scalar::Util qw(looks_like_number);
if (@ARGV != 5)
{
	die "need four parameters: ss_predictor, input_sequence, align_dir, output_file\n";
}

$fastafile		= "$ARGV[0]";
$alignfile		= "$ARGV[1]";
$alignments_dir		= "$ARGV[2]";
$ss_predictor = $ARGV[3]; #/msa2profile_sspro
$sspro_model = $ARGV[4]; #/sspro4.1/model/sspro.model

if (! -f $ss_predictor)
{
	die "can't find secondary structure predictor.\n"; 
}
if (! -f $sspro_model)
{
	die "can't find sspro model.\n"; 
}

if (! -d $alignments_dir)
{
	die "the alignment directory doesn't exists.\n";
}

if ( substr($alignments_dir, length($alignments_dir) - 1, 1) ne "/" )
{
	$alignments_dir .= "/"; 
}


#$fastafile		= "$fastadir/$id.fasta";
#$alignfile		= "$alignments_dir/$id.aln";
if(!(-e $fastafile) or !(-e $alignfile))
{
  print "Failed to find $alignfile or $fastafile\n";
}

if (! -f $fastafile)
{
	die "can't find the sequence file.\n"; 
}
open(SEQ_FILE, "$fastafile") || die "can't open sequence file.\n";
@content = <SEQ_FILE>;
close(SEQ_FILE);

$name = shift @content;
chomp $name; 
if(index($name,'|')>0) #1A8L:A|PDBID|CHAIN|SEQUENCE
{
@tmp = split(/\|/,$name);
$name = $tmp[0];
}
if(index($name,' ')>0) #1A8L:A|PDBID|CHAIN|SEQUENCE
{
@tmp = split(/\s/,$name);
$name = $tmp[0];
}
$name =~ s/\s//g;
if(substr($name ,0, 1) eq '>')
{
	$name = substr($name,1);
}
$name =~ s/\:/\-/g;
  
$align_dir = "$alignments_dir/$name";

if(-d "$align_dir")
{
  `rm -rf $align_dir`;
}else{
  `mkdir $align_dir`;
}



if(!(-d $align_dir))
{
  `mkdir $align_dir`;
}
#check if the alignment file exists
$alignment_file = "$align_dir/$name"; 

## new alignment 
open(TMP,"$alignfile") || die "Failed to open file $alignfile\n";
$total_align=0;
while(<TMP>)
{
	$li=$_;
	chomp $li;
	if (looks_like_number($li)) {
		next;
	}else{
		$total_align++;
	}
}
close TMP;
open(TMP,"$alignfile") || die "Failed to open file $alignfile\n";
open(TMPOUT,">$alignment_file") || die "Failed to open file $alignment_file\n";
print TMPOUT "$total_align\n";
$total_align=0;
while(<TMP>)
{
	$li=$_;
	chomp $li;
	if (looks_like_number($li)) {
		next;
	}else{
		print TMPOUT "$li\n";
	}
}
close TMP;
close TMPOUT;

if (! -f $alignment_file)
{
	die "Warning: alignment file: $alignment_file doesn't exist!\n";  
}
#the alignment file name contains space or dot, it will cause trouble. 
if ($name =~ /[\. ]/)
{
	die "sequence name shouldn't contain dot or space.\n";
}

$sequence="";
foreach $seq (@content)
{
  chomp $seq;
  if(substr($seq,0,1) eq '>')
  {
    next;
  }
  $seq =~ s/^\s+|\s+$//g;
  $sequence .=$seq;
}
#make up a pseudo ss string required by SSpro
$sstructure = ""; 
for ($ss_i = 0; $ss_i < length($sequence); $ss_i++)
{
	$sstructure .= "H"; 
}
if ( substr($align_dir, length($align_dir) - 1, 1) ne "/" )
{
	$align_dir .= "/"; 
}

open(TMPFILE, ">$align_dir/$name.tmp") || die "can't create temporary file $align_dir/$name.tmp.\n"; 
print TMPFILE "1 20 3\n$name\n$sequence\n$sstructure"; 
close TMPFILE; 

#predict the secondary structure with 3-line format
print("$ss_predictor $sspro_model $align_dir/$name.tmp $align_dir 0 > $align_dir/$name.msa2profile\n");
system("$ss_predictor $sspro_model $align_dir/$name.tmp $align_dir 0 > $align_dir/$name.msa2profile"); 

open(RES, "$align_dir/$name.msa2profile") || die "can't open prediction result file $align_dir/$name.msa2profile.\n"; 
@org_ss = <RES>;

#$seq_length = @org_ss;
if ( @org_ss != length($sequence) )	
{
	die "sequence length doesn't match with output from SSpro.".length($sequence)." != ".@org_ss."\n"; 
}

 open(OUT,">${alignfile}2profile") || die "Failed to open dir ${alignfile}2profile\n";
print OUT "idx\tpos\tA\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY\n";

foreach $line (@org_ss)
{
	chomp $line;
	@tmp = split(/\s/,$line);
	$id = shift @tmp;
	$res = substr($sequence,$id-1,1);
	print OUT "$id\t$res\t".join("\t",@tmp)."\n";
}
close OUT;

#remove the temporary files

