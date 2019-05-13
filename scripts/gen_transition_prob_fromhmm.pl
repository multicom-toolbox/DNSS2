#!/usr/bin/perl -w

$numArgs = @ARGV;
if($numArgs != 2)
{   
	print "the number of parameters is not correct!\n";
	exit(1);
}

$fastafile		= "$ARGV[0]";
$hmmfile		= "$ARGV[1]";
#$fastafile		= "$fastadir/$id.fasta";
#$hmmfile		= "$hmmdir/$id.hmm";
if(!(-e $fastafile) or !(-e $hmmfile))
{
  print "Failed to find $hmmfile or $fastafile\n";
}

open(OUT1,">${hmmfile}2freq") || die "Failed to open dir ${hmmfile}2freq\n";
open(OUT2,">${hmmfile}2prob") || die "Failed to open dir ${hmmfile}2prob\n";
open(IN,$fastafile) || die "Failed to open dir $fastafile\n";
@content = <IN>;
close(IN);
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




$start=0;

$c=0;

open(IN,$hmmfile) || die "Failed to open dir $hmmfile\n";
@content = <IN>;
close(IN);

foreach $line (@content)
{
  chomp $line;
  $line =~ s/^\s+|\s+$//g;
  if($line eq '//')
  {
    next;
  }
  if(substr($line,0,4) eq 'NULL')
  {
    $start=1;
  }
  
  if($start==0)
  {
    next;
  }
  if(substr($line,0,3) eq 'HMM')
  {
    $start=2;
  }
  if($start==2)
  {
    $c++;
  }
  
  
  if($c % 3 == 1)
  {
    
    if(substr($line,0,3) eq 'HMM')
    {
      @tmp = split(/\s+/,$line);
      shift @tmp;
      if(@tmp != 20)
      {
        print "Warning: only have ".@tmp." elements for emission probability, should be 20\n";
      }
      print OUT1 join("\t",@tmp);
      print OUT2 join("\t",@tmp);
    }else{
      @tmp = split(/\s+/,$line);
      pop @tmp;
      if(@tmp != 22)
      {
        print "Warning: only have ".@tmp." elements for emission probability, should be 22\n";
      }
      
      $amino_acid = $tmp[0];
      $position = $tmp[1];
      $amino_acid_in_seq = substr($sequence,$position-1,1);
      if( $amino_acid_in_seq ne $amino_acid)
      {
        die "Error, the amino acid in position $position is $amino_acid, but should be $amino_acid_in_seq\n";
      }
      
      $prob_array="";
      for($i=2;$i<22;$i++)
      {
        $freq = $tmp[$i]; #freq=  -1000 × log2 prob  https://github.com/soedinglab/hh-suite/blob/master/hhsuite-userguide.pdf
        if($freq eq '*')
        {
          $prob = 0; # 
        }else{
          $prob = sprintf("%.6f",2**(-$freq/1000)); # prob = 2^(-freq/1000)
        }
        $prob_array .= "\t".$prob;
      }
      $prob_array =~ s/^\s+|\s+$//g;
      print OUT1 join("\t",@tmp);
      print OUT2 "$amino_acid\t$position\t".$prob_array; 
    }
  }
  if($c % 3 == 2)
  {
    if(index($line,'Neff')>0)
    {
      @tmp = split(/\s+/,$line);
      pop @tmp;
      pop @tmp;
      pop @tmp;
      if(@tmp != 7)
      {
        print "Warning: only have ".@tmp." elements for transition probability, should be 7\n";
      }
      print OUT1 " ".join("\t",@tmp)."\n";
      print OUT2 " ".join("\t",@tmp)."\n";
    }else{
      @tmp = split(/\s+/,$line);
      pop @tmp;
      pop @tmp;
      pop @tmp;
      if(@tmp != 7)
      {
        print "Warning: only have ".@tmp." elements for transition probability, should be 7\n";
      }
      $prob_array="";
      for($i=0;$i<7;$i++)
      {
        $freq = $tmp[$i]; #freq=  -1000 × log2 prob  https://github.com/soedinglab/hh-suite/blob/master/hhsuite-userguide.pdf
        if($freq eq '*')
        {
          $prob = 0; # 
        }else{
          $prob = sprintf("%.6f",2**(-$freq/1000)); # prob = 2^(-freq/1000)
        }
        $prob_array .= "\t".$prob;
      }
      print OUT1 " ".join("\t",@tmp)."\n";
      print OUT2 " ".$prob_array."\n";
    }
    
  }
}
close OUT1;
close OUT2;
