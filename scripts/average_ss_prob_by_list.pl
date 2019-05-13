#!/usr/bin/env perl

#### written by Jie Hou to plot training figures
### made  2017/11/04

# perl /storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/paper_dataset/State1_Nets_On_OldFea/lib/average_ss_prob_by_list.pl   /storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/datasets/new20181005/adj_dncon-train.lst /storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/paper_dataset/State1_Nets_On_OldFea/emsemble/model_test.list  /storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/paper_dataset/State1_Nets_On_OldFea/outputs/ensemble_out/old/ /storage/htc/bdm/Collaboration/Zhiye/SSP/DNSS2/paper_dataset/State1_Nets_On_OldFea/outputs/ensemble_out/old/avg_out

$num = @ARGV;

if($num != 5)
{
	die "The number of parameter is not correct!\n";
}

$testfile = $ARGV[0]; #
$netfile = $ARGV[1]; #
$prediction_dir = $ARGV[2]; #
$out_dir = $ARGV[3]; #
$tag = $ARGV[4]; #

open(IN,"$netfile") || die "Failed to open file $netfile\n";

@network_array=();
while(<IN>)
{
  $netname=$_;
  chomp $netname;
  push @network_array,$netname;
}
close IN;

$net_num = @network_array;

open(IN,"$testfile") || die "Failed to open file $testfile\n";

while(<IN>)
{
  $target=$_;
  chomp $target;
  
  #print "$target\n";
  %prob_avg=();
  foreach $netname (@network_array)
  {
    $prob_file="$prediction_dir/$netname/test_prediction/$target.prob";
    #print "$prob_file\n";
    if(!(-e $prob_file))
    {
      die "Failed to find $prob_file\n";
    }
    open(TMP,"$prob_file") || die "Failed to find $prob_file\n";
    $c=0;
    while(<TMP>)
    {
      $li=$_;
      chomp $li;
      $c++;
      @tmp = split(/\s++/,$li);
      $h_prob = $tmp[0];
      $e_prob = $tmp[1];
      $c_prob = $tmp[2];
      if(exists($prob_avg{$c}))
      {
        @tmp2 =split(/\s/,$prob_avg{$c});
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
  $out_file1="$out_dir/$target.prob";
  $out_file2="$out_dir/$target.pred";
  open(OUT1,">$out_file1") || die "Failed to find $out_file1\n";
  open(OUT2,">$out_file2") || die "Failed to find $out_file2\n";
  foreach  $indx (sort { $a <=> $b } keys %prob_avg) {
    @tmp2 =split(/\s/,$prob_avg{$indx});
    $max_index = maxindex(\@tmp2);
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
}
close IN;

`touch $out_dir/${tag}_average.done`;

sub maxindex {
  my( $aref, $idx_max ) = ( shift, 0 );
  $aref->[$idx_max] > $aref->[$_] or $idx_max = $_ for 1 .. $#{$aref};
  return $idx_max;
}
