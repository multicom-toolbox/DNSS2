#!/usr/bin/perl -w
use Cwd;
use Cwd 'abs_path';

################ provide the path for databases
$uniref90db = '/faculty/jhou4/tools/DNSS2/database/uniref90/uniref90.fasta';
$uniclust30db = '/faculty/jhou4/tools/DNSS2/database/uniclust30_2017_10/uniclust30_2017_10/uniclust30_2017_10';



################Don't Change the code below##############

if (! -e "${uniclust30db}_a3m_db")
{
	die "can't find uniclust database <${uniclust30db}_a3m_db>.\n";
}

if (! -e "$uniref90db")
{
	die "can't find nr sequence database <$uniref90db>.\n";
}

$install_dir = getcwd;
$install_dir=abs_path($install_dir);

if (! -d $install_dir)
{
	die "can't find installation directory.\n";
}
if ( substr($install_dir, length($install_dir) - 1, 1) ne "/" )
{
	$install_dir .= "/"; 
}


if (prompt_yn("DNSS2 will be installed into <$install_dir> ")){

}else{
	die "The installation is cancelled!\n";
}
print "Start install DNSS2 into <$install_dir>\n\n";


$files		="scripts/PredictSS.pl,lib/DN_SSpred.pm,lib/DN_SSpred2.pm,./run_DNSS2.pl";

@updatelist		=split(/,/,$files);

foreach my $file (@updatelist) {
	$file2update=$install_dir.$file;
	
	$check_log ='$GLOBAL_PATH=';
	open(IN,$file2update) || die "Failed to open file $file2update\n";
	open(OUT,">$file2update.tmp") || die "Failed to open file $file2update.tmp\n";
	while(<IN>)
	{
		$line = $_;
		chomp $line;

		if(index($line,$check_log)>=0)
		{
			print "### Setting ".$file2update."\n";
			print "\t--- Current ".$line."\n";
			print "\t--- Change to ".substr($line,0,index($line, '=')+1)." \'".$install_dir."\';\n\n\n";
			print OUT substr($line,0,index($line, '=')+1)."\'".$install_dir."\';\n";
		}else{
			print OUT $line."\n";
		}
	}
	close IN;
	close OUT;
	system("mv $file2update.tmp $file2update");
	system("chmod 755  $file2update");


}


$files		="scripts/generate-pssm_withDB.sh,scripts/gen-pssm-less-stringent_withDB.sh,scripts/generate-hmm_withDB.sh";

@updatelist		=split(/,/,$files);

foreach my $file (@updatelist) {
	$file2update=$install_dir.$file;
	
	$check_log ='GLOBAL_PATH=';
	open(IN,$file2update) || die "Failed to open file $file2update\n";
	open(OUT,">$file2update.tmp") || die "Failed to open file $file2update.tmp\n";
	while(<IN>)
	{
		$line = $_;
		chomp $line;

		if(index($line,$check_log)>=0)
		{
			print "### Setting ".$file2update."\n";
			print "\t--- Current ".$line."\n";
			print "\t--- Change to ".substr($line,0,index($line, '=')+1)." \'".$install_dir."\';\n\n\n";
			print OUT substr($line,0,index($line, '=')+1)."\'".$install_dir."\';\n";
		}else{
			print OUT $line."\n";
		}
	}
	close IN;
	close OUT;
	system("mv $file2update.tmp $file2update");
	system("chmod 755  $file2update");
}


###### update database path
$files		="scripts/PredictSS.pl,./run_DNSS2.pl";

@updatelist		=split(/,/,$files);

foreach my $file (@updatelist) {
	$file2update=$install_dir.$file;
	
	$check_log1 ='$uniclust30db=';
	$check_log2 ='$uniref90db=';
	open(IN,$file2update) || die "Failed to open file $file2update\n";
	open(OUT,">$file2update.tmp") || die "Failed to open file $file2update.tmp\n";
	while(<IN>)
	{
		$line = $_;
		chomp $line;

		if(index($line,$check_log1)>=0)
		{
			print "### Setting ".$file2update."\n";
			print "\t--- Current ".$line."\n";
			print "\t--- Change to ".substr($line,0,index($line, '=')+1)." \'".$uniclust30db."\';\n\n\n";
			print OUT substr($line,0,index($line, '=')+1)."\'".$uniclust30db."\';\n";
		}elsif(index($line,$check_log2)>=0)
		{
			print "### Setting ".$file2update."\n";
			print "\t--- Current ".$line."\n";
			print "\t--- Change to ".substr($line,0,index($line, '=')+1)." \'".$uniref90db."\';\n\n\n";
			print OUT substr($line,0,index($line, '=')+1)."\'".$uniref90db."\';\n";
		}else{
			print OUT $line."\n";
		}
	}
	close IN;
	close OUT;
	system("mv $file2update.tmp $file2update");
	system("chmod 755  $file2update");


}

open(IN1,"$install_dir/scripts/sspro.model.default") || die "Failed to open file $install_dir/scripts/sspro.model.default\n";
open(OUT1,">$install_dir/scripts/sspro.model") || die "Failed to open file $install_dir/scripts/sspro.model\n";
while(<IN1>)
{
	$line = $_;
	chomp $line;

	if(index($line,'SOFTWARE_PATH')>=0)
	{
		$line =~ s/SOFTWARE_PATH/$install_dir/g;
		$line =~ s/\/\//\//g;
		print OUT1 $line."\n";
	}else{
		print OUT1 $line."\n";
	}
}
close IN1;
close OUT1;


open(IN1,"$install_dir/lib/DN_SSpred2.pm.default") || die "Failed to open file $install_dir/lib/DN_SSpred2.pm.default\n";
open(OUT1,">$install_dir/lib/DN_SSpred2.pm") || die "Failed to open file $install_dir/lib/DN_SSpred2.pm\n";
while(<IN1>)
{
	$line = $_;
	chomp $line;

	if(index($line,'SOFTWARE_PATH')>=0)
	{
		$line =~ s/SOFTWARE_PATH/$install_dir/g;
		$line =~ s/\/\//\//g;
		print OUT1 $line."\n";
	}else{
		print OUT1 $line."\n";
	}
}
close IN1;
close OUT1;

print "\n!!! Installation finished!\n";

sub prompt_yn {
  my ($query) = @_;
  my $answer = prompt("$query (Y/N): ");
  return lc($answer) eq 'y';
}
sub prompt {
  my ($query) = @_; # take a prompt string as argument
  local $| = 1; # activate autoflush to immediately show the prompt
  print $query;
  chomp(my $answer = <STDIN>);
  return $answer;
}


