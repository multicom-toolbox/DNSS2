package PDBUtils;

##### List of subs in PDBUtils #####
# get_res_name
# get_res_code
# get_atom
# get_res
# get_atom_num
# get_res_num
# get_coords
####################################

sub get_res_name{
	my ($res) = @_;
	%RES_NAMES = qw/A ALA C CYS D ASP E GLU F PHE G GLY H HIS
			I ILE K LYS L LEU M MET N ASN P PRO Q GLN
			R ARG S SER T TYR V VAL W TRP Y TYR/;
	return $RES_NAMES{$res};
}

sub get_res_code{
	my ($res) = @_;
	%RES_CODES = qw/ALA A CYS C ASP D GLU E PHE F GLY G HIS H 
			ILE I LYS K LEU L MET M ASN N PRO P GLN Q 
			ARG R SER S THR T VAL V TRP W TYR Y/;
	return $RES_CODES{$res};
}

sub get_atom {
	my ($pdbline) = @_;

	my $endline = substr($pdbline, 12);
	my ($atom) = $endline =~ /\s*(\S+)/;

	return $atom;
}

sub get_res {
	my ($pdbline) = @_;

	my $endline = substr($pdbline, 17);
	my ($res) = $endline =~ /\s*(\w+)/;

	return $res;
}

sub get_atom_num {
	my ($pdbline) = @_;

	my $endline = substr($pdbline, 6);
	my ($anum) = $endline =~ /\s*(\d+)/;

	return $anum;
}

sub get_res_num {
	my ($pdbline) = @_;

	my $endline = substr($pdbline, 20);
	my ($resnum) = $endline =~ /\s*(\d+)/;

	return $resnum;
}

sub get_coords {
	my ($pdbline) = @_;

	my $endline = substr($pdbline, 30);
	my ($xx, $yy, $zz) = $endline =~ /\s*(\S+)\s*(\S+)\s*(\S+)/;

	my @coords = ($xx, $yy, $zz);

	return @coords;
}

1;
