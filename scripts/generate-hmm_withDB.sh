#!/bin/bash
###############################################################################
# Name %n  : generate-pssm.sh
# Desc %d  : Generate and save a PSSM for a query using psi blast and 
#  suggested parameters from BMC Bioinformatics 2008, 9(Suppl 12):S12 
#  doi:10.1186/1471-2105-9-S12-S12
# Input %i : A fasta file and two output file names, database path
# Output %o: The pssm and a report 
#
# Author: Jie Hou, Jesse Eickholt
# URL: http://merit.oryo.us
# Date: Wed Jan Q25 2012
# Latest Modified Date: 11/02/2016
############################################################################### 
GLOBAL_PATH='/storage/htc/bdm/jh7x3/DNSS2_github/DNSS2/';

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <seq.fasta> <out_fname1(pssm)> <out_fname2(report)> <database>"; 
  exit;
fi

FASTA=$1
targetid=$2
hhblitdb=$3


hhblits3_dir="$GLOBAL_PATH/programs/hhsuite-3.0-beta.3-Linux/"; 

echo "Generate alignments from the nr database using hhblits...";

export HHLIB=$GLOBAL_PATH/programs/hhsuite-3.0-beta.3-Linux/
PATH=$PATH:$HHLIB/bin:$HHLIB/scripts
    
#echo "$hhblits3_dir/bin/hhblits -i $FASTA -d $hhblitdb -oa3m $targetid.a3m -n 3"
$hhblits3_dir/bin/hhblits -i $FASTA -d $hhblitdb -oa3m $targetid.a3m -n 3
$hhblits3_dir/bin/hhmake -i $targetid.a3m -o $targetid.hmm

### or try less strict
if [ -f "$targetid.hmm" ]
then
	echo "$targetid.hmm found."
fi


