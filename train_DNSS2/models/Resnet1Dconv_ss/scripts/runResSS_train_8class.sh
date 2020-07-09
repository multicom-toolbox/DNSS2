#!/bin/bash -l
#SBATCH -J  RESSS
#SBATCH -o RESSS-%j.out
#SBATCH -p Lewis,hpc4,hpc5
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 2-00:00
#SBATCH --mem 20G


export HDF5_USE_FILE_LOCKING=FALSE
temp_dir=$(pwd)
GLOBAL_PATH=${temp_dir%%DNSS2*}'DNSS2'

if [ -f "$GLOBAL_PATH/python_virtualenv_DNSS2/bin/activate" ]; then
    echo "Found virutal environment $GLOBAL_PATH/python_virtualenv_DNSS2/bin/activate."
else
    echo "Virtual environment ($GLOBAL_PATH/python_virtualenv_DNSS2/bin/activate) does not exist."
    exit
fi

source $GLOBAL_PATH/python_virtualenv_DNSS2/bin/activate
echo "Training 8-class secondary structure"

feature_dir=$GLOBAL_PATH/train_DNSS2/datasets/features_8class
output_dir=$GLOBAL_PATH/output/model_train_ResSS_8class
acclog_dir=$GLOBAL_PATH/output//evaluate_8class/

python $GLOBAL_PATH/train_DNSS2/models/Resnet1Dconv_ss/scripts/train_deepcovResnet_ss_8class.py  15 37 41 nadam '5'  100 3  $feature_dir $output_dir 25
#python $GLOBAL_PATH/lib/test_dnss.py $GLOBAL_PATH/datasets/dnss2_train.lst  15 37 41 nadam '5' $feature_dir $output_dir $acclog_dir 'deepss_1dResnet' 'train' 25
#python $GLOBAL_PATH/lib/test_dnss.py $GLOBAL_PATH/datasets/dnss2_val.lst  15 37 41 nadam '5' $feature_dir $output_dir $acclog_dir 'deepss_1dResnet' 'evalu' 25
#python $GLOBAL_PATH/lib/test_dnss.py $GLOBAL_PATH/datasets/dnss2-blind-test.lst  15 37 41 nadam '5' $feature_dir $output_dir $acclog_dir 'deepss_1dResnet' 'test' 25
