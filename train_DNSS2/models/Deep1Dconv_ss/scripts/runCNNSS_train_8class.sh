#!/bin/bash -l
#SBATCH -J  CNNSS
#SBATCH -o CNNSS-%j.out
#SBATCH -p Lewis,hpc4,hpc5
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 2-00:00
#SBATCH --mem 20G

source /storage/hpc/scratch/zggc9/keras_theano/keras_virtual_env/bin/activate
echo "Training secondary structure"
module load R/R-3.3.1

export HDF5_USE_FILE_LOCKING=FALSE
temp_dir=$(pwd)
GLOBAL_PATH=${temp_dir%%DNSS2*}'DNSS2'
feature_dir=$GLOBAL_PATH/datasets/newfeature/features_win1_with_atch_hmmEm_hmmTr_hhblitsMSA_no_aa_8class
output_dir=$GLOBAL_PATH/output/model_train_CNNSS_win1
acclog_dir=$GLOBAL_PATH/output/evaluate/paper_review

# python $GLOBAL_PATH/models/Deep1Dconv_ss/scripts/train_deepcov_ss_8class.py  15  40 5 nadam '6'  100 3  $feature_dir $output_dir 25
python $GLOBAL_PATH/lib/test_dnss_8class.py $GLOBAL_PATH/datasets/dnss2_train.lst  15 40 5 nadam '6' $feature_dir $output_dir $acclog_dir 'deepss_1dconv' 'train' 25
python $GLOBAL_PATH/lib/test_dnss_8class.py $GLOBAL_PATH/datasets/dnss2_val.lst  15 40 5 nadam '6' $feature_dir $output_dir $acclog_dir 'deepss_1dconv' 'evalu' 25
python $GLOBAL_PATH/lib/test_dnss_8class.py $GLOBAL_PATH/datasets/dnss2-blind-test.lst  15 40 5 nadam '6' $feature_dir $output_dir $acclog_dir 'deepss_1dconv' 'test' 25
