# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 21:37:04 2017

@author: Jie Hou, Zhiye
"""
import sys
import os
from shutil import copyfile
import platform

if len(sys.argv) != 11:
  print ('please input the right parameters')
  sys.exit(1)
current_os_name = platform.platform()
print ('%s' % current_os_name)


GLOBAL_PATH=os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

print (GLOBAL_PATH)
sys.path.insert(0, GLOBAL_PATH+'/lib/')
from Data_loading import *
from Model_training import *


inter=int(sys.argv[1]) #15
nb_filters=int(sys.argv[2]) #10
nb_layers=int(sys.argv[3]) #10
opt=sys.argv[4] #nadam
filtsize=sys.argv[5] #6_10
out_epoch=int(sys.argv[6]) #100
in_epoch=int(sys.argv[7]) #3
feature_dir = sys.argv[8]
outputdir = sys.argv[9]
batchsize = int(sys.argv[10])

test_datafile=GLOBAL_PATH+'/datasets/dnss2_val.lst'
train_datafile=GLOBAL_PATH+'/datasets/dnss2_train.lst'
val_datafile=GLOBAL_PATH+'/datasets/dnss2_val.lst'

CV_dir=outputdir+'/filter'+str(nb_filters)+'_layers'+str(nb_layers)+'_inter'+str(inter)+'_opt'+str(opt)+'_ftsize'+str(filtsize)+'_batchsize'+str(batchsize)

lib_dir=GLOBAL_PATH+'/lib/'

import tensorflow as tf
config = tf.ConfigProto(allow_soft_placement = True)
tf.GPUOptions(per_process_gpu_memory_fraction = 0.8)
config.gpu_options.allow_growth = True
sess= tf.Session(config = config)

filetsize_array = list(map(int,filtsize.split("_")))

if not os.path.exists(CV_dir):
    os.makedirs(CV_dir)
# else:
#     print 'File exit, quit!'
#     sys.exit(1)



import time
data_all_dict_padding = load_train_test_data_padding_with_interval(train_datafile, feature_dir, inter,5000, 'deepss_1dInception')
testdata_all_dict_padding = load_train_test_data_padding_with_interval(val_datafile, feature_dir, inter,5000, 'deepss_1dInception')

start_time = time.time()
DeepSS_1dconv_train_win_filter_layer_opt_fast(data_all_dict_padding,testdata_all_dict_padding,train_datafile,test_datafile,val_datafile,CV_dir, feature_dir,"deepss_1dInception",out_epoch,in_epoch,inter,5000,filetsize_array,True,'sigmoid',nb_filters,nb_layers,opt,lib_dir, batchsize)

print("--- %s seconds ---" % (time.time() - start_time))
