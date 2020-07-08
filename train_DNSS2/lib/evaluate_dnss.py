# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 21:47:26 2017

@author: Jie Hou
"""

import sys
import os
from Custom_class import remove_1d_padding
from Model_construct import DeepCov_SS_with_paras
from keras.models import model_from_json
import numpy as np
import time
import shutil
import shlex, subprocess
from subprocess import Popen, PIPE

def chkdirs(fn):
  dn = os.path.dirname(fn)
  if not os.path.exists(dn): os.makedirs(dn)

if len(sys.argv) != 4:
          print 'please input the right parameters'
          sys.exit(1)
GLOBAL_PATH='/mnt/data/zhiye/Python/DNSS2'
sys.path.insert(0, GLOBAL_PATH+'/lib/')

test_list=(sys.argv[1]) #15
tag=sys.argv[2]
CV_dir=(sys.argv[3]) #10

## start evaluate the dataset
sequence_file=open(test_list,'r').readlines() 
predir = CV_dir + '/test_prediction/'
chkdirs(predir)
dnssdir = CV_dir + '/test_prediction_dnss/'
chkdirs(dnssdir)
eva_dir = CV_dir + '/test_prediction_q3_sov_log_loss/'
chkdirs(eva_dir)

args_str ="perl " + GLOBAL_PATH + "/lib/evaluation_dnss_prediction.pl -pred "  + predir +  " -out " + dnssdir + " -list " + test_list + " -tag " + tag
#print "Running "+ args_str
args = shlex.split(args_str)
pipe = subprocess.Popen(args, stdin=subprocess.PIPE)

scorefile=dnssdir+'/'+tag+'.score'

found = 0
while (found == 0):
    print "Checking file ",scorefile
    time.sleep(10) 
    if os.path.exists(scorefile):
      found = 1

shutil.copy2(scorefile, eva_dir)
print "Score saved to file ",eva_dir
## clean for next iteration
#shutil.rmtree(predir)
#shutil.rmtree(dnssdir)
    

