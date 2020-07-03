# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 21:47:26 2017

@author: Jie Hou
"""
import sys
import os
import platform
from keras.models import model_from_json
import numpy as np
import time
import shutil
import shlex, subprocess
from subprocess import Popen, PIPE

import tensorflow as tf

if tf.test.gpu_device_name():
    print('GPU found')
else:
    print("No GPU found")

config = tf.ConfigProto(allow_soft_placement = True)
tf.GPUOptions(per_process_gpu_memory_fraction = 0.4)
config.gpu_options.allow_growth = True
sess= tf.Session(config = config)

def chkdirs(fn):
  dn = os.path.dirname(fn)
  if not os.path.exists(dn): os.makedirs(dn)

def chkfiles(fn):
  if os.path.exists(fn):
    return True 
  else:
    return False

if len(sys.argv) != 6:
          print 'please input the right parameters'
          sys.exit(1)

model_in=(sys.argv[1]) 
model_weight_in=sys.argv[2] #15
featurefile=sys.argv[3] #10
predfile=sys.argv[4] #10
probfile=sys.argv[5] #nadam

#model_in=CV_dir + '/model-train-%s.json' % model_prefix
#model_weight_in= CV_dir + '/model-train-weight-%s-best-val.h5' % model_prefix

if not os.path.isfile(model_in):
    print "model file not exists: ",model_in, " error!"
    exit(1)

if not os.path.isfile(model_weight_in):
    print "weight file not exists: ",model_weight_in, " error!"
    exit(1)

### Define the model 
if os.path.exists(model_in):
    print "######## Loading existing model ",model_in;
    # load json and create model
    json_file_model = open(model_in, 'r')
    loaded_model_json = json_file_model.read()
    json_file_model.close()
    
    print("######## Loaded model from disk")
    #DeepSS_CNN = model_from_json(loaded_model_json, custom_objects={'remove_1d_padding': remove_1d_padding}) 
    DeepSS_CNN = model_from_json(loaded_model_json)       
else:
    print "######## Failed to find model",model_in;
    sys.exit(1)

if os.path.exists(model_weight_in):
    print "######## Loading existing weights ",model_weight_in;
    DeepSS_CNN.load_weights(model_weight_in)
    DeepSS_CNN.compile(loss="categorical_crossentropy", metrics=['accuracy'], optimizer='nadam')
else:
    print "######## Failed to find weight",model_weight_in;
    sys.exit(1)
 

#featurefile = feature_dir + '/' + pdb_name + '.fea'
featuredata = np.loadtxt(featurefile) #(169, 51)
fea_len = featuredata.shape[0]
test_featuredata_all = featuredata.reshape(1,featuredata.shape[0],featuredata.shape[1])

predict_val= DeepSS_CNN.predict([test_featuredata_all])
targsize=8
predict_val= predict_val.reshape(predict_val.shape[1],predict_val.shape[2])
max_vals = np.reshape(np.repeat(predict_val.max(axis=1), targsize), (predict_val.shape[0], targsize))
#print "".format(predict_val[0], max_vals[0], (predict_val[0] >= max_vals[0]))
preds = 1 * (predict_val > max_vals - .0001)
#predfile = predir + pdb_name + ".pred";
#probfile = predir + pdb_name + ".prob";
np.savetxt(predfile, preds, fmt="%d")
np.savetxt(probfile, predict_val, fmt="%.6f")
    
 

