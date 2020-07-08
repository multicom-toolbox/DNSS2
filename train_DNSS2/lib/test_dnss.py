# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 21:47:26 2017

@author: Jie Hou
"""
import sys
import os
import platform
from Model_construct import DeepCov_SS_with_paras
from keras.models import model_from_json
import numpy as np
import time
import shutil
import shlex, subprocess
from subprocess import Popen, PIPE

import tensorflow as tf
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

if len(sys.argv) != 13:
          print ('please input the right parameters')
          sys.exit(1)
current_os_name = platform.platform()
print ('%s' % current_os_name)
GLOBAL_PATH=os.path.dirname(os.path.dirname(__file__))

sys.path.insert(0, GLOBAL_PATH+'/lib/')

test_list=(sys.argv[1]) 
inter=int(sys.argv[2]) #15
nb_filters=int(sys.argv[3]) #10
nb_layers=int(sys.argv[4]) #10
opt=sys.argv[5] #nadam
filtsize=sys.argv[6] #6_10
feature_dir=sys.argv[7]
output_dir=(sys.argv[8]) 
acclog_dir=(sys.argv[9])
model_prefix=(sys.argv[10])
tag=(sys.argv[11])
batchsize = int(sys.argv[12])

if (batchsize > 0):
    CV_dir=output_dir+'/filter'+str(nb_filters)+'_layers'+str(nb_layers)+'_inter'+str(inter)+'_opt'+str(opt)+'_ftsize'+str(filtsize)+'_batchsize'+str(batchsize)
else:
    CV_dir=output_dir+'/filter'+str(nb_filters)+'_layers'+str(nb_layers)+'_inter'+str(inter)+'_opt'+str(opt)+'_ftsize'+str(filtsize)



model_in=CV_dir + '/model-train-%s.json' % model_prefix
model_weight_in= CV_dir + '/model-train-weight-%s-best-val.h5' % model_prefix

if not os.path.isfile(model_in):
    print ("model file not exists: ",model_in, " error!")
    exit(1)

if not os.path.isfile(model_weight_in):
    print ("weight file not exists: ",model_weight_in, " error!")
    exit(1)


Testlist_data_keys = dict()
Testlist_targets_keys = dict()
sequence_file=open(test_list,'r').readlines() 
for i in xrange(len(sequence_file)):
    pdb_name = sequence_file[i].rstrip()
    featurefile = feature_dir + '/' + pdb_name + '.fea'
    if not os.path.isfile(featurefile):
                print ("feature file not exists: ",featurefile, " pass!")
                continue           
    
    featuredata = np.loadtxt(featurefile) #(169, 51)
    fea_len = featuredata.shape[0]
    test_labels = featuredata[:,0:3]#(169, 3)
    test_feature = featuredata[:,3:] #(169, 48)    
    if pdb_name in Testlist_data_keys:
        print ("Duplicate pdb name %s in Test list " % pdb_name)
    else:
        Testlist_data_keys[pdb_name]=test_feature.reshape(1,test_feature.shape[0],test_feature.shape[1])
    
    if pdb_name in Testlist_targets_keys:
        print ("Duplicate pdb name %s in Test list " % pdb_name)
    else:
        Testlist_targets_keys[pdb_name]=test_labels.reshape(1,test_labels.shape[0],test_labels.shape[1])

### Define the model 
if os.path.exists(model_in):
    print ("######## Loading existing model ",model_in)
    # load json and create model
    json_file_model = open(model_in, 'r')
    loaded_model_json = json_file_model.read()
    json_file_model.close()
    
    print("######## Loaded model from disk")
    DeepSS_CNN = model_from_json(loaded_model_json)       
else:
    print ("######## Failed to find model",model_in)
    sys.exit(1)

if os.path.exists(model_weight_in):
    print ("######## Loading existing weights ",model_weight_in)
    DeepSS_CNN.load_weights(model_weight_in)
    DeepSS_CNN.compile(loss="categorical_crossentropy", metrics=['accuracy'], optimizer='nadam')
else:
    print ("######## Failed to find weight",model_weight_in)
    sys.exit(1)
 

## start evaluate the dataset
sequence_file=open(test_list,'r').readlines() 
predir = CV_dir + '/test_prediction/'
chkdirs(predir)
dnssdir = CV_dir + '/test_prediction_dnss/'
chkdirs(dnssdir)
eva_dir = CV_dir + '/test_prediction_q3_sov_log_loss/'
chkdirs(eva_dir)

test_acc=0.0;
acc_num=0;
test_loss=0.0;
loss_num=0;
for i in xrange(len(sequence_file)):
    pdb_name = sequence_file[i].rstrip()
    test_featuredata_all=Testlist_data_keys[pdb_name]
    test_targets=Testlist_targets_keys[pdb_name]
    score, accuracy = DeepSS_CNN.evaluate([test_featuredata_all], test_targets, batch_size=10, verbose=0)
    
    test_acc += accuracy
    acc_num += 1
    
    test_loss += score
    loss_num += 1
    predict_val= DeepSS_CNN.predict([test_featuredata_all])
    targsize=3
    predict_val= predict_val.reshape(predict_val.shape[1],predict_val.shape[2])
    max_vals = np.reshape(np.repeat(predict_val.max(axis=1), targsize), (predict_val.shape[0], targsize))
    preds = 1 * (predict_val > max_vals - .0001)
    predfile = predir + pdb_name + ".pred";
    probfile = predir + pdb_name + ".prob";
    np.savetxt(predfile, preds, fmt="%d")
    np.savetxt(probfile, predict_val, fmt="%.6f")
    
    del test_featuredata_all
    del test_targets

test_acc /= acc_num 
test_loss /= loss_num 



############################  summarize the Q3/SOV information
# tag = 'abc'
print("###############\n")
args_str ="perl " + GLOBAL_PATH + "/lib/evaluation_dnss_prediction.pl -pred "  + predir +  " -out " + dnssdir + " -list " + test_list + " -tag " + tag
# #print "Running "+ args_str
args = shlex.split(args_str)
pipe = subprocess.Popen(args, stdin=subprocess.PIPE)

scorefile=dnssdir+'/'+tag+'.score'


time.sleep(10) 
found = 0
while (found == 0):
    print ("Checking file ",scorefile)
    time.sleep(10) 
    if os.path.exists(scorefile):
        found = 1

print ("Score saved to file ",scorefile)

Q3_acc=0
SOV_acc=0
score_content=open(scorefile,'r').readlines() 
for i in xrange(len(score_content)):
    score_line = score_content[i].rstrip()
    if '#Average Q3 score'  in score_line: 
        score_tmp = score_line.split(':')
        Q3_acc = score_tmp[1]
        Q3_acc = Q3_acc.replace(" ", "")
    if '#Average Sov score'  in score_line: 
        score_tmp = score_line.split(':')
        SOV_acc = score_tmp[1]
        SOV_acc = SOV_acc.replace(" ", "")

# CV_dir=output_dir+'/filter'+str(nb_filters)+'_layers'+str(nb_layers)+'_inter'+str(inter)+'_opt'+str(opt)+'_ftsize'+str(filtsize);
time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
acclog_dir_temp = acclog_dir.split('/')
test_list_temp = test_list.split('/')
acc_history_out = "%s/%s.acc_history" % (acclog_dir, acclog_dir_temp[len(acclog_dir_temp)-1])
print ('%s' % acc_history_out)
chkdirs(acc_history_out)
if chkfiles(acc_history_out):
    print ('####acc_file_exist,pass!')
    pass
else:
    print ('####create_acc_file!')
    with open(acc_history_out, "w") as myfile:
        myfile.write("time\t netname\t dataset\t filtnum\t layernum\t kernelsize\t Accuracy\t Loss\t Q3\t SOV\t batchsize\n")
if (batchsize > 0):
    acc_history_content = "%s\t %s\t %s\t %i\t %i\t %s\t %.4f\t %.4f\t %s\t %s\t %i\n" % (time_str, model_prefix, test_list_temp[len(test_list_temp)-1],nb_filters,nb_layers,filtsize,test_acc,test_loss,Q3_acc,SOV_acc, batchsize)
else:
    acc_history_content = "%s\t %s\t %s\t %i\t %i\t %s\t %.4f\t %.4f\t %s\t %s\n" % (time_str, model_prefix, test_list_temp[len(test_list_temp)-1],nb_filters,nb_layers,filtsize,test_acc,test_loss,Q3_acc,SOV_acc)
with open(acc_history_out, "a") as myfile: myfile.write(acc_history_content)  

print ("Test Accuracy %.5f" %(test_acc))
print ("Test Loss %.5f" %(test_loss))
print ("Q3 Accuracy %sf" %(Q3_acc))
print ("SOV Accuracy %sf" %(SOV_acc))




