# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 21:47:26 2017

@author: Jie Hou, Zhiye
"""
import os
from Model_construct import *
from keras.models import model_from_json,load_model, Sequential
import numpy as np
import time
import shutil
import shlex, subprocess
from subprocess import Popen, PIPE

from collections import defaultdict
import pickle
from PIL import Image

import keras.backend as K
from keras.datasets import mnist
from keras.layers import Input, Dense, Reshape, Flatten, Embedding, merge, Dropout
from keras.models import Sequential, Model
from keras.optimizers import Adam
from keras.utils.generic_utils import Progbar
from keras.constraints import maxnorm

from keras.models import Model
from keras.layers import Activation, Dense, Dropout, Flatten, Input, Merge, Convolution1D, Convolution2D
from keras.layers.normalization import BatchNormalization



def chkdirs(fn):
  dn = os.path.dirname(fn)
  if not os.path.exists(dn): os.makedirs(dn)    

def DeepSS_1dconv_train_win_filter_layer_opt_fast(data_all_dict_padding,testdata_all_dict_padding,train_list,test_list,val_list,CV_dir,feature_dir,model_prefix,epoch_outside,epoch_inside,interval_len,seq_end,win_array,use_bias,hidden_type,nb_filters,nb_layers,opt,lib_dir, batch_size_train): #/storage/htc/bdm/jh7x8/GANSS/Deep1Dconv_ss/lib/
    start=0
    end=seq_end
    import numpy as np
    Train_data_keys = dict()
    Train_targets_keys = dict()
    Test_data_keys = dict()
    Test_targets_keys = dict()
    
    feature_num=0; # the number of features for each residue
    for key in data_all_dict_padding.keys():
        if key <start: # run first model on 100 at most
            continue
        if key > end: # run first model on 100 at most
            continue
        print ('### Loading sequence length :', key)
        seq_len=key
        trainfeaturedata = data_all_dict_padding[key]
        train_labels = trainfeaturedata[:,:,0:3]
        train_feature = trainfeaturedata[:,:,3:]
        feature_num=train_feature.shape[2]
        if seq_len in testdata_all_dict_padding:
            testfeaturedata = testdata_all_dict_padding[seq_len]
            #print "Loading test dataset "
        else:
            testfeaturedata = trainfeaturedata
            print ("\n\n##Warning: Setting training dataset as testing dataset \n\n")
        
        
        test_labels = testfeaturedata[:,:,0:3]
        test_feature = testfeaturedata[:,:,3:]    
        sequence_length = seq_len
        
        if seq_len in Train_data_keys:
            raise Exception("Duplicate seq length %i in Train list, since it has been combined when loading data " % seq_len)
        else:
            Train_data_keys[seq_len]=(train_feature)
            
        if seq_len in Train_targets_keys:
            raise Exception("Duplicate seq length %i in Train list, since it has been combined when loading data " % seq_len)
        else:
            Train_targets_keys[seq_len]=train_labels        
        #processing test data 
        if seq_len in Test_data_keys:
            raise Exception("Duplicate seq length %i in Test list, since it has been combined when loading data " % seq_len)
        else:
            Test_data_keys[seq_len]=test_feature 
        
        if seq_len in Test_targets_keys:
            raise Exception("Duplicate seq length %i in Test list, since it has been combined when loading data " % seq_len)
        else:
            Test_targets_keys[seq_len]=test_labels
    
    Trainlist_data_keys = dict()
    Trainlist_targets_keys = dict()
    sequence_file=open(train_list,'r').readlines() 
    for i in range(len(sequence_file)):
        pdb_name = sequence_file[i].rstrip()
        featurefile = feature_dir + '/' + pdb_name + '.fea'
        if not os.path.isfile(featurefile):
                    print ("feature file not exists: ",featurefile, " pass!")
                    continue           
        
        featuredata = np.loadtxt(featurefile) #(169, 51)
        fea_len = featuredata.shape[0]
        train_labels = featuredata[:,0:3]#(169, 3)
        train_feature = featuredata[:,3:] #(169, 48)        
        if pdb_name in Trainlist_data_keys:
            print ("Duplicate pdb name %s in Train list " % pdb_name)
        else:
            Trainlist_data_keys[pdb_name]=train_feature.reshape(1,train_feature.shape[0],train_feature.shape[1])
        
        if pdb_name in Trainlist_targets_keys:
            print ("Duplicate pdb name %s in Train list " % pdb_name)
        else:
            Trainlist_targets_keys[pdb_name]=train_labels.reshape(1,train_labels.shape[0],train_labels.shape[1])
    
    Testlist_data_keys = dict()
    Testlist_targets_keys = dict()
    sequence_file=open(test_list,'r').readlines() 
    for i in range(len(sequence_file)):
        pdb_name = sequence_file[i].rstrip()
        #print "Loading ",pdb_name
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
    
    Vallist_data_keys = dict()
    Vallist_targets_keys = dict()
    sequence_file=open(val_list,'r').readlines() 
    for i in range(len(sequence_file)):
        pdb_name = sequence_file[i].rstrip()
        #print "Loading ",pdb_name
        featurefile = feature_dir + '/' + pdb_name + '.fea'
        if not os.path.isfile(featurefile):
                    print ("feature file not exists: ",featurefile, " pass!")
                    continue           
        
        featuredata = np.loadtxt(featurefile) #(169, 51)
        fea_len = featuredata.shape[0]
        val_labels = featuredata[:,0:3]#(169, 3)
        val_feature = featuredata[:,3:] #(169, 48)  
        if pdb_name in Vallist_data_keys:
            print ("Duplicate pdb name %s in Val list " % pdb_name)
        else:
            Vallist_data_keys[pdb_name]=val_feature.reshape(1,val_feature.shape[0],val_feature.shape[1])
        
        if pdb_name in Vallist_targets_keys:
            print ("Duplicate pdb name %s in Val list " % pdb_name)
        else:
            Vallist_targets_keys[pdb_name]=val_labels.reshape(1,val_labels.shape[0],val_labels.shape[1])
    
    ### Define the model 
    model_out= "%s/model-train-%s.json" % (CV_dir,model_prefix)
    model_weight_out = "%s/model-train-weight-%s.h5" % (CV_dir,model_prefix)
    model_weight_out_best = "%s/model-train-weight-%s-best-val.h5" % (CV_dir,model_prefix)
    
    if os.path.exists(model_out):
        print ("######## Loading existing model ",model_out)
        # load json and create model
        json_file_model = open(model_out, 'r')
        loaded_model_json = json_file_model.read()
        json_file_model.close()
        
        print("######## Loaded model from disk")
        #DeepSS_CNN = model_from_json(loaded_model_json, custom_objects={'remove_1d_padding': remove_1d_padding}) 
        DeepSS_CNN = model_from_json(loaded_model_json)       
    else:
        print ("######## Setting initial model")
        ## ktop_node is the length of input proteins
        if model_prefix == 'deepss_1dconv':
            # opt = Adam(lr=0.002, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.0)
            DeepSS_CNN = DeepCov_SS_with_paras(win_array,feature_num,use_bias,hidden_type,nb_filters,nb_layers,opt)
        elif model_prefix == 'deepss_1dRCNN':
            DeepSS_CNN = DeepCovRCNN_with_paras(win_array,feature_num,use_bias,hidden_type,nb_filters,nb_layers,opt)
        elif model_prefix == 'deepss_1dResnet':
            DeepSS_CNN = DeepResnet1D_with_paras(win_array,feature_num,use_bias,hidden_type,nb_filters,nb_layers,opt)
        elif model_prefix == 'deepss_1dInception':
            DeepSS_CNN = DeepInception1D_with_paras(win_array,feature_num,use_bias,hidden_type,nb_filters,nb_layers,opt)
        elif model_prefix == 'deepss_1dCRMN':
            DeepSS_CNN = DeepCRMN_SS_with_paras(win_array,feature_num,use_bias,hidden_type,nb_filters,nb_layers,opt)
        elif model_prefix == 'deepss_1dFrac':
            DeepSS_CNN = DeepFracNet_SS_with_paras(win_array,feature_num,use_bias,hidden_type,nb_filters,nb_layers,opt)
        else:
            DeepSS_CNN = DeepCov_SS_with_paras(win_array,feature_num,use_bias,hidden_type,nb_filters,nb_layers,opt)

    if os.path.exists(model_weight_out):
        print ("######## Loading existing weights ",model_weight_out)
        DeepSS_CNN.load_weights(model_weight_out)
        DeepSS_CNN.compile(loss="categorical_crossentropy", metrics=['accuracy'], optimizer=opt)
    else:
        print ("######## Setting initial weights")
        DeepSS_CNN.compile(loss="categorical_crossentropy", metrics=['accuracy'], optimizer=opt)
     
 
    train_acc_best = 0 
    test_acc_best = 0 
    val_acc_best = 0

    test_acc_history_out = "%s/testing.acc_history" % (CV_dir)
    chkdirs(test_acc_history_out)     
    with open(test_acc_history_out, "w") as myfile:
      myfile.write("Interval_len\tEpoch_outside\tEpoch_inside\tAccuracy\tLoss\n")
      
    train_acc_history_out = "%s/training.acc_history" % (CV_dir)
    chkdirs(train_acc_history_out)     
    with open(train_acc_history_out, "w") as myfile:
      myfile.write("Interval_len\tEpoch_outside\tEpoch_inside\tAccuracy\tLoss\n")
      
    val_acc_history_out = "%s/validation.acc_history" % (CV_dir)
    chkdirs(val_acc_history_out)     
    with open(val_acc_history_out, "w") as myfile:
      myfile.write("Interval_len\tEpoch_outside\tEpoch_inside\tAccuracy\tLoss\n")
    

    val_sequence_file=open(val_list,'r').readlines() 
    for epoch in range(0,epoch_outside):
        print ("\n############ Running epoch ", epoch )
        
        for key in data_all_dict_padding.keys():
            if key <start: # run first model on 100 at most
                continue
            if key > end: # run first model on 100 at most
                continue
            print ('### Loading sequence length :', key)
            seq_len=key
            
            train_featuredata_all=Train_data_keys[seq_len]
            train_targets=Train_targets_keys[seq_len]
            test_featuredata_all=Test_data_keys[seq_len]
            test_targets=Test_targets_keys[seq_len]
            print ("Train shape: ",train_featuredata_all.shape, " in outside epoch ", epoch)
            print ("Test shape: ",test_featuredata_all.shape, " in outside epoch ", epoch)

            DeepSS_CNN.fit([train_featuredata_all], train_targets, batch_size= batch_size_train, nb_epoch=epoch_inside,  validation_data=([test_featuredata_all], test_targets), verbose=1)

            del train_featuredata_all
            del train_targets
            del test_featuredata_all
            del test_targets
            
        
        ##### running validation
        print ("Now evaluate")
        
        val_acc=0.0;
        acc_num=0;
        val_loss=0.0;
        loss_num=0;
        for i in range(len(val_sequence_file)):
            pdb_name = val_sequence_file[i].rstrip()
            val_featuredata_all=Vallist_data_keys[pdb_name]
            val_targets=Vallist_targets_keys[pdb_name]
            score, accuracy = DeepSS_CNN.evaluate([val_featuredata_all], val_targets, batch_size=10, verbose=0)
            val_acc += accuracy
            acc_num += 1
            
            val_loss += score
            loss_num += 1                   
            del val_featuredata_all
            del val_targets
        
        val_acc /= acc_num 
        val_loss /= loss_num 
        
        val_acc_history_content = "%i\t%i\t%i\t%.4f\t%.4f\n" % (interval_len,epoch,epoch_inside,val_acc,val_loss)
        with open(val_acc_history_out, "a") as myfile:
                    myfile.write(val_acc_history_content)  
        
        if val_acc >= val_acc_best:
            val_acc_best = val_acc 
            score_imed = "Accuracy of Val: %.4f\t\n" % (val_acc_best)
            model_json = DeepSS_CNN.to_json()
            print("Saved model to disk")
            with open(model_out, "w") as json_file:
                json_file.write(model_json)
            print ("Saved best weight to disk, ", score_imed)
            DeepSS_CNN.save_weights(model_weight_out_best)
        print ('The val accuracy is %.5f' % (val_acc))

    print ("Training finished, best validation acc = ",val_acc_best)
    print ("Setting and saving best weights")
    DeepSS_CNN.load_weights(model_weight_out_best)
    DeepSS_CNN.save_weights(model_weight_out)

def DeepSS_1dconv_train_win_filter_layer_opt_fast_8class(data_all_dict_padding,testdata_all_dict_padding,train_list,test_list,val_list,CV_dir,feature_dir,model_prefix,epoch_outside,epoch_inside,interval_len,seq_end,win_array,use_bias,hidden_type,nb_filters,nb_layers,opt,lib_dir, batch_size_train): #/storage/htc/bdm/jh7x3/GANSS/Deep1Dconv_ss/lib/
    start=0
    end=seq_end
    import numpy as np
    Train_data_keys = dict()
    Train_targets_keys = dict()
    Test_data_keys = dict()
    Test_targets_keys = dict()
    
    feature_num=0; # the number of features for each residue
    for key in data_all_dict_padding.keys():
        if key <start: # run first model on 100 at most
            continue
        if key > end: # run first model on 100 at most
            continue
        print ('### Loading sequence length :', key)
        seq_len=key
        trainfeaturedata = data_all_dict_padding[key]
        train_labels = trainfeaturedata[:,:,0:8]
        train_feature = trainfeaturedata[:,:,8:]
        feature_num=train_feature.shape[2]
        if seq_len in testdata_all_dict_padding:
            testfeaturedata = testdata_all_dict_padding[seq_len]
            #print "Loading test dataset "
        else:
            testfeaturedata = trainfeaturedata
            print ("\n\n##Warning: Setting training dataset as testing dataset \n\n")
        
        
        test_labels = testfeaturedata[:,:,0:8]
        test_feature = testfeaturedata[:,:,8:]    
        sequence_length = seq_len
        
        if seq_len in Train_data_keys:
            raise Exception("Duplicate seq length %i in Train list, since it has been combined when loading data " % seq_len)
        else:
            Train_data_keys[seq_len]=(train_feature)
            
        if seq_len in Train_targets_keys:
            raise Exception("Duplicate seq length %i in Train list, since it has been combined when loading data " % seq_len)
        else:
            Train_targets_keys[seq_len]=train_labels        
        #processing test data 
        if seq_len in Test_data_keys:
            raise Exception("Duplicate seq length %i in Test list, since it has been combined when loading data " % seq_len)
        else:
            Test_data_keys[seq_len]=test_feature 
        
        if seq_len in Test_targets_keys:
            raise Exception("Duplicate seq length %i in Test list, since it has been combined when loading data " % seq_len)
        else:
            Test_targets_keys[seq_len]=test_labels
    
    Trainlist_data_keys = dict()
    Trainlist_targets_keys = dict()
    sequence_file=open(train_list,'r').readlines() 
    for i in range(len(sequence_file)):
        pdb_name = sequence_file[i].rstrip()
        featurefile = feature_dir + '/' + pdb_name + '.fea'
        if not os.path.isfile(featurefile):
                    print ("feature file not exists: ",featurefile, " pass!")
                    continue           
        
        featuredata = np.loadtxt(featurefile) #(169, 51)
        fea_len = featuredata.shape[0]
        train_labels = featuredata[:,0:8]#(169, 8)
        train_feature = featuredata[:,8:] #(169, 48)        
        if pdb_name in Trainlist_data_keys:
            print ("Duplicate pdb name %s in Train list " % pdb_name)
        else:
            Trainlist_data_keys[pdb_name]=train_feature.reshape(1,train_feature.shape[0],train_feature.shape[1])
        
        if pdb_name in Trainlist_targets_keys:
            print ("Duplicate pdb name %s in Train list " % pdb_name)
        else:
            Trainlist_targets_keys[pdb_name]=train_labels.reshape(1,train_labels.shape[0],train_labels.shape[1])
    
    Testlist_data_keys = dict()
    Testlist_targets_keys = dict()
    sequence_file=open(test_list,'r').readlines() 
    for i in range(len(sequence_file)):
        pdb_name = sequence_file[i].rstrip()
        #print "Loading ",pdb_name
        featurefile = feature_dir + '/' + pdb_name + '.fea'
        if not os.path.isfile(featurefile):
                    print ("feature file not exists: ",featurefile, " pass!")
                    continue           
        
        featuredata = np.loadtxt(featurefile) #(169, 51)
        fea_len = featuredata.shape[0]
        test_labels = featuredata[:,0:8]#(169, 8)
        test_feature = featuredata[:,8:] #(169, 48)    
        if pdb_name in Testlist_data_keys:
            print ("Duplicate pdb name %s in Test list " % pdb_name)
        else:
            Testlist_data_keys[pdb_name]=test_feature.reshape(1,test_feature.shape[0],test_feature.shape[1])
        
        if pdb_name in Testlist_targets_keys:
            print ("Duplicate pdb name %s in Test list " % pdb_name)
        else:
            Testlist_targets_keys[pdb_name]=test_labels.reshape(1,test_labels.shape[0],test_labels.shape[1])
    
    Vallist_data_keys = dict()
    Vallist_targets_keys = dict()
    sequence_file=open(val_list,'r').readlines() 
    for i in range(len(sequence_file)):
        pdb_name = sequence_file[i].rstrip()
        #print "Loading ",pdb_name
        featurefile = feature_dir + '/' + pdb_name + '.fea'
        if not os.path.isfile(featurefile):
                    print ("feature file not exists: ",featurefile, " pass!")
                    continue           
        
        featuredata = np.loadtxt(featurefile) #(169, 51)
        fea_len = featuredata.shape[0]
        val_labels = featuredata[:,0:8]#(169, 8)
        val_feature = featuredata[:,8:] #(169, 48)  
        if pdb_name in Vallist_data_keys:
            print ("Duplicate pdb name %s in Val list " % pdb_name)
        else:
            Vallist_data_keys[pdb_name]=val_feature.reshape(1,val_feature.shape[0],val_feature.shape[1])
        
        if pdb_name in Vallist_targets_keys:
            print ("Duplicate pdb name %s in Val list " % pdb_name)
        else:
            Vallist_targets_keys[pdb_name]=val_labels.reshape(1,val_labels.shape[0],val_labels.shape[1])
    
    ### Define the model 
    model_out= "%s/model-train-%s.json" % (CV_dir,model_prefix)
    model_weight_out = "%s/model-train-weight-%s.h5" % (CV_dir,model_prefix)
    model_weight_out_best = "%s/model-train-weight-%s-best-val.h5" % (CV_dir,model_prefix)
    
    if os.path.exists(model_out):
        print ("######## Loading existing model ",model_out)
        # load json and create model
        json_file_model = open(model_out, 'r')
        loaded_model_json = json_file_model.read()
        json_file_model.close()
        
        print("######## Loaded model from disk")
        #DeepSS_CNN = model_from_json(loaded_model_json, custom_objects={'remove_1d_padding': remove_1d_padding}) 
        DeepSS_CNN = model_from_json(loaded_model_json)       
    else:
        print ("######## Setting initial model")
        ## ktop_node is the length of input proteins
        if model_prefix == 'deepss_1dconv':
            # opt = Adam(lr=0.002, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.0)
            DeepSS_CNN = DeepCov_SS_with_paras(win_array,feature_num,use_bias,hidden_type,nb_filters,nb_layers,opt)
        elif model_prefix == 'deepss_1dRCNN':
            DeepSS_CNN = DeepCovRCNN_with_paras(win_array,feature_num,use_bias,hidden_type,nb_filters,nb_layers,opt)
        elif model_prefix == 'deepss_1dResnet':
            DeepSS_CNN = DeepResnet1D_with_paras(win_array,feature_num,use_bias,hidden_type,nb_filters,nb_layers,opt)
        elif model_prefix == 'deepss_1dInception':
            DeepSS_CNN = DeepInception1D_with_paras(win_array,feature_num,use_bias,hidden_type,nb_filters,nb_layers,opt)
        elif model_prefix == 'deepss_1dCRMN':
            DeepSS_CNN = DeepCRMN_SS_with_paras(win_array,feature_num,use_bias,hidden_type,nb_filters,nb_layers,opt)
        elif model_prefix == 'deepss_1dFrac':
            DeepSS_CNN = DeepFracNet_SS_with_paras(win_array,feature_num,use_bias,hidden_type,nb_filters,nb_layers,opt)
        else:
            DeepSS_CNN = DeepCov_SS_with_paras(win_array,feature_num,use_bias,hidden_type,nb_filters,nb_layers,opt)

    if os.path.exists(model_weight_out):
        print ("######## Loading existing weights ",model_weight_out)
        DeepSS_CNN.load_weights(model_weight_out)
        DeepSS_CNN.compile(loss="categorical_crossentropy", metrics=['accuracy'], optimizer=opt)
    else:
        print ("######## Setting initial weights")
        DeepSS_CNN.compile(loss="categorical_crossentropy", metrics=['accuracy'], optimizer=opt)
     
 
    train_acc_best = 0 
    test_acc_best = 0 
    val_acc_best = 0

    test_acc_history_out = "%s/testing.acc_history" % (CV_dir)
    chkdirs(test_acc_history_out)     
    with open(test_acc_history_out, "w") as myfile:
      myfile.write("Interval_len\tEpoch_outside\tEpoch_inside\tAccuracy\tLoss\n")
      
    train_acc_history_out = "%s/training.acc_history" % (CV_dir)
    chkdirs(train_acc_history_out)     
    with open(train_acc_history_out, "w") as myfile:
      myfile.write("Interval_len\tEpoch_outside\tEpoch_inside\tAccuracy\tLoss\n")
      
    val_acc_history_out = "%s/validation.acc_history" % (CV_dir)
    chkdirs(val_acc_history_out)     
    with open(val_acc_history_out, "w") as myfile:
      myfile.write("Interval_len\tEpoch_outside\tEpoch_inside\tAccuracy\tLoss\n")
    

    val_sequence_file=open(val_list,'r').readlines() 
    for epoch in range(0,epoch_outside):
        print ("\n############ Running epoch ", epoch )
        
        for key in data_all_dict_padding.keys():
            if key <start: # run first model on 100 at most
                continue
            if key > end: # run first model on 100 at most
                continue
            print ('### Loading sequence length :', key)
            seq_len=key
            
            train_featuredata_all=Train_data_keys[seq_len]
            train_targets=Train_targets_keys[seq_len]
            test_featuredata_all=Test_data_keys[seq_len]
            test_targets=Test_targets_keys[seq_len]
            print ("Train shape: ",train_featuredata_all.shape, " in outside epoch ", epoch)
            print ("Test shape: ",test_featuredata_all.shape, " in outside epoch ", epoch)

            DeepSS_CNN.fit([train_featuredata_all], train_targets, batch_size= batch_size_train, nb_epoch=epoch_inside,  validation_data=([test_featuredata_all], test_targets), verbose=1)

            del train_featuredata_all
            del train_targets
            del test_featuredata_all
            del test_targets
            
        
        ##### running validation
        print ("Now evaluate")
        
        val_acc=0.0;
        acc_num=0;
        val_loss=0.0;
        loss_num=0;
        for i in range(len(val_sequence_file)):
            pdb_name = val_sequence_file[i].rstrip()
            val_featuredata_all=Vallist_data_keys[pdb_name]
            val_targets=Vallist_targets_keys[pdb_name]
            score, accuracy = DeepSS_CNN.evaluate([val_featuredata_all], val_targets, batch_size=10, verbose=0)
            val_acc += accuracy
            acc_num += 1
            
            val_loss += score
            loss_num += 1                   
            del val_featuredata_all
            del val_targets
        
        val_acc /= acc_num 
        val_loss /= loss_num 
        
        val_acc_history_content = "%i\t%i\t%i\t%.4f\t%.4f\n" % (interval_len,epoch,epoch_inside,val_acc,val_loss)
        with open(val_acc_history_out, "a") as myfile:
                    myfile.write(val_acc_history_content)  
        
        if val_acc >= val_acc_best:
            val_acc_best = val_acc 
            score_imed = "Accuracy of Val: %.4f\t\n" % (val_acc_best)
            model_json = DeepSS_CNN.to_json()
            print("Saved model to disk")
            with open(model_out, "w") as json_file:
                json_file.write(model_json)
            print ("Saved best weight to disk, ", score_imed)
            DeepSS_CNN.save_weights(model_weight_out_best)
        print ('The val accuracy is %.5f' % (val_acc))

    print ("Training finished, best validation acc = ",val_acc_best)
    print ("Setting and saving best weights")
    DeepSS_CNN.load_weights(model_weight_out_best)
    DeepSS_CNN.save_weights(model_weight_out)