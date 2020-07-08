# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 21:40:80 2017

@author: Jie Hou
"""
import os
import numpy as np

def chkdirs(fn):
  dn = os.path.dirname(fn)
  if not os.path.exists(dn): os.makedirs(dn)

def chkfiles(fn):
  if os.path.exists(fn):
    return True 
  else:
    return False


def load_train_test_data_for_gan(data_list, feature_dir):
  import pickle
  
  #data_list ="/storage/htc/bdm/Collaboration/jh7x3/DeepCov_SS_SA_project/data/lists/adj_dncon-train.lst"
  #feature_dir='/storage/htc/bdm/Collaboration/jh7x3/DeepCov_SS_SA_project/Deep1Dconv_ss/features_win1/'
  sequence_file=open(data_list,'r').readlines() 
  data_all_dict = []
  seq_num=0;
  print ("######### Loading training data\n\t")
  for i in range(len(sequence_file)):
      pdb_name = sequence_file[i].rstrip()
      print (pdb_name, "..",)
      featurefile = feature_dir + '/' + pdb_name + '.fea'
      if not os.path.isfile(featurefile):
                  print ("feature file not exists: ",featurefile, " pass!")
                  continue           
      
      featuredata = np.loadtxt(featurefile) #(169, 23) # 3+ 20
      data_all_dict.append(featuredata)
      seq_num +=1
  data_all_dict =  np.concatenate(data_all_dict)
  print ("number of sequences: ",seq_num)
  print ("Shape of data_all_dict: ",data_all_dict.shape)
  return data_all_dict



def load_train_test_data_padding_with_interval(data_list, feature_dir,Interval,seq_end,prefix):
  import pickle
  
  #data_list ="/storage/htc/bdm/Collaboration/jh7x3/DeepCov_SS_SA_project/data/lists/adj_dncon-train.lst"
  #feature_dir='/storage/htc/bdm/Collaboration/jh7x3/DeepCov_SS_SA_project/Deep1Dconv_ss/features_win1/'
  sequence_file=open(data_list,'r').readlines() 
  data_all_dict = dict()
  print ("######### Loading training data\n\t")
  for i in range(len(sequence_file)):
      pdb_name = sequence_file[i].rstrip()
      # print (pdb_name, "..",)
      featurefile = feature_dir + '/' + pdb_name + '.fea'
      if not os.path.isfile(featurefile):
                  print ("feature file not exists: ",featurefile, " pass!")
                  continue           
      
      featuredata = np.loadtxt(featurefile) #(169, 51)
      
      fea_len = featuredata.shape[0]
      train_labels = featuredata[:,0:3]#(169, 3)
      train_feature = featuredata[:,3:] #(169, 48)
      for ran in range(0,seq_end,Interval):
          start_ran = ran
          end_ran = ran + Interval
          if end_ran > seq_end:
              end_ran = seq_end 
          if fea_len >start_ran and   fea_len <= end_ran:
              featuredata_all_pad = np.zeros((end_ran,featuredata.shape[1])) ## padd all zeros, including the target, not [1,0,0],  now is [0,0,0], but will be not used in the backproporgation
              featuredata_all_pad[:featuredata.shape[0],:featuredata.shape[1]] = featuredata
              
              #print "fea_len: ",fea_len
              fea_len_new=end_ran
              if fea_len_new in data_all_dict:
                  data_all_dict[fea_len_new].append(featuredata_all_pad)
              else:
                  data_all_dict[fea_len_new]=[]
                  data_all_dict[fea_len_new].append(featuredata_all_pad)               
          else:
              continue
  # list(data_all_dict.keys()).
  for key in data_all_dict.keys():
      myarray = np.asarray(data_all_dict[key])
      data_all_dict[key] = myarray.reshape(len(myarray),myarray.shape[1],myarray.shape[2])
      print ("keys: ", key, " shape: ", data_all_dict[key].shape)
  return data_all_dict

def load_train_test_data_padding_with_interval_8class(data_list, feature_dir,Interval,seq_end,prefix):
  import pickle
  #data_list ="/storage/htc/bdm/Collaboration/jh7x3/DeepCov_SS_SA_project/data/lists/adj_dncon-train.lst"
  #feature_dir='/storage/htc/bdm/Collaboration/jh7x3/DeepCov_SS_SA_project/Deep1Dconv_ss/features_win1/'
  sequence_file=open(data_list,'r').readlines() 
  data_all_dict = dict()
  print ("######### Loading training data\n\t")
  for i in range(len(sequence_file)):
      pdb_name = sequence_file[i].rstrip()
      # print (pdb_name, "..",)
      featurefile = feature_dir + '/' + pdb_name + '.fea'
      if not os.path.isfile(featurefile):
                  print ("feature file not exists: ",featurefile, " pass!")
                  continue           
      
      featuredata = np.loadtxt(featurefile) #(169, 51)
      
      fea_len = featuredata.shape[0]
      train_labels = featuredata[:,0:8]#(169, 8)
      train_feature = featuredata[:,8:] #(169, 48)
      for ran in range(0,seq_end,Interval):
          start_ran = ran
          end_ran = ran + Interval
          if end_ran > seq_end:
              end_ran = seq_end 
          if fea_len >start_ran and   fea_len <= end_ran:
              featuredata_all_pad = np.zeros((end_ran,featuredata.shape[1])) ## padd all zeros, including the target, not [1,0,0],  now is [0,0,0], but will be not used in the backproporgation
              featuredata_all_pad[:featuredata.shape[0],:featuredata.shape[1]] = featuredata
              
              #print "fea_len: ",fea_len
              fea_len_new=end_ran
              if fea_len_new in data_all_dict:
                  data_all_dict[fea_len_new].append(featuredata_all_pad)
              else:
                  data_all_dict[fea_len_new]=[]
                  data_all_dict[fea_len_new].append(featuredata_all_pad)               
          else:
              continue
  # list(data_all_dict.keys()).
  for key in data_all_dict.keys():
      myarray = np.asarray(data_all_dict[key])
      data_all_dict[key] = myarray.reshape(len(myarray),myarray.shape[1],myarray.shape[2])
      print ("keys: ", key, " shape: ", data_all_dict[key].shape)
  return data_all_dict