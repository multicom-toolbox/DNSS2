import time
# ID: %#2840d6996c 
import numpy as np
import util
import sys

def calc_hidden_probs(data, vh, hb):
   print "Calculating probs. of hidden layer, " + str(data.shape[0]) + " examples."
   probs = np.empty((data.shape[0], vh.shape[1]))

   print "Processing...",
   for i in range(0,data.shape[0]):
       probs[i,:] = 1. / (1 + np.exp(-(np.dot(data[i,:], vh) + hb)))
       if (i+1)%1000 == 0:
           print str(i+1) + "...",
           sys.stdout.flush()

   print "Done"
   return probs

#numpy version of rbm_cudamat
def calc_hidden_probs_batch(data, vh, hb, batchsize):
    """ 
    Calculate the probs in the next layer up given data, and weights

    """

    if data.shape[0] % batchsize != 0:
        print "WARNING!  Batchsize for calc_hidden_probs is not an even divisor of example cnt."

    print "Calculating probs. of hidden layer, " + str(data.shape[0]) + " examples."

    dev_data = data
    lrp_data = np.empty((data.shape[0], vh.shape[1]))
    
    cur_data = np.empty((batchsize, vh.shape[0]))
    nex_data = np.empty((batchsize, vh.shape[1]))

    vishid = vh
    hid_bias = hb
    
    num_batches = data.shape[0]/batchsize
                        
    for batch in range(0, num_batches):      
        row_idx = range(batch*batchsize, (batch+1)*batchsize)      
        cur_data =   dev_data[row_idx,:]
        nex_data = 1. / (1 + np.exp(-(np.dot(cur_data,vishid) + hid_bias)))
        lrp_data[batch*batchsize:(batch+1)*batchsize,:] = nex_data
  
    return lrp_data
