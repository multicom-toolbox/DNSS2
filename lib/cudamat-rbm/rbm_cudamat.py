import time
import numpy as np
from itertools import chain
import cudamat as cm
import asciihist as ascii

def calc_hidden_probs(data, vh, hb, batchsize):
    """ 
    Calculate the probs in the next layer up given data, and weights

    """

    if data.shape[0] % batchsize != 0:
        print "WARNING!  Batchsize for calc_hidden_probs is not an even divisor of example cnt."
    

    print "Calculating probs. of hidden layer, " + str(data.shape[0]) + " examples."

    dev_data = cm.CUDAMatrix(data)
    lrp_data = np.empty((data.shape[0], vh.shape[1]))

    cur_data = cm.empty((batchsize, vh.shape[0]))
    nex_data = cm.empty((batchsize, vh.shape[1]))

    vishid = cm.CUDAMatrix(vh)
    hid_bias = cm.CUDAMatrix(hb)
    
    num_batches = data.shape[0]/batchsize

    for batch in range(0, num_batches):
        cur_data = dev_data.get_row_slice(batch*batchsize, (batch+1)*batchsize)      
        cm.dot(cur_data, vishid, target = nex_data)
        nex_data.add_row_vec(hid_bias)
        nex_data.apply_sigmoid()

        nex_data.copy_to_host()
        lrp_data[batch*batchsize:(batch+1)*batchsize,:] = nex_data.asarray().copy()
  
    return lrp_data

