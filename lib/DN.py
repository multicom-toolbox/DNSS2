
import numpy as np
import pickle
import util
import gzip
class DN:
    """ Represents a Deep Network and related training and test procedures. 

    Attributes:
        layer_count    The number of layers in the network, including input
        weights        An array of weights for each layer based on node type
        params         An array of learning parameters for each layer 
        arch           An array of node_count and node_type for each layer

        Note that the weights and parameters begin indexing at 1, arch from 0
        with zero being the data/input layer.


    TODO:
       LEGACY - Set card to use
       Functions to save/load a model to a file in a binary format
       Training and testing functions that use an abstracted layer/interface
       Functions to reinitialize learning parameters based on architecture
       Function to randomly initialize weights
       Function to do adjust using backprob (for finetuning or learning NN)

    """ 

    layer_count = 0       
    weights = []
    params = []
    arch = []
    legacy_card_number = 0     

    def __init__(self, arch_specifier):
        """ Initialize a DN from an architecture specification. 

        The architecture specifier is a list with node counts by layer and
        optionally followed by node type ('S': sigmoid (default), '-': data).
        The learning parameters will also be set to default values based on 
        the node type.  The last layer should be the size of the target or
        target vector.


        """


        self.layer_count = 1
        self.arch=[]
        self.weights=[]
        self.params=[]
        self.arch.append({'node_type': '-', 'node_count': 0})
        self.weights.append(dict())  # Set with dummy value
        self.params.append(dict())   # Set with dummy value
        for value in arch_specifier:
            if type(value) is int:
                self.arch.append({'node_type': 'S', 'node_count': value})
                self.layer_count = self.layer_count + 1
            if type(value) is str:
                if value.upper() == 'S':
                    if self.layer_count > 0:
                        self.arch[(self.layer_count-1)]['node_type'] = value.upper()

                else:
                    print "Invalid node type encountered (%s).  Will " % value,
                    print "use a default node type for this layer ",
                    print "(%d)" % self.layer_count

        print "\n-----------------------------------------"
        for i in range(1, self.layer_count):  
            print "Setting up weights for layer %d," % i,
            print "node type is %s." % self.arch[i]['node_type']
            if self.arch[i]['node_type'] == 'S':
                self.weights.append({'w': 0, 'vb': 0, 'hb': 0})
                self.params.append({'batch_size': 1000, 'num_epochs': 20, 
                                    'final_momentum': 0.9, 
                                    'initial_momentum': 0.5, 
                                    'weight_cost': 0.0002, 'epsilon': 0.1 })
        print "-------------------------------------------\n"
    
    def calc_output_legacy(self, data, batch_size):
        """ Calculate the output (probababilies) for a set of data

        The purpose of this function is to calculate the output of a DN on 
        some set of data.  The values will calculated using rbm_cudamat
        on slices of data specified by the batch size

        """ 

        import cudamat as cm
        import rbm_numpy, rbm_cudamat

        # Initialize CUDA
        cm.cublas_init()
        cm.CUDAMatrix.init_random(1)

        if self.legacy_card_number != 0:
            cm.cuda_set_device(self.legacy_card_number)    

        # Create output, use the size of the last layer to do this
        output = np.empty((data.shape[0], 
                           self.arch[(self.layer_count-1)]['node_count']))

        # Slice up data, handling batches of batch_size. USE INT DIVISION
        processed = 0
        for j in range(data.shape[0] // batch_size):

            curr_data = data[j*batch_size:(j+1)*batch_size,:]

            for i in range(1, self.layer_count):

                # Handle a sigmoid node
                if self.arch[i]['node_type'] == 'S':
                    curr_data = \
                      rbm_cudamat.calc_hidden_probs(curr_data, 
                                                    self.weights[i]['w'], 
                                                    self.weights[i]['hb'], 
                                                    batch_size) 

            output[j*batch_size:(j+1)*batch_size,:] = curr_data[:,:]
            processed = processed + batch_size


        # Now handle anything that was left over i.e., what didn't fit in
        if processed != data.shape[0]:
            
            curr_data = data[processed:,:]

            for i in range(1, self.layer_count):

                # Handle a sigmoid node
                if self.arch[i]['node_type'] == 'S':
                    curr_data = \
                      rbm_numpy.calc_hidden_probs(curr_data, 
                                                  self.weights[i]['w'], 
                                                  self.weights[i]['hb'])

            output[processed:,:] = curr_data[:,:]

        cm.cublas_shutdown()

        return output 
    
    
    def convert_dnss_model(self,new_model):
        dict= {'l1_vh': self.weights[1]['w'], 'l1_hb': self.weights[1]['hb'],'l2_vh': self.weights[2]['w'], 'l2_hb': self.weights[2]['hb'],'l3_vh': self.weights[3]['w'], 'l3_hb': self.weights[3]['hb'],'l4_vh': self.weights[4]['w'], 'l4_hb': self.weights[4]['hb'],'l5_vh': self.weights[5]['w'], 'l5_hb': self.weights[5]['hb']}
        varlist = 'l1_vh l1_hb l2_vh l2_hb l3_vh l3_hb l4_vh l4_hb l5_vh l5_hb'
        util.save(new_model,varlist,dict)


def DN_load(filename):
    """ Load the deep network and associated architecture and layer count

    The deep network and values needed for classification are loaded from
    the specified file location.  This is done using the pickle module.
    Note the learning parameters may or may not be not saved.

    Usage: dn = DN_load('model.dat')

    """

    fh = open(filename, 'r')
    obj = pickle.load(fh)
        
    return obj

# TODO: do a check for the learning parameters and re-initalize if need be
def loadmodel(fname, target_dict, verbose = True):
    fo = gzip.GzipFile(fname, 'rb')
    var_list = pickle.load(fo)
    if verbose:
        print "Load architecture: ", var_list
    for var in var_list:
        target_dict[var] = pickle.load(fo)
    fo.close()