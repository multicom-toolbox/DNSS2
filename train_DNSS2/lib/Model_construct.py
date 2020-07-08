# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 21:41:28 2018

@author: Zhiye
"""
from collections import defaultdict
#import cPickle as pickle
import pickle
from PIL import Image

from six.moves import range

import keras.backend as K
from keras.datasets import mnist
from keras.engine.topology import Layer
from keras.layers import Input, Dense, Reshape, Flatten, Embedding, merge, Dropout, Lambda
from keras.layers.advanced_activations import LeakyReLU
from keras.layers.convolutional import UpSampling2D, Convolution2D
from keras.models import Sequential, Model
from keras.optimizers import Adam
from keras.utils.generic_utils import Progbar
from keras.constraints import maxnorm

from keras.layers.recurrent import time_distributed_dense
from keras import activations, initializations
from keras import utils
from keras.engine.topology import Layer
from keras.models import Model
from keras.layers import Activation, Dense, Dropout, Flatten, Input, Merge, MaxPooling1D, AveragePooling1D,UpSampling1D, Convolution1D, LSTM
from keras.layers.normalization import BatchNormalization
from keras.layers.advanced_activations import PReLU
from keras.activations import tanh, softmax
from keras.utils.visualize_util import plot
# from keras.utils import plot_model

# import numpy as np
# import tensorflow as tf
# import datetime
# random_enough_seed = int(datetime.datetime.utcnow().timestamp()*1000)
# np.random.seed(random_enough_seed)

# tf.InteractiveSession()
# tf.set_random_seed(random)

# Helper to build a conv -> BN -> relu block
def _conv_bn_relu1D(nb_filter, nb_row, subsample,use_bias=True):
    def f(input):
        conv = Convolution1D(nb_filter=nb_filter, filter_length=nb_row, subsample_length=subsample,bias=use_bias,
                             init="he_normal", activation='relu', border_mode="same")(input)
        ################## add dropout#################
        norm = BatchNormalization(mode=0, axis=2)(conv)
        return Activation("relu")(norm)
    return f


def _conv_relu1D(nb_filter, nb_row, subsample,use_bias=True):
    def f(input):
        conv = Convolution1D(nb_filter=nb_filter, filter_length=nb_row, subsample_length=subsample,bias=use_bias,
                             init="he_normal", activation='relu', border_mode="same")(input)
        return Activation("relu")(conv)
    return f

# Helper to build a conv -> BN -> softmax block
def _conv_bn_softmax1D(nb_filter, nb_row, subsample,name,use_bias=True):
    def f(input):
        conv = Convolution1D(nb_filter=nb_filter, filter_length=nb_row, subsample_length=subsample,bias=use_bias,
                             init="he_normal", activation='relu', border_mode="same",name="%s_conv" % name)(input)
        norm = BatchNormalization(mode=0, axis=2,name="%s_nor" % name)(conv)
        return Dense(output_dim=8, init="he_normal",name="%s_softmax" % name, activation="softmax")(norm) # change to predict 8 class ss
    
    return f


def _attention_layer(input_dim):
    def f(input):
        attention_probs = Dense(input_dim, activation='softmax')(input)
        attention_mul = merge([input, attention_probs],output_shape=input_dim, mode='mul')
        return attention_mul
    return f
    
def _conv_bn_relu1D_drop(nb_filter, nb_row, subsample,use_bias=True):
    def f(input):
        conv = Convolution1D(nb_filter=nb_filter, filter_length=nb_row, subsample_length=subsample,bias=use_bias,
                             init="he_normal", activation='relu', border_mode="same")(input)
        
        norm = BatchNormalization(mode=0, axis=2)(conv)
        acti = Activation("relu")(norm)
        return Dropout(0.2)(acti)
    return f
    
def DeepCov_SS_with_paras(win_array,feature_num,use_bias,hidden_type,nb_filters,nb_layers,opt):
    DeepSS_input_shape =(None,feature_num)
    filter_sizes=win_array
    DeepSS_input = Input(shape=DeepSS_input_shape)
    DeepSS_convs = []
    for fsz in filter_sizes:
        DeepSS_conv = DeepSS_input
        DeepSS_conv = BatchNormalization(mode=0, axis=2)(DeepSS_conv)
        for i in range(0,nb_layers):
            DeepSS_conv = _conv_bn_relu1D(nb_filter=nb_filters, nb_row=fsz, subsample=1,use_bias=use_bias)(DeepSS_conv)

        DeepSS_conv = _conv_bn_softmax1D(nb_filter=1, nb_row=fsz, subsample=1,use_bias=use_bias,name='local_start')(DeepSS_conv)
        #DeepSS_conv = remove_1d_padding(ktop=ktop_node)(DeepSS_conv) ## remove the padding rows because they don't have targets
        #no need here, because if target is 0, the cross-entropy is zero, error will be not passed
        
        DeepSS_convs.append(DeepSS_conv)
    
    if len(filter_sizes)>1:
        DeepSS_out = Merge(mode='average')(DeepSS_convs)
    else:
        DeepSS_out = DeepSS_convs[0]  
    
    DeepSS_CNN = Model(input=[DeepSS_input], output=DeepSS_out)
    DeepSS_CNN.compile(loss="categorical_crossentropy", metrics=['accuracy'], optimizer=opt)
    # plot(DeepSS_CNN, to_file='model.png')
    DeepSS_CNN.summary()
    return DeepSS_CNN

def RCL_block(input_num_filters,l,fsz):
    out_num_filters = input_num_filters
    
    conv1 = Convolution1D(nb_filter=out_num_filters, filter_length=fsz, border_mode='same')
    stack1 = conv1(l)       
    stack2 = BatchNormalization()(stack1)
    stack3 = Activation("relu")(stack2)
    
    conv2 = Convolution1D(nb_filter=out_num_filters, filter_length=fsz, border_mode='same', init = 'he_normal')
    stack4 = conv2(stack3)
    stack5 = merge([stack1, stack4], mode='sum')
    stack6 = BatchNormalization()(stack5)
    stack7 = Activation("relu")(stack6)
    
    conv3 = Convolution1D(nb_filter=out_num_filters, filter_length=fsz, border_mode='same', weights = conv2.get_weights())
    stack8 = conv3(stack7)
    stack9 = merge([stack1, stack8], mode='sum')
    stack10 = BatchNormalization()(stack9)
    stack11 = Activation("relu")(stack10)    
    
    conv4 = Convolution1D(nb_filter=out_num_filters, filter_length=fsz, border_mode='same', weights = conv2.get_weights())
    stack12 = conv4(stack11)
    stack13 = merge([stack1, stack12], mode='sum')
    stack14 = BatchNormalization(mode=0, axis=2)(stack13)
    stack15 = Activation("relu")(stack14)    
    stack16 = Dropout(0.1)(stack15)
    
    return stack16

def DeepCovRCNN_with_paras(win_array,feature_num,use_bias,hidden_type,nb_filters,nb_layers,opt):
    
    #Build Network
    filter_sizes=win_array
    DeepSS_input_shape =(None,feature_num)
    DeepSS_input = Input(shape=DeepSS_input_shape)
    DeepSS_convs = []
    for fsz in filter_sizes:
        DeepSS_conv = DeepSS_input
        conv_l = Convolution1D(nb_filter=nb_filters, filter_length=fsz, border_mode='same', activation='relu')
        DeepSS_conv = conv_l(DeepSS_conv)
        
        for n in range(nb_layers):
            DeepSS_conv = RCL_block(nb_filters, DeepSS_conv,fsz)
        DeepSS_conv = _conv_bn_softmax1D(nb_filter=1, nb_row=fsz, subsample=1,use_bias=use_bias,name='local_start')(DeepSS_conv)
        
        DeepSS_convs.append(DeepSS_conv)
    
    if len(filter_sizes)>1:
        DeepSS_out = Merge(mode='average')(DeepSS_convs)
    else:
        DeepSS_out = DeepSS_convs[0]  
    DeepSS_RCNN = Model(input=[DeepSS_input], output=DeepSS_out)
    DeepSS_RCNN.compile(loss="categorical_crossentropy", metrics=['accuracy'], optimizer=opt)
    DeepSS_RCNN.summary()
    return DeepSS_RCNN

def identity_Block_deep_noBN(input, nb_filter, kernel_size, with_conv_shortcut=False,use_bias=True, mode='sum'):
    x = _conv_relu1D(nb_filter=nb_filter/4, nb_row=1, subsample=1,use_bias=use_bias)(input)
    x = _conv_relu1D(nb_filter=nb_filter/4, nb_row=kernel_size, subsample=1,use_bias=use_bias)(x)
    x = _conv_relu1D(nb_filter=nb_filter, nb_row=1, subsample=1,use_bias=use_bias)(x)
    if with_conv_shortcut:
        shortcut = _conv_relu1D(nb_filter=nb_filter, nb_row=kernel_size, subsample=1,use_bias=use_bias)(input)
        x = merge([x, shortcut], mode=mode)
        return x
    else:
        x = merge([x, input], mode=mode)
        return x

def identity_Block_deep(input, nb_filter, kernel_size, with_conv_shortcut=False,use_bias=True, mode='sum'):
    x = _conv_relu1D(nb_filter=nb_filter, nb_row=1, subsample=1,use_bias=use_bias)(input)
    x = _conv_relu1D(nb_filter=nb_filter*2, nb_row=kernel_size, subsample=1,use_bias=use_bias)(x)
    x = _conv_relu1D(nb_filter=nb_filter, nb_row=1, subsample=1,use_bias=use_bias)(x)
    if with_conv_shortcut:
        shortcut = _conv_relu1D(nb_filter=nb_filter, nb_row=kernel_size, subsample=1,use_bias=use_bias)(input)
        x = merge([x, shortcut], mode=mode)
        return x
    else:
        x = merge([x, input], mode=mode)
        return x

def DeepResnet1D_with_paras(win_array, feature_num, use_bias, hidden_type, nb_filters, nb_layers, opt):
    filter_sizes = win_array
    DeepSS_input_shape = (None, feature_num)
    DeepSS_input = Input(shape=DeepSS_input_shape)

    DeepSS_convs = []
    for fsz in filter_sizes:
        DeepSS_conv = DeepSS_input
        DeepSS_conv = BatchNormalization(mode=0, axis=2)(DeepSS_conv)
        DeepSS_conv = _conv_relu1D(nb_filter=nb_filters, nb_row=fsz, subsample=1, use_bias=use_bias)(DeepSS_conv)
        # for i in range(0, nb_layers):
        #     DeepSS_conv = identity_Block(DeepSS_conv, nb_filter=nb_filters, kernel_size=fsz,with_conv_shortcut=False,use_bias=True)
        DeepSS_conv = identity_Block_deep(DeepSS_conv, nb_filter=nb_filters, kernel_size=fsz,with_conv_shortcut=False,use_bias=True, mode='sum')
        DeepSS_conv = identity_Block_deep(DeepSS_conv, nb_filter=nb_filters, kernel_size=fsz,with_conv_shortcut=False,use_bias=True, mode='sum')
        DeepSS_conv = identity_Block_deep(DeepSS_conv, nb_filter=nb_filters, kernel_size=fsz,with_conv_shortcut=False,use_bias=True, mode='concat')
        DeepSS_conv = Dropout(0.1)(DeepSS_conv)
        DeepSS_conv = BatchNormalization(mode=0, axis=2)(DeepSS_conv)

        DeepSS_conv = identity_Block_deep(DeepSS_conv, nb_filter=nb_filters*2, kernel_size=fsz,with_conv_shortcut=False,use_bias=True, mode='sum')
        DeepSS_conv = identity_Block_deep(DeepSS_conv, nb_filter=nb_filters*2, kernel_size=fsz,with_conv_shortcut=False,use_bias=True, mode='sum')
        DeepSS_conv = identity_Block_deep(DeepSS_conv, nb_filter=nb_filters*2, kernel_size=fsz,with_conv_shortcut=False,use_bias=True, mode='sum')
        DeepSS_conv = identity_Block_deep(DeepSS_conv, nb_filter=nb_filters*2, kernel_size=fsz,with_conv_shortcut=False,use_bias=True, mode='concat')
        DeepSS_conv = Dropout(0.2)(DeepSS_conv)
        DeepSS_conv = BatchNormalization(mode=0, axis=2)(DeepSS_conv)
        
        DeepSS_conv = identity_Block_deep(DeepSS_conv, nb_filter=nb_filters*4, kernel_size=fsz,with_conv_shortcut=False,use_bias=True, mode='sum')
        DeepSS_conv = identity_Block_deep(DeepSS_conv, nb_filter=nb_filters*4, kernel_size=fsz,with_conv_shortcut=False,use_bias=True, mode='sum')
        DeepSS_conv = identity_Block_deep(DeepSS_conv, nb_filter=nb_filters*4, kernel_size=fsz,with_conv_shortcut=False,use_bias=True, mode='sum')
        DeepSS_conv = identity_Block_deep(DeepSS_conv, nb_filter=nb_filters*4, kernel_size=fsz,with_conv_shortcut=False,use_bias=True, mode='sum')
        DeepSS_conv = identity_Block_deep(DeepSS_conv, nb_filter=nb_filters*4, kernel_size=fsz,with_conv_shortcut=False,use_bias=True, mode='sum')
        DeepSS_conv = identity_Block_deep(DeepSS_conv, nb_filter=nb_filters, kernel_size=fsz,with_conv_shortcut=False,use_bias=True, mode='concat')
        DeepSS_conv = Dropout(0.3)(DeepSS_conv)
        DeepSS_conv = BatchNormalization(mode=0, axis=2)(DeepSS_conv)
        

        DeepSS_conv = _conv_bn_softmax1D(nb_filter=1, nb_row=fsz, subsample=1, use_bias=use_bias, name='local_start')(DeepSS_conv)
        DeepSS_convs.append(DeepSS_conv)

    if len(filter_sizes) > 1:
        DeepSS_out = Merge(mode='average')(DeepSS_convs)
    else:
        DeepSS_out = DeepSS_convs[0]

    DeepSS_RES = Model(input=[DeepSS_input], output=DeepSS_out)
    DeepSS_RES.compile(loss="categorical_crossentropy", metrics=['accuracy'], optimizer=opt)
    DeepSS_RES.summary()
    return DeepSS_RES

def block_inception_a(input,nb_filters,kernel_size,use_bias=True):
    branch_0 = _conv_bn_relu1D(nb_filter=nb_filters, nb_row=kernel_size, subsample=1, use_bias=use_bias)(input)
    branch_1 = _conv_bn_relu1D(nb_filter=nb_filters+8, nb_row=kernel_size, subsample=1, use_bias=use_bias)(input)
    branch_2 = _conv_bn_relu1D(nb_filter=nb_filters+16, nb_row=kernel_size, subsample=1, use_bias=use_bias)(input)
    branch_3 = _conv_bn_relu1D(nb_filter=nb_filters+24, nb_row=kernel_size, subsample=1, use_bias=use_bias)(input)
    #x = concatenate([branch_0, branch_1, branch_2, branch_3], axis=channel_axis)
    net = merge([branch_0, branch_1, branch_2, branch_3], mode='concat')
    return net

def DeepInception1D_with_paras(win_array, feature_num, use_bias, hidden_type, nb_filters, nb_layers, opt):
    filter_sizes = win_array
    DeepSS_input_shape = (None, feature_num)
    DeepSS_input = Input(shape=DeepSS_input_shape)

    DeepSS_convs = []
    for fsz in filter_sizes:
        net = DeepSS_input
        net = _conv_bn_relu1D(nb_filter=nb_filters, nb_row=fsz, subsample=1, use_bias=use_bias)(net)

        ## start inception 1
        branch_0 = _conv_bn_relu1D(nb_filter=nb_filters, nb_row=fsz, subsample=1, use_bias=use_bias)(net)
        branch_1 = _conv_bn_relu1D(nb_filter=nb_filters*2, nb_row=fsz, subsample=1, use_bias=use_bias)(net)
        net = merge([branch_0, branch_1], mode='concat')
        net = Dropout(0.3)(net)#0.3

        ## start inception 2
        branch_0 = _conv_bn_relu1D(nb_filter=nb_filters, nb_row=fsz, subsample=1, use_bias=use_bias)(net)
        branch_0 = _conv_bn_relu1D(nb_filter=nb_filters*2, nb_row=fsz, subsample=1, use_bias=use_bias)(branch_0)

        branch_1 = _conv_bn_relu1D(nb_filter=nb_filters, nb_row=fsz, subsample=1, use_bias=use_bias)(net)
        branch_1 = _conv_bn_relu1D(nb_filter=nb_filters*2, nb_row=1, subsample=1, use_bias=use_bias)(branch_1)
        branch_1 = _conv_bn_relu1D(nb_filter=nb_filters*4, nb_row=fsz, subsample=1, use_bias=use_bias)(branch_1)

        net = merge([branch_0, branch_1], mode='concat')
        net = Dropout(0.3)(net)#0.3

        # 35 x 35 x 384
        # 4 x Inception-A blocks
        for idx in range(nb_layers):
            net = block_inception_a(net,nb_filters,fsz)
            net = Dropout(0.3)(net)#0.3
            
        DeepSS_conv = _conv_bn_softmax1D(nb_filter=1, nb_row=fsz, subsample=1, use_bias=use_bias,
                                             name='local_start')(net)
        DeepSS_convs.append(DeepSS_conv)

    if len(filter_sizes) > 1:
        DeepSS_out = Merge(mode='average')(DeepSS_convs)
    else:
        DeepSS_out = DeepSS_convs[0]

    DeepSS_INCEP = Model(input=[DeepSS_input], output=DeepSS_out)
    DeepSS_INCEP.compile(loss="categorical_crossentropy", metrics=['accuracy'], optimizer=opt)
    DeepSS_INCEP.summary()
    return DeepSS_INCEP

def identity_Block_CRMN(input, nb_filter, kernel_size, with_conv_shortcut=False,use_bias=True, mode='sum'):
    x = _conv_bn_relu1D(nb_filter=nb_filter, nb_row=kernel_size, subsample=1,use_bias=use_bias)(input)
    x = _conv_bn_relu1D(nb_filter=nb_filter, nb_row=kernel_size, subsample=1,use_bias=use_bias)(x)
    if with_conv_shortcut:
        shortcut = _conv_bn_relu1D(nb_filter=nb_filter, nb_row=kernel_size, subsample=1,use_bias=use_bias)(input)
        x = merge([x, shortcut], mode=mode)
        return x
    else:
        x = merge([x, input], mode=mode)
        return x
def DeepCRMN_SS_with_paras(win_array, feature_num, use_bias, hidden_type, nb_filters, nb_layers, opt):

    DeepSS_input_shape = (None, feature_num)
    filter_sizes = win_array
    DeepSS_input = Input(shape=DeepSS_input_shape)
    DeepSS_convs = []
    for fsz in filter_sizes:
        DeepSS_conv = DeepSS_input
        for i in range(0, nb_layers):
            cnn = _conv_bn_relu1D(nb_filter=nb_filters, nb_row=fsz, subsample=1, use_bias=use_bias)(DeepSS_conv)
            res = identity_Block_CRMN(cnn, nb_filter=nb_filters, kernel_size=fsz, with_conv_shortcut=False, use_bias=True)
            cnnres = merge([cnn, res], mode = 'sum')
            cnnres = Dropout(0.2)(cnnres)

            lstm = LSTM(nb_filters, return_sequences=True)(cnnres)
            Lstmlayer = LSTM(nb_filters, return_sequences=True)(lstm)
            DeepSS_conv = cnnres

        DeepSS_conv = merge([DeepSS_conv, Lstmlayer], mode='concat')
        DeepSS_conv = Dropout(0.2)(DeepSS_conv)
        DeepSS_conv = _conv_bn_softmax1D(nb_filter=1, nb_row=fsz, subsample=1, use_bias=use_bias, name='local_start')(DeepSS_conv)
        DeepSS_convs.append(DeepSS_conv)

    if len(filter_sizes) > 1:
        DeepSS_out = Merge(mode='average')(DeepSS_convs)
    else:
        DeepSS_out = DeepSS_convs[0]

    DeepSS_CRMN = Model(input=[DeepSS_input], output=DeepSS_out)
    DeepSS_CRMN.compile(loss="categorical_crossentropy", metrics=['accuracy'], optimizer=opt)
    DeepSS_CRMN.summary()
    return DeepSS_CRMN

def fractal_block(nb_filter, nb_row, subsample=1, use_bias = True):
    def f(input):
        c1 = _conv_relu1D(nb_filter=nb_filter, nb_row=nb_row, subsample=subsample, use_bias=use_bias)(input)
        c2 = _conv_relu1D(nb_filter=nb_filter, nb_row=nb_row, subsample=subsample, use_bias=use_bias)(input)
        c3 = _conv_relu1D(nb_filter=nb_filter, nb_row=nb_row, subsample=subsample, use_bias=use_bias)(input)
        c4 = _conv_relu1D(nb_filter=nb_filter, nb_row=nb_row, subsample=subsample, use_bias=use_bias)(input)
        c4 = _conv_relu1D(nb_filter=nb_filter, nb_row=nb_row, subsample=subsample, use_bias=use_bias)(c4)
        M1 = merge([c3, c4], mode = 'concat')
        c3 = _conv_relu1D(nb_filter=nb_filter, nb_row=nb_row, subsample=subsample, use_bias=use_bias)(M1)
        c4 = _conv_relu1D(nb_filter=nb_filter, nb_row=nb_row, subsample=subsample, use_bias=use_bias)(M1)
        c4 = _conv_relu1D(nb_filter=nb_filter, nb_row=nb_row, subsample=subsample, use_bias=use_bias)(c4)
        M2 = merge([c2, c3, c4], mode = 'concat')
        c2 = _conv_relu1D(nb_filter=nb_filter, nb_row=nb_row, subsample=subsample, use_bias=use_bias)(M2)
        c3 = _conv_relu1D(nb_filter=nb_filter, nb_row=nb_row, subsample=subsample, use_bias=use_bias)(M2)
        c4 = _conv_relu1D(nb_filter=nb_filter, nb_row=nb_row, subsample=subsample, use_bias=use_bias)(M2)
        c4 = _conv_relu1D(nb_filter=nb_filter, nb_row=nb_row, subsample=subsample, use_bias=use_bias)(c4)
        M3 = merge([c3, c4], mode = 'concat')
        c3 = _conv_relu1D(nb_filter=nb_filter, nb_row=nb_row, subsample=subsample, use_bias=use_bias)(M3)
        c4 = _conv_relu1D(nb_filter=nb_filter, nb_row=nb_row, subsample=subsample, use_bias=use_bias)(M3)
        c4 = _conv_relu1D(nb_filter=nb_filter, nb_row=nb_row, subsample=subsample, use_bias=use_bias)(c4)
        M4 = merge([c1, c2, c3, c4], mode = 'concat')
        return M4
    return f

def DeepFracNet_SS_with_paras(win_array, feature_num, use_bias, hidden_type, nb_filters, nb_layers, opt):
    DeepSS_input_shape = (None, feature_num)
    filter_sizes = win_array
    DeepSS_input = Input(shape=DeepSS_input_shape)
    DeepSS_convs = []
    for fsz in filter_sizes:
        DeepSS_conv = DeepSS_input
        for i in range(0, nb_layers):
            DeepSS_conv = fractal_block(nb_filter=nb_filters, nb_row=fsz, subsample=1, use_bias=use_bias)(DeepSS_conv)
            DeepSS_conv = BatchNormalization(mode=0, axis=2)(DeepSS_conv)
            DeepSS_conv = Dropout(0.35)(DeepSS_conv)

        DeepSS_conv = _conv_bn_softmax1D(nb_filter=1, nb_row=fsz, subsample=1, use_bias=use_bias, name='local_start')(DeepSS_conv)
        DeepSS_convs.append(DeepSS_conv)

    if len(filter_sizes) > 1:
        DeepSS_out = Merge(mode='average')(DeepSS_convs)
    else:
        DeepSS_out = DeepSS_convs[0]

    DeepSS_FRAC = Model(input=[DeepSS_input], output=DeepSS_out)
    DeepSS_FRAC.compile(loss="categorical_crossentropy", metrics=['accuracy'], optimizer=opt)
    DeepSS_FRAC.summary()
    return DeepSS_FRAC
