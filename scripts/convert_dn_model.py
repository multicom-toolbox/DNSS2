#python2 convert_dn_model.py  Best_model_toRenzhi0324_9_fea/Arch_10_10_1_TopHid_0_epsilon_0.1_weight_cost_0.0002_im_0.5_fm_0.9_epvb_0.1_ephb_0.1_max_iter_1000/models/DBNmodel_4.dat  Best_model_toRenzhi0324_9_fea/Arch_10_10_1_TopHid_0_epsilon_0.1_weight_cost_0.0002_im_0.5_fm_0.9_epvb_0.1_ephb_0.1_max_iter_1000/models/DBNmodel_final.dat

import sys
GLOBAL_PATH='/space1/jh7x3/DNSS_work/';
sys.path.insert(0, GLOBAL_PATH+'/lib/')

from DN import DN, DN_load 

import numpy as np
import os


if len(sys.argv) != 3:
   sys.stderr.write('Usage: sys.argv[0] Data_prediction(SVM format) model_file addr_output Device')
   print "\n"
   print "For example\npython predict_DN.py Combined.SVM Experiment/Arch_20_10_1_TopHid_0_epsilon_0.00011_weight_cost_0.007_im_0.8_fm_0.9_epvb_0.051_ephb_0.051_max_iter_3/models/DBNmodel_0.dat result/test_output CPU\n"
   sys.exit(1)

model_file = sys.argv[1]  
output_file = sys.argv[2]


dn = DN_load(model_file)

dn.convert_dnss_model(output_file)
