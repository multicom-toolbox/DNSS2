import minimize as mi
# ID: %#6429ee5728 
import check_grad as cg
import numpy as np


min = np.zeros(4)

r_val = mi.minimize(min, cg.test_f, cg.test_grad, [], maxnumlinesearch=80)
min = r_val[0]

print min
print r_val[1]
