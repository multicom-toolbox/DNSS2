import pylab
# ID: %#b42a2b9c32 
import numpy as np

def montage(X, colormap=pylab.cm.gist_gray, filename=''):

    num_blocks, width = np.shape(X)
   
    mont_width = int(np.ceil(np.sqrt(num_blocks)))
    block_size = np.sqrt(width)

    M = np.zeros((mont_width*block_size+1, mont_width*block_size+1))

    blk_count = 0
    for j in range(mont_width):
        for i in range(mont_width):
            if blk_count >= num_blocks:
                break
           
            sliceM, sliceN = j * block_size, i * block_size
            M[sliceM:sliceM + block_size, sliceN:sliceN + block_size] = np.reshape(X[j*mont_width+i],(block_size,block_size))
            blk_count += 1

    if len(filename) == 0:

      pylab.imshow(M, cmap=colormap)
      pylab.show()

    else:

      pylab.imsave(filename, M)


