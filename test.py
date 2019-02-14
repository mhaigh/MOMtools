import numpy as np
import matplotlib.pyplot as plt

import os

from diagnostics import diagnostics as dg
from diagnostics import MOMnc





path = 'test/'

uh = np.load(path+'uh.npy')
vh = np.load(path+'vh.npy')

uhDiv = np.load(path+'uhDiv.npy')

var1 = np.load(path+'var_uhDiv.npy')

var2 = np.load(path+'var_H.npy')

cov = np.load(path+'cov.npy')

var1inv = 1. / np.sqrt(var1)
var2inv = 1. / np.sqrt(var2)

lim1 = 200

var1inv = np.where(var1inv>lim1, lim1, var1inv)

corr = cov * var1inv * var2inv


#var1 = dg.smooth(var1)


plt.subplot(221)
plt.contourf(var1inv[:,:,0]); plt.colorbar();

plt.subplot(222)
plt.contourf(var2inv[:,:,0]); plt.colorbar();

plt.subplot(223)
plt.contourf(corr[:,:,0]); plt.colorbar()

plt.subplot(224)
plt.contourf(dg.smooth(corr[:,:,0],2)); plt.colorbar()

plt.tight_layout()
plt.show()
