# geoBal.pys
#==============================================================

# Calculate large-scale velocity field in geostrophic balance with LS field H.
# Load H from LS nc file.

#==============================================================

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import diagnostics.diagnostics as dg
import diagnostics.MOMnc as MOMnc
import diagnostics.plotting as plot


#==============================================================

# Files
str1 = 'N513'
str2 = 'relax'
str3 = 'pre-mean_240'

path_xy = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/' + str1 + '/' + str2 + '/'
path = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/' + str1 + '/' + str2 + '/' + str3 + '/'

path = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/N513/relax/' 

files_ss = [f for f in os.listdir(path) if f.startswith('ss__')]
files.sort()

nf = len(files)
fs = 0
fe = nf

# Domain
x, y = MOMnc.readVars(path_xy+'prog__0050_359.nc')

dx = x[1] - x[0]
dy = y[1] - y[0]

beta = 2.0e-11
f0 = 0.88e-4

# Make sure f is correct. Is it the same as as MOM_input? Does it compare reasonably with Q?

f = f0 + beta * y

plt.plot(f); plt.show()

#==============================================================

# Loop through each file, read LS H data, take appropriate derivatives to calculate Ug, Vg. 

# Given H, geo. bal. implies: Ug = - g ( dH/dy ) / f and Vg = g ( dH/dx ) / f, where H is SSH (not thickness). 

H = MOMnc.readField(path+files_ls[fs],'LSh')

# Redefine H as sum over each layer.
H = np.sum(H, axis=3)

U, V = dg.geoBalUV(H,dx,dy,f,g)









#==============================================================
#==============================================================

















