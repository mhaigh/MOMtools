# avSSfield.py

#=======================================================

# A function for time-averaging small-scale field, and shifting the 
# time-average over to the large-scale field

#=======================================================

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

from filters import spectralFilter, windowFilter
import MOMnc

#=======================================================

str1 = 'N513'
str2 = 'relax'

#path = '~/cluster/MOM3/'
path = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/' + str1 + '/' + str2 + '/pre-mean_240/'
#path = './'
path_orig = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/' + str1 + '/' + str2 + '/'

files = [f for f in os.listdir(path) if f.startswith('ss__')]
files.sort()

Nf = len(files)

N = 512

qav = np.zeros((N,N,3)); uav = np.zeros((N,N,3)) 	
vav = np.zeros((N,N,3)); hav = np.zeros((N,N,3))

Nt_tot = 0

#=======================================================

for fi in range(0,Nf):
	# For each file in decomposition
	file_ = files[fi]
	print(file_)
	file_id = file_[2::]

	SS = MOMnc.readField(path + 'ss' + file_id)
	
	SSq = np.array(SS[0]); SSu = np.array(SS[1])
	SSv = np.array(SS[2]); SSh = np.array(SS[3])
	
	N, Nx, Nt, Nl = np.shape(SSq)

	Nt_tot += Nt

	qav += np.sum(SSq,2); uav += np.sum(SSu,2)
	vav += np.sum(SSv,2); hav += np.sum(SSh,2)

qav = qav / Nt_tot; uav = uav / Nt_tot
vav = vav / Nt_tot; hav = hav / Nt_tot

AV = np.array([qav,uav,vav,hav])
np.save('AV.npy',AV)



AV = np.load('AV.npy')
qav = AV[0]
print(np.shape(qav))

sum_ = np.zeros((N,N,3))
for ti in range(0,Nt):
	SSq[:,:,ti,:] = np.array(SSq[:,:,ti,:]) - qav
	sum_ += SSq[:,:,ti,:]


sum_ = sum_ / Nt


plt.contourf(qav[:,:,0])
plt.colorbar()
plt.show()
	
plt.contourf(sum_[:,:,0])
plt.colorbar()
plt.show()
		















