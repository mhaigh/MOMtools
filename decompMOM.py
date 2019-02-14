# decompMOM.py
#=======================================================

# At every snapshot in time, this module decomposes the PV, velocities and thicknesses into large-scale and small-scale components (saving the LS component only)
# Removing the time-averaged part of the small-scale flow and adding it to the large-scale flow is done in a separate module.

#=======================================================

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

from diagnostics.filters import spectralFilter, windowFilter
from diagnostics import MOMnc

#=======================================================

str1 = 'N513'
str2 = 'no_relax'

path = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/' + str1 + '/' + str2 + '/'
path_out = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/' + str1 + '/' + str2 + '/pre-mean_240/'

files = [f for f in os.listdir(path) if f.startswith('prog__')]
files.sort()
Nf = len(files)

x, y = MOMnc.readVars(path+files[0])

print('Read from ' + files[22])

N = len(x)

Lx = 3840000.
Ly = 3840000.

dx = 1. / (N-1)
dy = 1. / (N-1)

yg, xg = np.mgrid[slice(-1./2,1./2+dy,dy),slice(-1./2,1./2+dx,dx)]

ts = 0

# Filter width (m)
fw = 240.e3
str3 = str(int(fw/1000))

#=======================================================

fs = 22 # file 22 where we begin decomposition from

fs = fs
fe = fs+1
for fi in range(fs,fe):
	print(fi)
	q, u, v, h = MOMnc.readFields(path+files[fi])


	print(np.shape(q))	
	
	q = np.transpose(q,[2,3,0,1])
	u = np.transpose(u,[2,3,0,1])
	v = np.transpose(v,[2,3,0,1])
	h = np.transpose(h,[2,3,0,1])

	N, Nx, Nt, Nl = np.shape(q)
	
	N = N-1

	q = q[0:N,0:N,:,:]
	u = u[0:N,0:N,:,:]
	v = v[0:N,0:N,:,:]
	h = h[0:N,0:N,:,:]	

	#=======================================================

	# Initialise small-scale and large-scale flow components
	SSq = np.zeros((N,N,Nt,Nl)); SSu = np.zeros((N,N,Nt,Nl))
	SSv = np.zeros((N,N,Nt,Nl)); SSh = np.zeros((N,N,Nt,Nl))
	LSq = np.zeros((N,N,Nt,Nl)); LSu = np.zeros((N,N,Nt,Nl))
	LSv = np.zeros((N,N,Nt,Nl)); LSh = np.zeros((N,N,Nt,Nl))

	print('Decompose q')
	SSq, LSq = windowFilter(q,fw)
	print('Decompose u')
	SSu, LSu = windowFilter(u,fw)
	print('Decompose v')
	SSv, LSv = windowFilter(v,fw)
	print('Decompose h')
	SSh, LSh = windowFilter(h,fw)
	print(' ')

	LS = [LSq,LSu,LSv,LSh]
	SS = [SSq,SSu,SSv,SSh]
	MOMnc.writeFields(path_out,files[fi],LS,'ls')
	MOMnc.writeFields(path_out,files[fi],SS,'ss')


#=======================================================

sys.exit()

fi = 30
SS = MOMnc.readField(path_out + 'ss' + files[fi][4::])
SSq = SS[0]

q, u, v, h = MOMnc.readFields(path+files[fi])
q = np.transpose(q,[2,3,0,1])
q = q[0:N-1,0:N-1,:,:]

LSq = q[:,:,:,0] - SSq[:,:,:,0]

#sys.exit()

cm = 'seismic'
fs = 18
ss_lim = np.max(np.abs(SSq[:,:,0]))

plt.figure(figsize=[19,6])
plt.subplot(131)
plt.pcolormesh(xg,yg,q[:,:,0,0],cmap=cm,vmin=1.5e-7,vmax=6.e-7)
plt.colorbar()
plt.text(-0.4,0.4,'Full flow',fontsize=fs)
plt.axis([xg.min(), xg.max(), yg.min(), yg.max()])
plt.grid()

plt.subplot(132)
plt.pcolormesh(xg,yg,LSq[:,:,0],cmap=cm)
plt.colorbar()
plt.text(-0.4,0.4,'LS',fontsize=fs)
plt.axis([xg.min(), xg.max(), yg.min(), yg.max()])
plt.grid()

plt.subplot(133)
plt.pcolormesh(xg,yg,SSq[:,:,0,0],cmap=cm,vmin=-ss_lim,vmax=ss_lim)
plt.colorbar()
plt.axis([xg.min(), xg.max(), yg.min(), yg.max()])
plt.text(-0.4,0.4,'SS',fontsize=fs)
plt.grid()

plt.tight_layout()
plt.savefig('decomp_' + str1 + '_' + str2 + '_' + str3 + '.png')
plt.show()

sys.exit()














