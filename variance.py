# variance.py
#=======================================================

# Calculate variance of small-scale fields.
# After subtracting their average (stored in AV.py),
# they have zero time-mean, so that we can sum their
# squares at each timestep.

#=======================================================

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import diagnostics
import spectralFilter
import MOMnc

#=======================================================



str1 = 'N513'
str2 = 'relax'
str3 = 'pre-mean_240'

path_orig = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/' + str1 + '/' + str2 + '/'
path = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/' + str1 + '/' + str2 + '/' + str3 + '/'

files_ss = [f for f in os.listdir(path) if f.startswith('ss__')]
files_ls = [f for f in os.listdir(path) if f.startswith('ls__')]
files_ss.sort()
files_ls.sort()

nf = len(files_ss)
fs = 0 
fe = nf

# Read all files
x, y = MOMnc.readVars(path_orig+'prog__0050_359.nc')

L = 3840000.

dx = x[1] - x[0]
dy = y[1] - y[0]

N = len(x)

dx_ = 1./512
dy_ = 1./512

yg, xg = np.mgrid[slice(-1./2,1./2+dy_,dy_),slice(-1./2,1./2+dx_,dx_)]

var = np.load('qvar.npy')
print(np.shape(var))
print(np.shape(xg))

plt.pcolormesh(xg, yg, np.sqrt(var[:,:,0])); plt.colorbar()
plt.xlim(-.5,.5); plt.ylim(-.5,.5);
plt.grid()
plt.text(-.4,-.4,'q st. dev.',fontsize=14,color='w')
plt.show()


q, u, v, h = MOMnc.readFieldsLS(path+files_ss[fs])	

q = np.array(q); u = np.array(u)
v = np.array(v); h = np.array(h)

N, Nx, Nt, Nl = np.shape(q)

AV = np.load(path+'AV.npy')

qav = np.array(AV[0]); uav = np.array(AV[1])
vav = np.array(AV[2]); hav = np.array(AV[3])

for ti in range(0,Nt):
	q[:,:,ti,:] -= qav
	u[:,:,ti,:] -= uav
	v[:,:,ti,:] -= vav
	h[:,:,ti,:] -= hav


# Start calculation of variance.
uvar = np.sum(u**2,axis=2)
vvar = np.sum(v**2,axis=2)
hvar = np.sum(h**2,axis=2)
qvar = np.sum(q**2,axis=2)


quit()

Ttot = Nt
for fi in range(fs+1,fe):
	print(fi)

	q, u, v, h = MOMnc.readFieldsLS(path+files_ss[fi])	

	print('Read data')

	q = np.array(q); u = np.array(u)
	v = np.array(v); h = np.array(h)
	
	N, Nx, Nt, Nl = np.shape(q)	
	Ttot += Nt
	#N = N-1
	for ti in range(0,Nt):
		q[:,:,ti,:] -= np.array(qav)
		u[:,:,ti,:] -= np.array(uav)
		v[:,:,ti,:] -= np.array(vav)
		h[:,:,ti,:] -= np.array(hav)

	# Cumulatively add variances.
	uvar += np.sum(u**2,axis=2)
	vvar += np.sum(v**2,axis=2)
	hvar += np.sum(h**2,axis=2)
	qvar += np.sum(q**2,axis=2)

uvar /= Ttot
vvar /= Ttot
hvar /= Ttot
qvar /= Ttot

np.save('uvar',uvar)
np.save('vvar',vvar)
np.save('hvar',hvar)
np.save('qvar',qvar)
		

plt.contourf(hvar[:,:,0]); plt.colorbar(); plt.show()

quit()


	
uq = np.load('uq.npy')
vq = np.load('vq.npy')
Uq = np.load('Uq.npy')
Vq = np.load('Vq.npy')

qconv = - diagnostics.Div(uq[:,:,0],vq[:,:,0],dx,dy)

lim = 1.e-9
lim2 = .1e-11

qconv = np.where(np.abs(qconv)>lim2,np.sign(qconv)*lim2,qconv)
uq = np.where(np.abs(uq)>lim,np.sign(uq)*lim,uq)
vq = np.where(np.abs(vq)>lim,np.sign(vq)*lim,vq)


plt.contourf(qconv[:,:]); plt.colorbar(); plt.show()

plt.subplot(121)
plt.contourf(uq[:,:,0]); plt.colorbar(); plt.grid()

plt.subplot(122);
plt.contourf(Uq[:,:,0]); plt.colorbar(); plt.grid()

plt.show()

## We need to remember that terms like Uq do not average to zero like they do in my simple SW model.


#q = np.transpose(q,[2,3,0,1])

