# divRV.py

#==============================================================

import os

import numpy as np
import matplotlib.pyplot as plt

import MOMnc
import diagnostics

#==============================================================

str1 = 'N513'
str2 = 'relax'
str3 = 'pre-mean'
path_orig = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/' + str1 + '/' + str2 + '/'

path = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/N513/relax/' + str3 + '/' 

files = [f for f in os.listdir(path) if f.startswith('ss__')]
files.sort()

nf = len(files)
fi = 28

# Read all files
x, y = MOMnc.readVars(path_orig+'prog__0050_359.nc')
dx = x[1] - x[0]; dy = y[1] - y[0]

q, u, v, h = MOMnc.readFieldsSS(path+files[fi])

AV = np.load('AV.npy')

qav = AV[0]; uav = AV[1]; vav = AV[2]; hav = AV[3]

#q = np.transpose(q,[2,3,0,1])

#=======================================================

N, Nx, Nt, Nl = np.shape(q)

N = N-1

#q = q[0:N,0:N,:,:]

Lx = 3840000.
Ly = 3840000.

dx = 1. / (N)
dy = 1. / (N)

yg, xg = np.mgrid[slice(-1./2,1./2+dy,dy),slice(-1./2,1./2+dx,dx)]
yg_low, xg_low = np.mgrid[slice(-1./2,1./2+dy*4,4*dy),slice(-1./2,1./2+4*dx,4*dx)]

#=======================================================


u_ = u[:,:,0,:] - uav[:,:,:] 
v_ = v[:,:,0,:] - vav[:,:,:]
h_ = h[:,:,0,:] - hav[:,:,:]
q_ = q[:,:,0,:] - qav[:,:,:]

RV = np.zeros((N,N,Nt,Nl))
Div = np.zeros((N,N,Nt,Nl))
for li in range(0,3):
	RV[:,:,0,li] = diagnostics.RV(u_[0:N,0:N,li],v_[0:N,0:N,li],dx,dy)
	Div[:,:,0,li] = diagnostics.Div(u_[0:N,0:N,li],v_[0:N,0:N,li],dx,dy)



#=======================================================

cm = 'bwr'
fs = 16

if True:
	for li in range(0,3):

		RVlim = 0.8*np.max(np.abs(RV[1:N-1,1:N-1,0,li]))
		Divlim = 0.8*np.max(np.abs(Div[1:N-1,1:N-1,0,li]))
		hlim = 0.8*np.max(np.abs(h_[1:N-1,1:N-1,li]))
		
		plt.figure(figsize=[16,5])

		plt.subplot(131)
		plt.pcolormesh(xg,yg,RV[:,:,0,li],cmap=cm,vmin=-RVlim,vmax=RVlim)
		plt.colorbar()
		plt.text(-0.4,0.4,'RV',fontsize=fs)
		plt.axis([xg.min(), xg.max(), yg.min(), yg.max()])
		plt.grid()

		plt.subplot(132)
		plt.pcolormesh(xg,yg,Div[:,:,0,li],cmap=cm,vmin=-Divlim,vmax=Divlim)
		plt.colorbar()
		plt.text(-0.4,0.4,'Div',fontsize=fs)
		plt.axis([xg.min(), xg.max(), yg.min(), yg.max()])
		plt.grid()

		plt.subplot(133)
		plt.pcolormesh(xg,yg,h_[:,:,li],cmap=cm,vmin=-hlim,vmax=hlim)
		plt.colorbar()
		plt.text(-0.4,0.4,'h',fontsize=fs)
		plt.axis([xg.min(), xg.max(), yg.min(), yg.max()])
		plt.grid()

		plt.tight_layout()
		#plt.show()
		plt.savefig('plots/all'+str(li+1))

		plt.close()

		qlim = np.max(np.abs(q_[:,:,li]))
		plt.pcolormesh(xg,yg,q_[:,:,li],cmap=cm,vmin=-qlim,vmax=qlim)
		plt.colorbar()
		plt.text(-0.4,0.4,'q',fontsize=fs)
		plt.axis([xg.min(), xg.max(), yg.min(), yg.max()])
		plt.grid()
		#plt.show()
		plt.savefig('plots/q'+str(li+1))


		plt.close()
		


