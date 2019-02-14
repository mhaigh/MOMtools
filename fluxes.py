# fluxes.py
#=======================================================

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from diagnostics import diagnostics, MOMnc

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

dx = x[1] - x[0]
dy = y[1] - y[0]
print(dx)

q, u, v, h = MOMnc.readFieldsLS(path+files_ss[fs])	
Q, U, V, H = MOMnc.readFieldsLS(path+files_ls[fs])

q = np.array(q); u = np.array(u)
v = np.array(v); h = np.array(h)
Q = np.array(Q); U = np.array(U)
V = np.array(V); H = np.array(H)

N, Nx, Nt, Nl = np.shape(q)

AV = np.load(path+'AV.npy')

qav = AV[0]; uav = AV[1]; vav = AV[2]; hav = AV[3]

for ti in range(0,Nt):
	q[:,:,ti,:] -= np.array(qav)
	u[:,:,ti,:] -= np.array(uav)
	v[:,:,ti,:] -= np.array(vav)
	h[:,:,ti,:] -= np.array(hav)

	Q[:,:,ti,:] += np.array(qav)
	U[:,:,ti,:] += np.array(uav)
	V[:,:,ti,:] += np.array(vav)
	H[:,:,ti,:] += np.array(hav)

Umean = np.sum(U,axis=2)
Vmean = np.sum(V,axis=2)
Hmean = np.sum(Q,axis=2)
Qmean = np.sum(H,axis=2)

# Need to change the last two inputs to this function
uq, vq, Uq, Vq, uQ, vQ, UQ, VQ = diagnostics.flux(u,U,v,V,q,Q)

plt.contourf(uq[:,:,0]); plt.show()

qt = np.sum(q,axis=2)


Ttot = Nt
for fi in range(fs+1,fe):
	print(fi)

	q, u, v, h = MOMnc.readFieldsLS(path+files_ss[fi])	
	Q, U, V, H = MOMnc.readFieldsLS(path+files_ls[fi])	

	print('Read data')

	q = np.array(q); u = np.array(u)
	v = np.array(v); h = np.array(h)
	Q = np.array(Q); U = np.array(U)
	V = np.array(V); H = np.array(H)
	
	N, Nx, Nt, Nl = np.shape(Q)	
	Ttot += Nt
	#N = N-1
	for ti in range(0,Nt):
		q[:,:,ti,:] -= np.array(qav)
		u[:,:,ti,:] -= np.array(uav)
		v[:,:,ti,:] -= np.array(vav)
		h[:,:,ti,:] -= np.array(hav)

		Q[:,:,ti,:] += np.array(qav)
		U[:,:,ti,:] += np.array(uav)
		V[:,:,ti,:] += np.array(vav)
		H[:,:,ti,:] += np.array(hav)

	Umean += np.sum(U,axis=2)
	Vmean += np.sum(V,axis=2)
	Hmean += np.sum(H,axis=2)
	Qmean += np.sum(Q,axis=2)


	qt += np.sum(q, axis=2)

	# Need to change the last two inputs to this function
	uf, vf, Uf, Vf, uF, vF, UF, VF = diagnostics.flux(u,U,v,V,q,Q)
	uq += uf
	vq += vf	
	Uq += Uf
	Vq += Vf
	uQ += uF
	vQ += vF	
	UQ += UF
	VQ += VF


Umean /= Ttot
Vmean /= Ttot
Hmean /= Ttot
Qmean /= Ttot

np.save('Umean',Umean)
np.save('Vmean',Vmean)
np.save('Hmean',Hmean)
np.save('Qmean',Qmean)
		
qt /= Ttot

plt.contourf(qt[:,:,0]); plt.colorbar(); plt.show()
uq /= Ttot
vq /= Ttot
Uq /= Ttot
Vq /= Ttot
uQ /= Ttot
vQ /= Ttot
UQ /= Ttot
VQ /= Ttot


np.save('uq',uq)
np.save('vq',vq)
np.save('Uq',Uq)
np.save('Vq',Vq)
np.save('uQ',uQ)
np.save('vQ',vQ)
np.save('UQ',UQ)
np.save('VQ',VQ)

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

