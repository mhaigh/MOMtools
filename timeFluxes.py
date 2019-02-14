# timeFluxes.py
#=======================================================

# Output time-dependent sample of fluxes or convergences for animation.

#=======================================================

import os
import sys

import numpy as np

import diagnostics
import MOMnc

#=======================================================

def fluxes(path):

	files_ss = [f for f in os.listdir(path) if f.startswith('ss__')]
	files_ls = [f for f in os.listdir(path) if f.startswith('ls__')]
	files_ss.sort()
	files_ls.sort()

	nf = len(files_ss)
	fs = nf-1
	fe = nf

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

	Ttot = Nt

	# Loop through remaining files.
	for fi in range(fs+1,fe):
		print(fi)
	
		qtmp, utmp, vtmp, htmp = MOMnc.readFieldsLS(path+files_ss[fi])	
		Qtmp, Utmp, Vtmp, Htmp = MOMnc.readFieldsLS(path+files_ls[fi])	

		print('Read data')

		qtmp = np.array(qtmp); utmp = np.array(utmp); vtmp = np.array(vtmp); htmp = np.array(htmp)
		Qtmp = np.array(Qtmp); Utmp = np.array(Utmp); Qtmp = np.array(Vtmp); Htmp = np.array(Htmp)
	
		N, Nx, Nt, Nl = np.shape(Qtmp)	
		Ttot += Nt
		#N = N-1
		for ti in range(0,Nt):
			qtmp[:,:,ti,:] -= np.array(qav)
			utmp[:,:,ti,:] -= np.array(uav)
			vtmp[:,:,ti,:] -= np.array(vav)
			htmp[:,:,ti,:] -= np.array(hav)

			Qtmp[:,:,ti,:] += np.array(qav); Utmp[:,:,ti,:] += np.array(uav)
			Vtmp[:,:,ti,:] += np.array(vav); Htmp[:,:,ti,:] += np.array(hav)

		u = np.concatenate((u,utmp),axis=2)
		v = np.concatenate((v,vtmp),axis=2)
		q = np.concatenate((q,qtmp),axis=2)
		U = np.concatenate((U,Utmp),axis=2)
		V = np.concatenate((V,Vtmp),axis=2)
		Q = np.concatenate((Q,Qtmp),axis=2)

	# Make fluxes
	uf = (u*h)[:,:,:,0]# + (U*h)[:,:,:,0] + (u*H)[:,:,:,0]
	vf = (v*h)[:,:,:,0]# + (V*h)[:,:,:,0] + (v*H)[:,:,:,0]
	#uq = (U*Q)[:,:,:,0] 
	#vq = (V*Q)[:,:,:,0]	
	
	dx = 3840000. / N
	conv = np.zeros((N,N,Nt))
	for ti in range(0,Nt):
		conv[:,:,ti] = -diagnostics.Div(uf[:,:,ti],vf[:,:,ti],dx,dx)
		

	return conv



