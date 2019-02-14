# plotFluxes.py

#==================================================================

import numpy as np
import matplotlib.pyplot as plt

import diagnostics.diagnostics as dg
import diagnostics.MOMnc as MOMnc
import diagnostics.plotting as plot

#==================================================================

N = 511

Lx = 3840000.
Ly = 3840000.

dx = 1. / (N)
dy = 1. / (N)

yg, xg = np.mgrid[slice(-1./2,1./2+dy,dy),slice(-1./2,1./2+dx,dx)]
y = np.linspace(-1./2,1./2+dy,N)
x = np.linspace(-1./2,1./2+dx,N)

dx = Lx / N
dy = Ly / N

#==================================================================
#==================================================================

field = 'q'
conv = True
fw = '240km'

qpath = '/home/mike/Documents/GulfStream/MOM/Code/mean/' + fw + '/fluxes/q/'
upath = '/home/mike/Documents/GulfStream/MOM/Code/mean/' + fw + '/fluxes/u/'
vpath = '/home/mike/Documents/GulfStream/MOM/Code/mean/' + fw + '/fluxes/v/'
hpath = '/home/mike/Documents/GulfStream/MOM/Code/mean/' + fw + '/fluxes/h/'

#qpath = '/home/mike/Documents/GulfStream/MOM/Code/'
#upath = '/home/mike/Documents/GulfStream/MOM/Code/'
#vpath = '/home/mike/Documents/GulfStream/MOM/Code/'
#hpath = '/home/mike/Documents/GulfStream/MOM/Code/'


mpath = '/home/mike/Documents/GulfStream/MOM/Code/mean/'+fw+'/'

Umean = np.load(mpath+'Umean.npy')
Vmean = np.load(mpath+'Vmean.npy')

U2 = Umean[:,:,0]**2 + Vmean[:,:,0]**2

yd = 100

jet = np.zeros(N)
jeti = np.zeros(N)
for i in range(0,N):
	jeti[i] = int(np.argmax(U2[yd:N-yd,i]) + yd)
	jet[i] = y[jeti[i]]

jet[0:10] = float('nan')
jet[400:N] = float('nan')

#plt.plot(jet); plt.show()

if field == 'q':

	uq = np.load(qpath+'uq.npy')
	vq = np.load(qpath+'vq.npy')
	
	Uq = np.load(qpath+'Uq.npy')
	Vq = np.load(qpath+'Vq.npy')

	uQ = np.load(qpath+'uQ.npy')
	vQ = np.load(qpath+'vQ.npy')

	UQ = np.load(qpath+'UQ.npy')
	VQ = np.load(qpath+'VQ.npy')

	lim1 = 1.e-7 /2
	lim2 = 1.e-7 /2 
	lim3 = 1.e-7 /2
	lim4 = 1.e-7 /2

	str1 = 'uq'
	str2 = 'Uq'
	str3 = 'uQ'
	str4 = 'UQ'

elif field == 'u':

	uq = np.load(upath+'uu.npy')
	vq = np.load(upath+'vu.npy')

	Uq = np.load(upath+'Uu.npy')
	Vq = np.load(upath+'Vu.npy')

	uQ = np.load(upath+'uU.npy')
	vQ = np.load(upath+'vU.npy')

	UQ = np.load(upath+'UU.npy')
	VQ = np.load(upath+'VU.npy')

	lim1 = 1.e-1
	lim2 = 1.e-1
	lim3 = 1.e-1
	lim4 = 1.e-1

	str1 = 'uu'
	str2 = 'Uu'
	str3 = 'uU'
	str4 = 'UU'

elif field == 'v':

	uq = np.load(vpath+'uv.npy')
	vq = np.load(vpath+'vv.npy')

	Uq = np.load(vpath+'Uv.npy')
	Vq = np.load(vpath+'Vv.npy')

	uQ = np.load(vpath+'uV.npy')
	vQ = np.load(vpath+'vV.npy')

	UQ = np.load(vpath+'UV.npy')
	VQ = np.load(vpath+'VV.npy')

	lim1 = 1.e-1
	lim2 = 1.e-1
	lim3 = 1.e-1
	lim4 = 1.e-1

	str1 = 'uv'
	str2 = 'Uv'
	str3 = 'uV'
	str4 = 'UV'

elif field == 'h':

	uq = np.load(hpath+'uh.npy')
	vq = np.load(hpath+'vh.npy')

	Uq = np.load(hpath+'Uh.npy')
	Vq = np.load(hpath+'Vh.npy')

	uQ = np.load(hpath+'uH.npy')
	vQ = np.load(hpath+'vH.npy')

	UQ = np.load(hpath+'UH.npy')
	VQ = np.load(hpath+'VH.npy')

	lim1 = 1.e1
	lim2 = 1.e1
	lim3 = 1.e1
	lim4 = 1.e2

	str1 = 'uh'
	str2 = 'Uh'
	str3 = 'uH'
	str4 = 'UH'



if conv:
	uq = - dg.divLT(uq[:,:,0]*1,vq[:,:,0]*1,dx,dy)
	Uq = - dg.divLT(Uq[:,:,0]*1,Vq[:,:,0]*1,dx,dy)
	uQ = - dg.divLT(uQ[:,:,0]*1,vQ[:,:,0]*1,dx,dy)
	UQ = - dg.divLT(UQ[:,:,0]*1,VQ[:,:,0]*1,dx,dy)

	lim1 = np.max(np.abs(uq)) / 8
	lim2 = lim1
	lim3 = np.max(np.abs(uQ)) / 4
	lim4 = lim3

else:
	
	uq = vq[:,:,0]
	Uq = Vq[:,:,0]
	uQ = vQ[:,:,0]
	UQ = VQ[:,:,0]

	lim1 = np.max(np.max(uq)) / 2
	lim2 = lim1
	lim3 = lim1
	lim4 = np.max(np.max(UQ)) / 2
	

plt.subplot(221)
plt.pcolor(xg,yg,uq,vmin=-lim1,vmax=lim1)
plt.colorbar()
plt.plot(x,jet,linewidth=2.0, color='k')
plt.xlim(np.min(xg[0,]), np.max(xg[0,]))
plt.ylim(np.min(yg[:,0]), np.max(yg[:,0]))
plt.text(0.4,0.4,str1)

plt.subplot(222)
plt.pcolor(xg,yg,Uq,vmin=-lim2,vmax=lim2)
plt.colorbar()
plt.plot(x,jet,linewidth=2.0, color='k')
plt.xlim(np.min(xg[0,]), np.max(xg[0,]))
plt.ylim(np.min(yg[:,0]), np.max(yg[:,0]))
plt.text(0.4,0.4,str2)

plt.subplot(223)
plt.pcolor(xg,yg,uQ,vmin=-lim3,vmax=lim3)
plt.colorbar()
plt.plot(x,jet,linewidth=2.0, color='k')
plt.xlim(np.min(xg[0,]), np.max(xg[0,]))
plt.ylim(np.min(yg[:,0]), np.max(yg[:,0]))
plt.text(0.4,0.4,str3)

plt.subplot(224)
plt.pcolor(xg,yg,UQ,vmin=-lim4,vmax=lim4)
plt.colorbar()
plt.plot(x,jet,linewidth=2.0, color='k')
plt.xlim(np.min(xg[0,]), np.max(xg[0,]))
plt.ylim(np.min(yg[:,0]), np.max(yg[:,0]))
plt.text(0.4,0.4,str4)

plt.tight_layout()
plt.show()

#==

F = uq + uQ + uQ
Flim = np.max(np.abs(F))/2

plt.subplot(121)
plt.pcolor(xg,yg,F,vmin=-Flim,vmax=Flim)
plt.colorbar()
plt.plot(x,jet,linewidth=2.0, color='k')
plt.xlim(np.min(xg[0,]), np.max(xg[0,]))
plt.ylim(np.min(yg[:,0]), np.max(yg[:,0]))
plt.text(0.4,0.4,str1)

plt.subplot(122)
plt.pcolor(xg,yg,UQ,vmin=-lim4,vmax=lim4)
plt.colorbar()
plt.plot(x,jet,linewidth=2.0, color='k')
plt.xlim(np.min(xg[0,]), np.max(xg[0,]))
plt.ylim(np.min(yg[:,0]), np.max(yg[:,0]))
plt.text(0.4,0.4,str2)

plt.tight_layout()
plt.show()

quit()










plt.subplot(221)
plt.contourf(uq[:,:,0]); plt.colorbar()

plt.subplot(222)
plt.contourf(vq[:,:,0]); plt.colorbar()

plt.subplot(223)
plt.contourf(Uq[:,:,0]); plt.colorbar()

plt.subplot(224)
plt.contourf(Vq[:,:,0]); plt.colorbar()

plt.tight_layout()
plt.show()



plt.subplot(221)
plt.contourf(uQ[:,:,0]); plt.colorbar()

plt.subplot(222)
plt.contourf(vQ[:,:,0]); plt.colorbar()

plt.subplot(223)
plt.contourf(UQ[:,:,0]); plt.colorbar()

plt.subplot(224)
plt.contourf(VQ[:,:,0]); plt.colorbar()

plt.tight_layout()
plt.show()
