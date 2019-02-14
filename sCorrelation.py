# sCorrelation.py

#==================================================================

# Calculatute lagged autocorrelation and horiztonal correlation in field.

#==================================================================

import os 

import numpy as np
import matplotlib.pyplot as plt

import diagnostics
import MOMnc

#==================================================================

# How does it work?

# Spatial correlation first.
# Select a grid point x0, at which we want its spatial correlation statistics.
# We need the variance (st. dev.) of the field at ALL grid points x, which we may already have.
# Take the covariance of the field evaluated at x0 with same field evaluated at x, divide by their standard deviations.

# If finding corr for fluxes, we can take the mean in saved numpy files.

#==================================================================

def covloop(cov,f,x0,y0,N):

	for i in range(0,N):
		for j in range(0,N):
			cov[j,i,] += np.sum((f[y0,x0,] * f[j,i,]))

	return cov

#==================================================================

str1 = 'N513'
str2 = 'relax'
str3 = 'pre-mean_240'

path = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/' + str1 + '/' + str2 + '/' + str3 + '/'

files_ss = [f for f in os.listdir(path) if f.startswith('ss__')]
files_ss.sort()

nf = len(files_ss)
fs = 0 
fe = nf

dx = 7.5; dy = 7.5


x0 = 512-30
y0 = 300



#==================================================================

# Read the first field.

# Only read one field.
f = np.array(MOMnc.readField(path+files_ss[fs],'LSh'))

N, Nx, Nt, Nl = np.shape(f)
Ntot = Nt

AV = np.load(path+'AV.npy')

qav = AV[0]; uav = AV[1]; vav = AV[2]; hav = AV[3]

# Subtract time-mean
for ti in range(0,Nt):
	f[:,:,ti,:] -= np.array(hav)


# Restrict to upper layer only.
f = f[:,:,:,0]

#==================================================================

# Covariance - this works for data with zero time-mean.

var = np.zeros((N,N))
cov = np.zeros((N,N))

cov = covloop(cov,f,x0,y0,N)
var += np.sum(f**2,axis=2)

#==================================================================

# Now loop through remaining files.

for fi in range(fs+1,fe):

	print(fi)
	# Calculate f again.
	f = np.array(MOMnc.readField(path+files_ss[fi],'LSh'))

	N, Nx, Nt, Nl = np.shape(f)
	Ntot += Nt

	# Subtract time-mean
	for ti in range(0,Nt):
		f[:,:,ti,:] -= np.array(hav)

	# Upper layer only.
	f = f[:,:,:,0]	

	# Update covariance.
	cov = covloop(cov,f,x0,y0,N)
	var += np.sum(f**2,axis=2)

# End loop.
#==================================================================

print(Ntot)

cov /= Ntot
var /= Ntot

# Find horizontal correlation - divide by standard deviations.
corr = cov / np.sqrt(var[y0,x0] * var)

name = 'hCorr_' + str(x0) + '_' + str(y0)
namev = 'hVar_' + str(x0) + '_' + str(y0)

np.save(name,corr)
np.save(namev,var)

#hCorr[y0,x0]=1

e,w,n,s = diagnostics.corrLength(corr,0.3,x0,y0,dx,dy)

print((e+w)/2)
print((n+s)/2)

title = 'h corr (' + str(x0) + ', ' + str(y0) + ')'

plt.pcolor(corr,vmin=-1,vmax=1)
plt.xlim(0,512)
plt.ylim(0,512)
plt.grid(); plt.colorbar(); plt.title(title); plt.show()

plt.pcolor(corr,vmin=-.5,vmax=.5)
plt.xlim(0,512)
plt.ylim(0,512)
plt.grid(); plt.colorbar(); plt.title(title); plt.show()




















