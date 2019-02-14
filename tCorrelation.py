# sCorrelation.py

#==================================================================

# Calculatute lagged autocorrelation and horiztonal correlation in field.

#==================================================================

import os
import sys 

import numpy as np
import matplotlib.pyplot as plt

import diagnostics.diagnostics as dg
import diagnostics.MOMnc as MOMnc
import diagnostics.plotting as plot

#==================================================================

# How does it work?

# Spatial correlation first.
# Select a grid point x0, at which we want its spatial correlation statistics.
# We need the variance (st. dev.) of the field at ALL grid points x, which we may already have.
# Take the covariance of the field evaluated at x0 with same field evaluated at x, divide by their standard deviations.

# If finding corr for fluxes, we can take the mean in saved numpy files.

#==================================================================

def corr1():

	corr = np.zeros((N,Nx,Nt))

	# Find first lag, i.e. variance.
	corr[:,:,0,] = np.sum(f**2,axis=2)

	# All other lags, normalise by variance
	for ti in range(1,Nt):

		# Correlation.
		corr[:,:,ti,] = np.sum(f[:,:,0:Nt-ti,] * f[:,:,ti:Nt,], axis=2) / corr[:,:,0,]

	corr[:,:,0,] = 1

	return corr

#==

def corrWSS():

	corr = np.zeros((N,Nx,Nt))

	# Find first lag, i.e. variance.
	corr[:,:,0,] = np.mean(f[:,:,0:Nt]**2,axis=2)

	# All other lags, normalise by variance
	for ti in range(1,Nt):

		# Calculate Lagged data, subtract mean, find variance.
		fLag = f[:,:,ti:Nt+ti,]

		# Correlation.
		corr[:,:,ti,] = np.mean(f[:,:,0:Nt,] * fLag, axis=2) / corr[:,:,0,]

	corr[:,:,0,] = 1

	return corr

#==================================================================

str1 = 'N513'
str2 = 'relax'
str3 = 'pre-mean_240'

path_xy = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/' + str1 + '/' + str2 + '/'
path = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/' + str1 + '/' + str2 + '/' + str3 + '/'

files_ss = [f for f in os.listdir(path) if f.startswith('ss__')]
files_ls = [f for f in os.listdir(path) if f.startswith('ls__')]
files_ss.sort()
files_ls.sort()

nf = len(files_ss)
fs = 15

x, y = MOMnc.readVars(path_xy+'prog__0050_359.nc')
AV = np.load(path+'AV.npy')

dx = x[1] - x[0]
dy = y[1] - y[0]

layer = 0

#==================================================================

print("Reading data...")

# Read from two nc files.
u = MOMnc.readField(path+files_ss[fs],'LSu')[:,:,:,layer]
v = MOMnc.readField(path+files_ss[fs],'LSv')[:,:,:,layer]
h = MOMnc.readField(path+files_ss[fs],'LSh')[:,:,:,layer]

U = MOMnc.readField(path+files_ls[fs],'LSu')[:,:,:,layer]
V = MOMnc.readField(path+files_ls[fs],'LSv')[:,:,:,layer]
H = MOMnc.readField(path+files_ls[fs],'LSh')[:,:,:,layer]

Nt = np.shape(u)[2]

if True:

	utmp = MOMnc.readField(path+files_ss[fs+1],'LSu')[:,:,:,layer]
	vtmp = MOMnc.readField(path+files_ss[fs+1],'LSv')[:,:,:,layer]
	htmp = MOMnc.readField(path+files_ss[fs+1],'LSh')[:,:,:,layer]

	Utmp = MOMnc.readField(path+files_ls[fs+1],'LSu')[:,:,:,layer]
	Vtmp = MOMnc.readField(path+files_ls[fs+1],'LSv')[:,:,:,layer]
	Htmp = MOMnc.readField(path+files_ls[fs+1],'LSh')[:,:,:,layer]

	# Concatentate the data.
	u = np.concatenate((u,utmp),axis=2); v = np.concatenate((v,vtmp),axis=2); h = np.concatenate((h,htmp),axis=2)
	U = np.concatenate((U,Utmp),axis=2); V = np.concatenate((V,Vtmp),axis=2); H = np.concatenate((H,Htmp),axis=2)

# Dimensions.
N, Nx, Nt2 = np.shape(u)

print('No. of time samples = ' +str(np.shape(u)[2]))

print("Averaging data...")
# Load time-mean data to make corrections to LS and SS field. Another option is to ensure SS field has zero time-mean over the active nc files.
AV = np.load(path+'AV.npy')[:,:,:,layer]

# Subtract time-mean from each SS field...
u = dg.subtractAv(u,AV[1],Nt); v = dg.subtractAv(v,AV[2],Nt); h = dg.subtractAv(h,AV[3],Nt)
# and add time-mean to LS field.
U = dg.addAv(U,AV[1],Nt); V = dg.addAv(V,AV[2],Nt); H = dg.addAv(H,AV[3],Nt)


print("Preparing field...")
# What data do we want time-correlations of? Store it in f.
if True:
	uh = u*h + U*h + u*H
	vh = v*h + V*h + v*H

	f = dg.divLT(uh*0,vh,dx,dy)

if False:
	f = H.copy()

# Require average
fav = np.mean(f,axis=2)
f = dg.subtractAv(f,fav,np.shape(f)[2])

# Define lags
lags = [5.*ti for ti in range(Nt)]

#==================================================================

print('Calculating correlations...')

corr = corrWSS()
#corr = corr1()

np.save('tCorr_FH',corr)
	
#==================================================================

# Plotting

xp = [N-30,30,20,20]
yp = [300,300,400,200]

for i in range(4):
	plt.plot(lags,corr[yp[i],xp[i],:],label='('+str(xp[i])+','+str(yp[i])+')')

plt.legend()
plt.grid()
plt.xlabel('Lag (days)')
plt.ylabel('Autocorrelation')

plt.show()


sys.exit()

# Code below here is for partial autocorrelation function.




