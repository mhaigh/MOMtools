# fCorrelation.py

#==================================================================

# Calculatute correlation between two fields. 
# This code is built to find corr between eddy forcing on a large-scale field F,
# but can be edited so that correlations are between any fields.

#==================================================================

import os 

import numpy as np
import matplotlib.pyplot as plt

import diagnostics.diagnostics as dg
import diagnostics.MOMnc as MOMnc
import diagnostics.plotting as plot

#==================================================================

# How does it work?

# Start with first input file, consider SSH for this example, which we want to correlate with SSH eddy forcing.
# Read large-scale (LS) H, add time-mean component of small-scale (SS) field. Now subtract its full time-mean. (MAY WANT TO RECALCULATE THESE TO MAKE SURE DATA IS CORRECT.)
# Read SS fields (u, v, h) and subtract their time-mean, stored in AV.py.
# Now we have data from this nc file which all has zero time-mean (evaluated over all nc files).
# Calculate flux terms uh, Uh, uH and same for v. Subtract their time-mean, also already stored.
# Add flux terms and take divergence of flux vector.
# Sum over all time-steps to calculate the covariance (F with H) and variances (of F and H).

# Repeat process for all other nc files, adding covariance and variance sums already calculated.
# When done, correlation = cov / sqrt( var(F) var(H) ). 

#==================================================================

str1 = 'N513'
str2 = 'relax'
str3 = 'pre-mean_240'

path_xy = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/' + str1 + '/' + str2 + '/'
path = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/' + str1 + '/' + str2 + '/' + str3 + '/'
path_mean = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/' + str1 + '/' + str2 + '/mean/240km/'
#path_mean = '/home/mike/Documents/GulfStream/MOM/Code/mean/240km/'
path_flux = path_mean + 'fluxes/h/'

# Define nc files and sort.
files_ss = [f for f in os.listdir(path) if f.startswith('ss__')]
files_ls = [f for f in os.listdir(path) if f.startswith('ls__')]
files_ss.sort()
files_ls.sort()

nf = len(files_ss)
fs = 0
fe = nf

# Domain
x, y = MOMnc.readVars(path_xy+'prog__0050_359.nc')

dx = x[1] - x[0]
dy = y[1] - y[0]

lag = 0 # Equivalent to 5 days per index.


#==================================================================

# Read the first field.

# Only read one field.
q, u, v, h = MOMnc.readFieldsLS(path+files_ss[fs])	
Q, U, V, H = MOMnc.readFieldsLS(path+files_ls[fs])

#plot.twoFields([H[:,:,0,0], U[:,:,0,0]], ['H', 'U'])

N, Nx, Nt, Nl = np.shape(h)
Ntot = Nt

# Number of time samples with lag accounted for
Ntl = Nt - lag

# Time index limits accounting for lag.
if lag>=0:
	sl = 0; sr = Ntl
	ll = lag; lr = Nt
else:
	sl = lag; sr = Nt
	ll = 0; lr = Ntl

print(sl,sr)
print(ll,lr)


print('No. of time samples = ' + str(Ntl))

MEAN = 'FILE'
# Get all time-mean data. If looping through all nc files, set MEAN = FILE
if MEAN == 'FILE':

	# Read time-mean
	AV = np.load(path+'AV.npy')
	qav = np.array(AV[0]); uav = np.array(AV[1]); vav = np.array(AV[2]); hav = np.array(AV[3])

	# Small-scale fields have zero time-mean, large-scale fields have nonzero mean.
	Fav = np.load(path_mean+'Hmean.npy')

	# Also need time-mean fluxes; load them here.
	ufav = np.load(path_flux+'uh.npy'); Ufav = np.load(path_flux+'Uh.npy'); uFav = np.load(path_flux+'uH.npy')
	vfav = np.load(path_flux+'vh.npy'); Vfav = np.load(path_flux+'Vh.npy'); vFav = np.load(path_flux+'vH.npy')

# If not reading all files, take mean over one nc file. Make sure to take means of fluxes and LS fields appropriately.
else:
	print('Calculating time-mean small-scale fields...')
	#qav = np.mean(q,axis=2); uav = np.mean(u,axis=2); vav = np.mean(v,axis=2); hav = np.mean(h,axis=2)

	# Read time-mean
	AV = np.load(path+'AV.npy')
	qav = np.array(AV[0]); uav = np.array(AV[1]); vav = np.array(AV[2]); hav = np.array(AV[3])

#==================================================================

# Subtract time-mean from each SS field...
q = dg.subtractAv(q,qav,Nt); u = dg.subtractAv(u,uav,Nt); v = dg.subtractAv(v,vav,Nt); h = dg.subtractAv(h,hav,Nt)
# and add time-mean to LS field.
Q = dg.addAv(Q,qav,Nt); U = dg.addAv(U,uav,Nt); V = dg.addAv(V,vav,Nt); H = dg.addAv(H,hav,Nt);
# DON'T ACCOUNT FOR LAG HERE. THESE LINES JUST MAKE CORRECTIONS TO THE SS AND LS FIELDS.

# Change SS and LS fields here.
f = h.copy()
F = H.copy()

plt.pcolor(h[:,:,0,0]); plt.colorbar(); plt.show()

plot.oneField(h[:,:,0,0],'H')


# Fluxes
uf = u*f; Uf = U*f; uF = u*F
vf = v*f; Vf = V*f; vF = v*F

# If MEAN == 'FILE', already have correct time-mean fluxes and LS fields. Otherwise, calculate them here 
if MEAN != 'FILE':
	print('Calculating time-mean fluxes...')
	ufav = np.mean(uf[:,:,sl:sr,:],axis=2); Ufav = np.mean(Uf[:,:,sl:sr,:],axis=2); uFav = np.mean(uF[:,:,sl:sr,:],axis=2)
	vfav = np.mean(vf[:,:,sl:sr,:],axis=2); Vfav = np.mean(Vf[:,:,sl:sr,:],axis=2); vFav = np.mean(vF[:,:,sl:sr,:],axis=2)
	Fav = np.mean(F[:,:,ll:lr,:],axis=2)

# Subtract mean (need to do this for calculation of cov and var) and account for time lag.
uf = dg.subtractAv(uf[:,:,sl:sr,:],ufav,Ntl); Uf = dg.subtractAv(Uf[:,:,sl:sr,:],Ufav,Ntl); uF = dg.subtractAv(uF[:,:,sl:sr,:],uFav,Ntl)
vf = dg.subtractAv(vf[:,:,sl:sr,:],vfav,Ntl); Vf = dg.subtractAv(Vf[:,:,sl:sr,:],Vfav,Ntl); vF = dg.subtractAv(vF[:,:,sl:sr,:],vFav,Ntl)
F = dg.subtractAv(F[:,:,ll:lr,:],Fav,Ntl)

# Put all unresolved fluxes into one term
uf = uf + Uf + uF
vf = vf + Vf + vF

# Divergence
ufDiv = dg.divLT(uf,vf,dx,dy)
#F = dg.divLT(uf,0*vf,dx,dy)
ufRV = dg.RVLT(uf,vf,dx,dy)
#uhDiv = uh

#plot.twoFields([ufRV[:,:,0,0], ufDiv[:,:,0,0]], ['RV uv', 'Div uv'])

#plot.oneField(np.sqrt(np.mean(uhDiv[:,:,:,0]**2,axis=2)),'std dev FH')
#plot.oneField(1./np.sqrt(np.mean(uhDiv[:,:,:,0]**2,axis=2)),'1 / std dev FH')

# Covariance and variances
cov = np.sum(ufDiv * F, axis=2)
var1 = np.sum(ufDiv**2, axis=2)
var2 = np.sum(F**2, axis=2)

plot.oneField(var1[:,:,0],'std.dev.')

# Save files if True
if False:
	np.save('uh',uf); np.save('vh',vf)
	np.save('uhDiv',ufDiv)	
	np.save('cov',cov)
	np.save('var_uhDiv',var1)
	np.save('var_H',var2)

#==================================================================

# Now loop through remaining files.

for fi in range(fs+1,fe):

	print(fi)

	# Load next nc file
	q, u, v, h = MOMnc.readFieldsLS(path+files_ss[fs])	
	Q, U, V, H = MOMnc.readFieldsLS(path+files_ls[fs])

	f = h.copy()
	F = H.copy()

	Nt = q.shape[2]	
	
	# Subtract and add mean to/from raw data
	q = dg.subtractAv(q,qav,Nt); u = dg.subtractAv(u,uav,Nt); v = dg.subtractAv(v,vav,Nt); h = dg.subtractAv(h,hav,Nt)
	Q = dg.addAv(Q,qav,Nt); U = dg.addAv(U,uav,Nt); V = dg.addAv(V,vav,Nt); H = dg.addAv(H,hav,Nt)

	# Calculate fluxes
	uf = u*f; Uf = U*f; uF = u*F
	vf = v*f; Vf = V*f; vF = v*F

	# Subtract mean from fluxes and large-scale field
	uf = dg.subtractAv(uf,ufav,Nt); Uf = dg.subtractAv(Uf,Ufav,Nt); uF = dg.subtractAv(uF,uFav,Nt)
	vf = dg.subtractAv(vf,vfav,Nt); Vf = dg.subtractAv(Vf,Vfav,Nt); vF = dg.subtractAv(vF,vFav,Nt)
	F = dg.subtractAv(F,Fav,Nt)

	# Collect unresolved fluxes into one term
	uf = uf + Uf + uF
	vf = vf + Vf + vF

	# Calculate divergence
	ufDiv = dg.divLT(uf,vf,dx,dy)

	# Add to covariance and variance sums
	cov += np.sum(ufDiv * F, axis=2)
	var1 += np.sum(ufDiv**2, axis=2)
	var2 += np.sum(F**2, axis=2)

# End loop.
#==================================================================

#uhDiv = dg.smooth(uhDiv,2) # fix smoothing function.

var1inv = 1. / np.sqrt(var1)
var2inv = 1. / np.sqrt(var2)

lim1 = 200
lim2 = 200

var1inv = np.where(var1inv>lim1, lim1, var1inv)
var2inv = np.where(var2inv>lim2, lim2, var2inv)

# Correlation
corr = cov * var1inv * var2inv

plot.oneField(corr[:,:,0],'Corr(F_Hv, F_Hu)')

# [200:400,0:300,0]

corrfile = 'corr_H_Fh'
np.save(corrfile,corr)





