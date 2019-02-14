# filters

#====================================================

import numpy as np
import matplotlib.pyplot as plt

#====================================================

# spectralFilter
def spectralFilter(F,fw):

	L = 3840000.

	ny, nx = np.shape(F)

	# Large-scale, small-scale outputs.
	LS = np.zeros((ny,nx))
	SS = np.zeros((ny,nx))

	# Define a doubly periodic field
	dp_field = np.zeros((2*ny,2*nx))

	# Bottom-left corner
	dp_field[0:ny,0:nx] = F
	# Top-left corner
	dp_field[ny:2*ny,0:nx] = F[::-1,:]
	# Right-hand side
	dp_field[:,nx:2*nx] = dp_field[:,0:nx][:,::-1]

	# Now pick the wavenumbers in i,j directions
	dx = L / (nx)
	
	# May need to fix this as we have removed final grid point
	k = np.fft.fftfreq(2*nx,dx)
	l = np.fft.fftfreq(2*ny,dx)
	
	# Small-scale
	kl_ss = [(j,i) for i in range(0,2*nx) for j in range(0,2*ny) if ((k[i]**2 + l[j]**2)**0.5 < 1./fw)]
	#print(len(kl_ss)
	ft_ss = np.fft.fft2(dp_field)
	ft_ss[kl_ss] = 0.
	SS = np.real(np.fft.ifft2(ft_ss)[0:ny,0:nx])
	
	LS = F - SS
	
	return SS, LS

#====================================================

# windowFilter
def windowFilter(F,fw):

	shape = np.shape(F)
	ny = shape[0]
	nx = shape[1]	

	# Intialise large-scale decomposition
	LS = np.zeros(shape)

	dx = 3840000. / ny
	# Number of grid points in each direction in the filter.
	fwgp = int(fw / dx) + 1
	#print(fwgp)
	
	# dx = 7.5 km for high-res (512) simulations
	# fwgp = n --> (n-1) * 7.5 km
	# fwgp = 5 --> 30 km ~ 1 * Rd
	# fwgp = 7 --> 60 km ~ 2 * Rd
	# fwgp = 9 --> 90 km ~ 3 * Rd

	# Number of grid points either side of (i,j)
	f2 = (fwgp - 1) / 2

	
	for j in range(0,ny):
		# Set left and right filter end points
		if j < f2:
			lj = 0
			rj = j + f2 + 1
		elif j > ny - f2 - 1:
			lj = j - f2
			rj = ny
		else:
			lj = j - f2
			rj = j + f2 + 1
		nj = rj - lj
		
		for i in range(0,nx):
			# Set left and right filter end points
			if i < f2:
				li = 0
				ri = i + f2 + 1
			elif i > nx - f2 - 1:
				li = i - f2
				ri = nx
			else:
				li = i - f2
				ri = i + f2 + 1
			ni = ri - li	

			# Filter
			LS[j,i,] = np.sum(F[lj:rj,li:ri,],axis=(0,1)) / (ni * nj)
			
	SS = F - LS
	
	return SS, LS

