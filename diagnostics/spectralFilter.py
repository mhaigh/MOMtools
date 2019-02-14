# spectralFilter

#====================================================

import numpy as np
import matplotlib.pyplot as plt

#====================================================

# filter
def filter(F,fw):

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
	#print(len(kl_ss))
	ft_ss = np.fft.fft2(dp_field)
	ft_ss[kl_ss] = 0.
	SS = np.real(np.fft.ifft2(ft_ss)[0:ny,0:nx])
	
	LS = F - SS
	
	return SS, LS


