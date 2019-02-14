# windowFilter

#====================================================

import numpy as np
import matplotlib.pyplot as plt

#====================================================

# filter
def filter(F,fw):

	ny, nx = np.shape(F)
	LS = np.zeros((ny,nx))

	fw_gp = 7
	f2 = (fw_gp - 1) / 2
	print(f2)
	for j in range(0,ny):
		if j < f2:
			lj = 0
			rj = j + f2 + 1
		elif j > ny - f2:
			lj = j - f2
			rj = ny	
		else:
			lj = j - f2
			rj = j + f2 + 1
		nj = rj - lj
		print(nj)
		for i in range(0,nx):
			if i < f2:
				li = 0
				ri = i + f2 + 1
			elif i > ny - f2:
				li = i - f2
				ri = nx	
			else:
				li = i - f2
				ri = i + f2 + 1
			ni = ri - li			
		
			LS[j,i] = np.sum(F[lj:rj,li:ri]) / (ni * nj)
			
	SS = F - LS
	
	return SS, LS


