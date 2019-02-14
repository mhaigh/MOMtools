# diagnostics.py

#===============================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#===============================================================

# corrLength
def corrLength(corr,lim,x0,y0,dx,dy):
	'''Find the correlation lengthscale for a given threshold, c.'''

	east = x0
	west = x0
	north = y0
	south = y0	
	
	while (corr[y0,west] > lim):
		west -= 1

	while (corr[y0,east] > lim):
		east += 1

	while (corr[north,x0] > lim):
		north += 1

	while (corr[south,x0] > lim):
		south -= 1

	# Return distance in each direction.
	return (east-x0)*dx, (x0-west)*dx, (north-y0)*dy, (y0-south)*dy		

#===============================================================

# subtractTimeMean
def subtractAv(field, av, Nt):
	'''Subtract time-mean component from field.'''

	for ti in range(0,Nt):
		field[:,:,ti,] -= av

	return field

#===============================================================

# addAv
def addAv(field, av, Nt):
	'''Add time-mean component onto field.'''

	for ti in range(0,Nt):
		field[:,:,ti,] += av

	return field
	

#================================================================

# animate
def animate(u,xg,yg,Nt,vmin,vmax,cm,string,name):
	'''Animate field u. Calls animatei.'''

	def animatei(i):
		'''A function to loop through, used by animate.'''

		plt.cla()
		s = 'day = ' + str(i*5)
		cax.set_array(u[:,:, i])
		
		plt.text(0.20,0.35,s,fontsize=18)
	
	fig = plt.figure()
	ax = plt.subplot(111)

	cax = ax.pcolormesh(u[:,:,0],cmap=cm)#,vmin=vmin,vmax=vmax)
	ax.text(-0.4,-.4,string,fontsize=18)
	ax.axis([xg.min(), xg.max(), yg.min(), yg.max()])
	plt.grid()
	plt.show()

	anim = animation.FuncAnimation(fig, animatei, interval=100, frames=Nt)
 
	plt.tight_layout()
	plt.draw()
	
	anim.save(name+'.mp4',metadata={'artist':'Guido'},writer='ffmpeg',fps=8,bitrate=500)

	plt.show()

#================================================================

# geoBalUV
def geoBalUV(H,dx,dy,f,g):
	'''Calculate geostrophically balanced velocities U and V from SSH field H.
	H has to be SSH, i.e. sum of thickness in each layer.'''

	U = np.zeros(H.shape)
	V = np.zeros(V.shape)

	# Given H, geo. bal. implies: Ug = - g ( dH/dy ) / f and Vg = g ( dH/dx ) / f.	
	U = - g * ddy(H, dy)
	V = g * ddx(H, dx)

	# Now divide by Coriolis.
	for j in range(0,len(f)):
		U[j,] = U[j,] / f[j]
		V[j,] = V[j,] / f[j]
	
	return U, V


#================================================================

# RV
def RV(u,v,dx,dy):
	"Calculate relative vorticity associated with velocity fields u,v"

	# Assume that u,v have only two dimensions
	u_y = diff(u,0,0,dy)
	v_x = diff(v,1,1,dx)

	return v_x - u_y

#================================================================

# RVLT
def RVLT(u,v,dx,dy):
	"Calculate relative vorticity associated with velocity fields u,v"

	return ddx(v,dx) - ddy(u,dy)

#================================================================

# div
def div(u,v,dx,dy):
	"Calculate divergence associated with velocity fields u,v"

	# Assume that u,v have only two dimensions
	u_x = diff(u,1,0,dx)
	v_y = diff(v,0,0,dy)

	return u_x + v_y

#================================================================

# divLT
def divLT(u,v,dx,dy):
	'''Calculate divergence associated with velocity fields u,v.
	General function for arbitrary number of time-steps and layers.'''

	return ddx(u,dx) + ddy(v,dy)

#================================================================

def flux(uss,uls,vss,vls,fss,fls):

	uf = np.sum(uss * fss, axis=2)
	vf = np.sum(vss * fss, axis=2)

	Uf = np.sum(uls * fss, axis=2)
	Vf = np.sum(vls * fss, axis=2)

	uF = np.sum(uss * fls, axis=2)
	vF = np.sum(vss * fls, axis=2)

	UF = np.sum(uls * fls, axis=2)
	VF = np.sum(vls * fls, axis=2)

	return uf, vf, Uf, Vf, uF, vF, UF, VF

#================================================================

# smooth
def smooth(f,f2):
	'''Smooth a field f using a window filter'''

	shape = np.shape(f)
	Ny = shape[0]
	Nx = shape[1]

	fsmooth = np.zeros((shape))

	for j in range(0,Ny):
		# Set left and right filter end points
		if j < f2:
			lj = 0
			rj = j + f2 + 1
		elif j > Ny - f2 - 1:
			lj = j - f2
			rj = Ny
		else:
			lj = j - f2
			rj = j + f2 + 1
		nj = rj - lj
		
		for i in range(0,Nx):
			# Set left and right filter end points
			if i < f2:
				li = 0
				ri = i + f2 + 1
			elif i > Nx - f2 - 1:
				li = i - f2
				ri = Nx
			else:
				li = i - f2
				ri = i + f2 + 1
			ni = ri - li			
		
			# Filter
			a=np.sum(f[lj:rj,li:ri,],axis=(0,1)) / (ni * nj)
			fsmooth[j,i,] = a			

	return fsmooth


#================================================================

# ddx
def ddx(f,dx):
	'''Derivative w.r.t x of field f.
	Written for div and RV calculations to speed up algorithm
	by taking time and layer dependence outside of loops.'''

	# Initialise derivative
	Ny = np.shape(f)[0]
	Nx = np.shape(f)[1]
	df = np.zeros(f.shape,dtype=f.dtype)

	# Find derivative
	df[:,1:Nx-1,] = f[:,2:Nx,] - f[:,0:Nx-2,]

	# Boundary terms
	df[:,0,] = 2. * (f[:,1,] - f[:,0,])
	df[:,Nx-1,] = 2. * (f[:,Nx-1,] - f[:,Nx-2,])

	# Divide by space step
	df = 0.5 * df / dx

	return df

#================================================================

# ddy
def ddy(f,dy):
	'''Derivative w.r.t x of field f.
	Written for div and RV calculations to speed up algorithm
	by taking time and layer dependence outside of loops.'''

	Ny = np.shape(f)[0]
	Nx = np.shape(f)[1]
	df = np.zeros(f.shape,dtype=f.dtype)

	# Find derivative
	df[1:Ny-1,] = f[2:Ny,] - f[0:Ny-2,]

	# Boundary terms
	df[0,] = 2. * (f[1,] - f[0,])
	df[Ny-1,] = 2. * (f[Ny-1,] - f[Ny-2,]);

	# Divide by space step
	df = 0.5 * df / dy

	return df
	
#================================================================

# diff
def diff(f,d,p,delta):
# Function for differentiating a vector
# f[y,x] is the function to be differentiated.
# d is the direction in which the differentiation is to be taken:
# d=0 for differentiation over the first index, d=1 for the second.
# d=2 for a 1D vector
# p is a periodic switch:
# p=1 calculates periodic derivative.
# Need to be careful with periodic derivatives, output depends on whether f[0]=f[dim-1] or =f[dim].
	
	if d != 2:
		dimx = np.shape(f)[1];		# Finds the number of gridpoints in the x and y directions
		dimy = np.shape(f)[0];
		df = np.zeros((dimy,dimx),dtype=f.dtype);
	else:
		dimy = np.shape(f)[0];
		df = np.zeros(dimy,dtype=f.dtype);
	
	if p == 0:
		# Solid boundary derivative.
		# Note multiplication by 2 of boundary terms are to
		# invert the division by 2 at the end of the module.
		if d == 0:
		
			df[1:dimy-1,:] = f[2:dimy,:] - f[0:dimy-2,:];
		
			df[0,:] = 2 * (f[1,:] - f[0,:]);
			df[dimy-1,:] = 2 * (f[dimy-1,:] - f[dimy-2,:]);
	
		elif d == 1:

			df[:,1:dimx-1] = f[:,2:dimx] - f[:,0:dimx-2];

			df[:,0] = 2 * (f[:,1] - f[:,0]);
			df[:,dimx-1] = 2 * (f[:,dimx-1] - f[:,dimx-2]);
		
		elif d == 2:

			df[1:dimy-1] = f[2:dimy] - f[0:dimy-2];

			df[0] = 2 * (f[1] - f[0]);
			df[dimy-1] = 2 * (f[dimy-1] - f[dimy-2]);

		else:
			print('error')

	elif p == 1:
		# Periodic option

		if d == 0:

			df[1:dimy-1,:] = f[2:dimy,:] - f[0:dimy-2,:];	

			df[0,:] = f[1,:] - f[dimy-1,:];
			df[dimy-1,:] = f[0,:] - f[dimy-2,:];
	
		elif d == 1:

			df[:,1:dimx-1] = f[:,2:dimx]-f[:,0:dimx-2];

			df[:,0] = f[:,1] - f[:,dimx-1];
			df[:,dimx-1] = f[:,0] - f[:,dimx-2];

		elif d == 2:

			df[1:dimy-1] = f[2:dimy] - f[0:dimy-2];

			df[0] = f[1] - f[dimy-2];
			df[dimy-1] = f[1] - f[dimy-2];
			
			#print(str(df[0])+'='+str(f[1])+'+'+str(f[dimy-1]));
			#print(str(df[dimy-1])+'='+str(f[0])+'+'+str(f[dimy-2]));
		
		else:
			print('error')

	else:
		print('error')

	df = 0.5 * df / delta;

	return df
