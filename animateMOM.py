# decompMOM.py
#=======================================================

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import timeFluxes
import spectralFilter
import MOMnc

#=======================================================

# Parameters

panels = 1

str1 = 'N513'
str2 = 'relax'
str3 = 'pre-mean_240'

path_orig = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/' + str1 + '/' + str2 + '/'
path = '/media/mike/Seagate Expansion Drive/Documents/MOM_data/benchmark/' + str1 + '/' + str2 + '/' + str3 + '/'

#============================================================


# Read all files
x, y = MOMnc.readVars(path_orig+'prog__0050_359.nc')

q = timeFluxes.fluxes(path)

#=======================================================
#=======================================================

#LS = MOMnc.readField(path + 'ls' + files[fi][4::])
#SS = q[:,:,:,0] - LS[:,:,:]

if panels == 1:

	shape = np.shape(q)
	print(shape)
	print(shape[0])
	N = shape[0]; Nx = shape[1]; Nt = shape[2]

	N = N-1
	
	#q = q[0:N,0:N,:,0]
	q = q[0:N,0:N,:]

	Lx = 3840000.
	Ly = 3840000.

	dx = 1. / (N)
	dy = 1. / (N)

	yg, xg = np.mgrid[slice(-1./2,1./2+dy,dy),slice(-1./2,1./2+dx,dx)]

	#str1 = 'conv(uh+Uh+uH)'	
	str1 = 'conv(uh)'

	cm = 'bwr'
	#cm = 'jet'
	#ss_lim = np.max(np.abs(SS))
	vmax = np.max(np.abs(q)) / 20.
	vmin = -vmax

	fig = plt.figure()

	ax1 = plt.subplot(111)

	cax1 = ax1.pcolormesh(xg,yg,q[:,:,0],cmap=cm,vmin=vmin,vmax=vmax)
	ax1.axis([xg.min(), xg.max(), yg.min(), yg.max()])
	ax1.text(-.4,-.4,str1,fontsize=18)
	fig.colorbar(cax1)
	plt.grid()

	def animate(i):
		#plt.cla()		
		s = 'day = ' + str(i*5)
		cax1.set_array(q[:,:, i].flatten())
		#cax1.set_array(q[:-1,:-1,i].flatten())
		
	anim = animation.FuncAnimation(fig, animate, interval=100, frames=Nt)
 
	plt.tight_layout()
	plt.draw()

	anim.save('PV.mp4',metadata={'artist':'Guido'},writer='ffmpeg',fps=8,bitrate=500)

	plt.show()


elif panels == 2:

	q1 = np.load('q.npy'); N,Nx,Nt = np.shape(q1); q2=q1.copy()
	#N, Nx, Nt, Nl = np.shape(q)

	N = N-1
	
	#q = q[0:N,0:N,:,:]

	Lx = 3840000.
	Ly = 3840000.

	dx = 1. / (N)
	dy = 1. / (N)

	yg, xg = np.mgrid[slice(-1./2,1./2+dy,dy),slice(-1./2,1./2+dx,dx)]
	yg_low, xg_low = np.mgrid[slice(-1./2,1./2+dy*4,4*dy),slice(-1./2,1./2+4*dx,4*dx)]


	str1 = 'Uq'
	str2 = 'Vq'

	cm = 'bwr'
	#cm = 'jet'
	#ss_lim = np.max(np.abs(SS))
	vmax = np.max(np.abs(q1)) / 2.
	vmin = -vmax

	fig = plt.figure(figsize=[13,5])

	ax1 = fig.add_subplot(121)
	ax2 = fig.add_subplot(122)

	#cax1 = ax1.pcolormesh(xg_low,yg_low,q_low[:,:,0],cmap=cm,vmin=vmin,vmax=vmax)
	#ax1.axis([xg_low.min(), xg_low.max(), yg_low.min(), yg_low.max()])
	#fig.colorbar(cax1)
	cax1 = ax1.pcolormesh(xg,yg,q1[:,:,0],cmap=cm,vmin=vmin,vmax=vmax)
	ax1.text(-0.4,-.4,str1,fontsize=18)
	ax1.axis([xg.min(), xg.max(), yg.min(), yg.max()])
	plt.grid()

	cax2 = ax2.pcolormesh(xg,yg,q2[:,:,0],cmap=cm,vmin=vmin,vmax=vmax)
	ax2.axis([xg.min(), xg.max(), yg.min(), yg.max()])
	fig.colorbar(cax2)
	plt.grid()

	def animate(i):
		plt.cla()
		cax2 = ax2.pcolormesh(xg,yg,q2[:,:,0],cmap=cm,vmin=vmin,vmax=vmax)
		ax2.axis([xg.min(), xg.max(), yg.min(), yg.max()])
		s = 'day = ' + str(i*5)
		cax1.set_array(q1[:-1,:-1, i].flatten())
		
		#cax2.set_array(q2[:-1,:-1,i].flatten())
		cax2.set_array(q2[:-1,:-1,i].flatten())
		plt.text(0.20,0.35,s,fontsize=18)
		plt.text(-0.4,-0.4,str2,fontsize=18)
		#cax2.set_array(q_low[:-1,:-1, i].flatten())
	anim = animation.FuncAnimation(fig, animate, interval=100, frames=Nt)
 
	plt.tight_layout()
	plt.draw()

	anim.save('PV.mp4',metadata={'artist':'Guido'},writer='ffmpeg',fps=8,bitrate=500)

	plt.show()









