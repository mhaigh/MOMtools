# readMOMnc.py
#=======================================================

# This files read netcdf data, as produced by MOM6. 

#=======================================================

import numpy as np
import netCDF4 as nc

#=======================================================

def readFields(ncFile):

	f = nc.Dataset(ncFile)	# Read the NETCDF file

	q = np.array(f.variables['PV'])
	u = np.array(f.variables['u'])
	v = np.array(f.variables['v'])
	h = np.array(f.variables['h'])

	return q, u, v, h

#=======================================================

def readFieldsLS(ncFile):
	'''Small-scale and large-scale components both use LS... identifier.
	Just make sure to read from correct directory.'''

	f = nc.Dataset(ncFile)	# Read the NETCDF file

	# Confusingly, small-scale components use the identifier LS...
	q = np.array(f.variables['LSq'])
	u = np.array(f.variables['LSu'])
	v = np.array(f.variables['LSv'])
	h = np.array(f.variables['LSh'])

	return q, u, v, h

#=======================================================

def readField(ncFile,field):

	# field should be string identifier for field in ncFile.

	f = nc.Dataset(ncFile)	# Read the NETCDF file

	# Confusingly, small-scale components use the identifier LS...
	return np.array(f.variables[field])
	

#=======================================================

def readVars(ncFile):

	f = nc.Dataset(ncFile)	# Read the NETCDF file

	x = f.variables[str('xq')]
	y = f.variables[str('yq')]

	return x, y

#=======================================================

def writeFields(path,prog,LS,fname):
	''' Save decomposed fields, ls or ss'''

	filename = path + fname + prog[4::];

	print('Writing to nc file' + filename)

	sz = np.shape(LS[0])

	# Initialise file
	ls_file = nc.Dataset(filename,'w',format='NETCDF4')

	# Dimensions
	yd = ls_file.createDimension('yd',sz[0])
	xd = ls_file.createDimension('xd',sz[1])
	td = ls_file.createDimension('td',sz[2])
	ld = ls_file.createDimension('ld',sz[3])

	# Create dependent variables
	#y = RSW1L_Eigenmodes.createVariable('y','f8',('yd',)); 
	#x = RSW1L_Eigenmodes.createVariable('x','f8',('xd',)); 
	#t = RSW1L_Eigenmodes.createVariable('t','f8',('td',)); 
	

	# Assign values (dx=1 for now)
	#y[:] = np.linspace(0,sz[0]-1,sz[0])
	#x[:] = np.linspace(0,sz[1]-1,sz[1])
	#t[:] = np.linspace(0,sz[2]-1,sz[2])

	# Initialise fields
	lsq = ls_file.createVariable('LSq','f8',('yd','xd','td','ld',))
	lsu = ls_file.createVariable('LSu','f8',('yd','xd','td','ld',))
	lsv = ls_file.createVariable('LSv','f8',('yd','xd','td','ld',))
	lsh = ls_file.createVariable('LSh','f8',('yd','xd','td','ld',))

	# Assign values
	lsq[:,:,:,:] = LS[0]; lsu[:,:,:,:] = LS[1]
	lsv[:,:,:,:] = LS[2]; lsh[:,:,:] = LS[3]
	
	ls_file.close()

	return

#=======================================================

def writeFieldsSS(path,file_id,SS):
# Save large-scale PV
	
	file_name = path + 'ss' + file_id;

	sz = np.shape(SS[0])

	# Initialise file
	ss_file = nc.Dataset(file_name,'w',format='NETCDF4')

	# Dimensions
	yd = ss_file.createDimension('yd',sz[0])
	xd = ss_file.createDimension('xd',sz[1])
	td = ss_file.createDimension('td',sz[2])
	ld = ss_file.createDimension('ld',sz[3])

	# Initialise fields
	ssq = ss_file.createVariable('SSq','f8',('yd','xd','td','ld',))
	ssu = ss_file.createVariable('SSu','f8',('yd','xd','td','ld',))
	ssv = ss_file.createVariable('SSv','f8',('yd','xd','td','ld',))
	ssh = ss_file.createVariable('SSh','f8',('yd','xd','td','ld',))

	# Assign values
	ssq[:,:,:,:] = SS[0]; ssu[:,:,:,:] = SS[1]
	ssv[:,:,:,:] = SS[2]; ssh[:,:,:] = SS[3]
	
	ss_file.close()

	return

#=======================================================

def readFieldAll(ncFile):

	f = nc.Dataset(ncFile)	# Read the NETCDF file

	# Extract the eigenmodes and reconstruct from real and imag. parts.
	LS = [f.variables['LSq'],f.variables['LSu'],f.variables['LSv'],f.variables['LSh']]
	
	return LS

#=======================================================

def a():
	# Initialise the nc file
	RSW1L_Eigenmodes = nc.Dataset(file_name,'w',format='NETCDF4');
		
	# Create dimensions
	y3_dim = RSW1L_Eigenmodes.createDimension('y3_dim',dim);
	omega_dim = RSW1L_Eigenmodes.createDimension('omega_dim',dim);
	real_imag = RSW1L_Eigenmodes.createDimension('real_imag',2);

	# Initialise dimension variables...
	y3 = RSW1L_Eigenmodes.createVariable('y3','f8',('y3_dim',)); # Here y3_dim is 3 times the domain, minus some endpoints
	omega = RSW1L_Eigenmodes.createVariable('omega','f8',('omega_dim','real_imag',));
	# ...and assign the data.
	if BC == 'NO-SLIP':		
		y3[0:N-2] = y_nd[1:N-1];
		y3[N-2:2*N-4] = y_nd[1:N-1];
		y3[2*N-4:3*N-4] = y_nd;
	if BC == 'FREE-SLIP':
		y3[0:N] = y_nd[0:N];
		y3[N:2*N-2] = y_nd[1:N-1];
		y3[2*N-2:3*N-2] = y_nd;		
	omega[:,0] = np.real(val); omega[:,1] = np.imag(val);
	
	# Initialise solution variables...
	vec = RSW1L_Eigenmodes.createVariable('vec','f8',('y3_dim','omega_dim','real_imag',));
	count = RSW1L_Eigenmodes.createVariable('count','i4',('y3_dim',));
	
	# ...and assign the data to the variables.
	vec[:,:,0] = np.real(modes); vec[:,:,1] = np.imag(modes);
	count[:] = zero_count;

	RSW1L_Eigenmodes.close();
	

	
