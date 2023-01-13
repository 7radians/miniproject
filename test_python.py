import netCDF4 as NC
import numpy as np
# Read netcdf file
dat = NC.Dataset('output_data', 'r', format='NETCDF4')

print(str(dat.data_writer))
print(str(dat.nx))
print(str(dat.ny))
print(str(dat.problem))

var = dat.variables['x_axis']
print(np.array(var[:]))
var = dat.variables['y_axis']
print(np.array(var[:]))
var = dat.variables['t_axis']
print(np.array(var[:]))

var = dat.variables['rho']
print(np.array(var[:]))
var = dat.variables['phi']
print(np.array(var[:]))
var = dat.variables['Ex']
print(np.array(var[:]))
var = dat.variables['Ey']
print(np.array(var[:]))

var = dat.variables['x']
print(np.array(var[:]))
var = dat.variables['y']
print(np.array(var[:]))
var = dat.variables['vx']
print(np.array(var[:]))
var = dat.variables['vy']
print(np.array(var[:]))
var = dat.variables['ax']
print(np.array(var[:]))
var = dat.variables['ay']
print(np.array(var[:]))