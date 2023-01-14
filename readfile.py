import netCDF4
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as NC

#get data
dat = NC.Dataset('output_data', 'r', format='NETCDF4')

#extract variables

ex = dat.variables['Ex'][:,:]
xpos = dat.variables['x'][:]
ypos = dat.variables['y'][:]

#a pseudocolour plot
plt.plot()
plt.xlabel("x")
plt.ylabel("y")
plt.imshow(ex,  interpolation='none')
plt.title("Electric field Ex")
plt.colorbar()
plt.show()
    
#a scatter plot
plt.plot()
plt.scatter(xpos, ypos)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Particle position")
plt.show()


