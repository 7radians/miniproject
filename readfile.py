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
xax = dat.variables['x_axis'][:]
yax = dat.variables['y_axis'][:]


# set [-1. 1] axes for the plots and make them look neat
xaxis = xax[0::10]
xaxis[-1] = round(xaxis[-1], 1)
xaxis[0] = round(xaxis[0], 1)

yaxis = yax[0::10]
yaxis[-1] = round(yaxis[-1], 1)
yaxis[0] = round(yaxis[0], 1)

xnum = np.linspace(0, len(xax), len(xaxis))
ynum = np.linspace(0, len(yax), len(yaxis))


#a pseudocolour plot
plt.plot()
plt.xlabel('x position')
plt.ylabel('y position')
plt.imshow(ex,  interpolation='none')
plt.title('Electric field Ex')
plt.colorbar(cmap='jet', label ='Electric field')
plt.xticks(xnum, xaxis)
plt.yticks(ynum, yaxis)
plt.show()


#a scatter plot
plt.plot()
plt.scatter(xpos, ypos)
plt.xlabel("x position")
plt.ylabel("y position")
plt.title("Particle position")
plt.xticks(xaxis)
plt.yticks(yaxis)
plt.show()


