import netCDF4
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os

# read the variables
path = os.getcwd()
path2 = path + '/output.txt'
print(path2)
file = netCDF4.Dataset(path2, 'r')
#edit variable names
varname1 = file.variables['varname1'][:, :]
varname2 = file.variables['varname2'][:]
varname3 = file.variables['varname3'][:]

#a pseudocolour plot
plt.plot()
#plt.figure()
plt.xlabel("x")
plt.ylabel("y")
plt.imshow(varname1,  interpolation='none')
plt.title("Plot title")
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax)
plt.show()
    
#a scatter plot
plt.scatter(varname2, varname3)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Plot title")
plt.show()

