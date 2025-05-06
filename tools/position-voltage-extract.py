# author Dmitrii Maksimov
# modifications to python3 compatibility only: Mariana Rossi 2025
from sys import argv
import numpy as np

bohrtoang=0.52917721
e0=8.854187817e-12 #SI units

cube = open(argv[1])

trash=cube.readline()
trash=cube.readline()

line=cube.readline()
natoms=int(line.split()[0])
x0, y0, z0 =list(map(float, line.split()[1:])) #these are in bohr

line=cube.readline()
nx=int(line.split()[0])
dx=np.array(list(map(float, line.split()[1:])))*bohrtoang
line=cube.readline()
ny=int(line.split()[0])
dy=np.array(list(map(float, line.split()[1:])))*bohrtoang
line=cube.readline()
nz=int(line.split()[0])
dz=float(line.split()[3])*bohrtoang
zdens=np.zeros(nz)

for i in range(natoms):
    trash=cube.readline()

nlines=int(nz/6.0)+1 #takes care of the convention of cube files to output 6 fields per line

X0 = np.linspace(0, dx*nx, nx)[:,0]
Y0 = np.linspace(0, dy*ny, ny)[:,1]
Z0 = np.linspace(0, dz*nz, nz)
# following loop assumes the convention of cube files that z direction is the innermost loop.

All_values = np.array([])

for ix in range(nx):
   for iy in range(ny):
      whole = [list(map(float, cube.readline().split())) for nl in range(nlines)]
      whole=np.hstack(whole)
      #zdens+=whole
      All_values = np.append(All_values, np.array(whole))

with open('XYZ_cube.dat', 'w') as xyzV:
    xyzV.write('{:<20}{:<20}{:<20}{:<20}\n'.format('X', 'Y', 'Z', 'V'))
    k=-1      
    for ix in X0:
        for iy in Y0:
            for iz in Z0:
                k+=1
                xyzV.write('{:<20}{:<20}{:<20}{:<20}\n'.format(ix, iy, iz, All_values[k]))
            
      
      
