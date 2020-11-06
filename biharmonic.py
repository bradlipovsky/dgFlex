#
# Solves the beam flexure equations in a way meant to reflect the drainage of supraglacial lakes on ice shelves.
#
# Uses a DG approach. 
#
# based off this fenics demo:  https://fenicsproject.org/docs/dolfin/latest/python/demos/biharmonic/demo_biharmonic.py.html
#
# I set up the environment using, conda create -n fenicsproject -c conda-forge mshr=2019.1.0=py38h2af9582_2 scipy netcdf4 fenics matplotlib 
#

import math
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import milne_boundaries as bdd
from dolfin import *
from scipy.interpolate import interp2d
import netCDF4 as nc


rhow = 1000
rho = 916
g = 9.8
hw = 1
nu=0.3
E = 8.7e9

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

# Make mesh ghosted for evaluation of DG terms
parameters["ghost_mode"] = "shared_facet"

# Create mesh and define function space
mesh = Mesh('milne_v22_500-5m.xml')
V = FunctionSpace(mesh, "CG", 2)

def get_thickness(coords,thickinterp,r):
	xCoord = coords[:,0]
	yCoord = coords[:,1]
	n = len(xCoord)
	hCoord = np.empty_like(xCoord)
	for i in range(0,n):
		# Note correction for surface height -> ice thickness
		hCoord[i] = thickinterp(xCoord[i],yCoord[i]) *r
		if hCoord[i] < 1:
			hCoord[i] = 1
	return hCoord

ds = nc.Dataset('milne_thickness_ArcticDEM_32m.nc')
thickness_interpolator= interp2d( ds['X'][:],ds['Y'][:], ds['H'][:])
h = Function(V)
MeshCoordinates=V.tabulate_dof_coordinates()
h.vector()[:] = get_thickness(MeshCoordinates,thickness_interpolator,
			1/(1-rho/rhow))

EI = E/(1-nu**2)/12 * h**3
h_ref = 100
EI_ref = E/(1-nu**2)/12 * h_ref**3
#lambd = (EI/rho/g)**(1/4)
#print (lambd)

class Source(UserExpression):
    def eval(self, values, x):
	# Plus sign because the lake drained
        values[0] = +rhow*g*hw * FractureZone(x)


def FractureZone(x):
	return np.sqrt(  (x[1] + 6.35e5)**2 + (x[0] + 4.72e5)**2  )  < 1e3


class DirichletBoundary(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary		

# Define boundary condition
u0 = Constant(0.0)
#bc = DirichletBC(V, u0, DirichletBoundary() )
bc = DirichletBC(V, u0, bdd.clamped_margins())

# Define trial and test functions
u = TrialFunction(V)
v = TestFunction(V)

# Define normal component, mesh size and right-hand side
h = CellDiameter(mesh)
h_avg = (h('+') + h('-'))/2.0
n = FacetNormal(mesh)
f = Source(degree=2)

# Penalty parameter
alpha = Constant(EI_ref*8.0)

# Define bilinear form
a = inner(EI*div(grad(u)), div(grad(v)))*dx \
  - inner(avg(EI*div(grad(u))), jump(grad(v), n))*dS \
  - inner(jump(grad(u), n), avg(EI*div(grad(v))))*dS \
  + alpha/h_avg*inner(jump(grad(u),n), jump(grad(v),n))*dS\
  + rho*g*u*v*dx

# Define linear form
L = f*v*dx

# Solve variational problem
u = Function(V)
solve(a == L, u, bc)

fig,ax=plt.subplots(figsize=(5,5))
cw = cm.get_cmap('coolwarm',512)
vm = 1
c=plot(u*10,mode='color',cmap=cw,vmin=-vm, vmax=vm)
plt.colorbar(c)
plt.show()
