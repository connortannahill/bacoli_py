# Solving the Reaction-Convection-Diffusion (RCD) system. This
# problem is also called the Catalytic Surface Reaction Model 
# (CSRM). It consists of four PDE's.
# ------------------------------------------------------------------
# This example is based off a FORTRAN analogue for the original 
# BACOLI. The original can be found at:
#   http://cs.stmarys.ca/~muir/BACOLI95-3_Source/3-Problems/RCDsys.f
# ------------------------------------------------------------------

import bacoli_py
import numpy

solver = bacoli_py.BacoliPy()
npde = 4

# Initialize problem-dependent parameters.
a1 = 30.0
a2 = 30.0
d1 = 1.50
d2 = 1.20
r = 1000.0
c = .96
n = 1.0
pe1 = 1.0 * 10**4
pe2 = 1.0 * 10**4

# Set t0.
initial_time = 0.0

# Function defining the PDE system.
def f(t, x, u, ux, uxx, fval):
    fval[0] = -ux[0]+n*(d1*u[2]-a1*u[0]*(1.0-u[2]-u[3])) + (1.0/pe1)*uxx[0]
    fval[1] = -ux[1]+n*(d2*u[3]-a2*u[1]*(1.0-u[2]-u[3]))+(1.0/pe1)*uxx[1]
    fval[2] = a1*u[0]*(1.0-u[2]-u[3])-d1*u[2]-r*u[2]*u[3]*(1.0-u[2]-u[3])**2 \
              +(1.0/pe2)*uxx[2]
    fval[3] = a2*u[1]*(1.0-u[2]-u[3])-d2*u[3]-r*u[2]*u[3]*(1.0-u[2]-u[3])**2 \
              +(1.0/pe2)*uxx[3]
    return fval

# Function defining the left spatial boundary condition.
def bndxa(t, u, ux, bval):
    bval[0] = ux[0]+pe1*(2.0-c-u[0])
    bval[1] = ux[1]+pe1*(c-u[1])
    bval[2] = ux[2]
    bval[3] = ux[3]
    return bval

# Function defining the right spatial boundary condition.
def bndxb(t, u, ux, bval):
    bval[0] = ux[0]
    bval[1] = ux[1]
    bval[2] = ux[2]
    bval[3] = ux[3]
    return bval

# Function defining the initial conditions.
def uinit(x, u):
    u[0] = 2.0-c
    u[1] = c
    u[2] = 0.0
    u[3] = 0.0
    return u 

# Instantiate problem definition object.
problem_definition = bacoli_py.ProblemDefinition(npde, f=f, 
                                            bndxa=bndxa, 
                                            bndxb=bndxb,
                                            uinit=uinit)

# Initial time and mesh.
initial_time = 0
initial_mesh = numpy.linspace(0, 1, 101)

# Output points.
tspan = numpy.linspace(0, 18, 101)
xspan = numpy.linspace(0, 1, 101)

# Set a high level of error control.
atol = 1.0e-6
rtol = atol

solution = solver.solve(problem_definition, initial_time, initial_mesh,
                           tspan, xspan, atol, rtol)

# Get the approximate solution from the Solution object.
u = solution.u

# Plotting these numerical results in 3D.
import matplotlib as mpl
mpl.use('AGG')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

styling = {
    'cmap': cm.coolwarm,
    'linewidth': 0,
    'antialiased': True
}

# Convert xspan and tspan into coordinate arrays for plotting.
T, X = numpy.meshgrid(tspan, xspan)

# Extract the solution for the first PDE in the solved system.
Z = u[0,:,:]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, T, numpy.transpose(Z), **styling)

ax.set_xlabel('$x$')
ax.set_ylabel('$t$')
ax.set_zlabel('$u_1(t,x)$')

ax.view_init(azim=30)

plt.savefig('U1.png')