# Solving 4 One Layer Burgers Equations.
# PDE: u_t = eps*u_xx - u*ux, with initial and boundary conditions
# defined from the exact solution.
# ------------------------------------------------------------------
# This example is based off a FORTRAN analogue for the original 
# BACOLI. The original can be found at:
#   http://cs.stmarys.ca/~muir/BACOLI95-3_Source/3-Problems/burg1.f
# ------------------------------------------------------------------

import bacoli_py
import numpy
from numpy import tanh

# Initialize the BacoliPy object.
solver = bacoli_py.BacoliPy()

# Specify the number of PDE's in this system.
npde = 1

# Initialize problem-dependent parameters.
eps = 1.0e-3

# Function defining the PDE system.
def f(t, x, u, ux, uxx, fval):
    fval[0] = eps*uxx[0] - u[0]*ux[0]

    return fval

# Function defining the left spatial boundary condition.
def bndxa(t, u, ux, bval):
    bval[0] = u[0] - 0.5 + 0.5*tanh( (-0.5*t-0.25) / (4.0*eps) )

    return bval

# Function defining the right spatial boundary condition.
def bndxb(t, u, ux, bval):
    bval[0] = 0.5*tanh((0.75-0.5*t)/(4.0*eps)) - 0.5 + u[0]

    return bval

# Function defining the initial conditions.
def uinit(x, u):
    u[0] = 0.5 - 0.5 * tanh((x - 0.25) / (4.0*eps))

    return u

# Pack all of these callbacks and the number of PDE's into a 
# ProblemDefinition object.
problem_definition = bacoli_py.ProblemDefinition(npde, f=f, 
                                            bndxa=bndxa, 
                                            bndxb=bndxb,
                                            uinit=uinit)

# Specify initial mesh, output_points and output_times.

# Set t0.
initial_time = 0.0

# Define the initial spatial mesh.
initial_mesh = numpy.linspace(0, 1, 11)

# Choose output times and points.
tspan = numpy.linspace(0, 1, 100)
xspan = numpy.linspace(0, 1, 100)

# Solve this problem.
solution = solver.solve(problem_definition, initial_time, initial_mesh,
                           tspan, xspan, atol=1e-6, rtol=1e-6)

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
Z = solution.u[0,:]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(T, X, Z, **styling)

ax.set_xlabel('$x$')
ax.set_ylabel('$t$')
ax.set_zlabel('$u(t,x)$')

plt.savefig('trimesh.png')
