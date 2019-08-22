# Two Layer Burgers Equation
# ------------------------------------------------------------------
# This example is based off a FORTRAN analogue for the original 
# BACOLI. The original can be found at:
#   http://cs.stmarys.ca/~muir/BACOLI95-3_Source/3-Problems/burg1.f
# ------------------------------------------------------------------

import bacoli_py
import numpy
from numpy import tanh, exp

# Initialize the Solver object.
solver = bacoli_py.Solver()

# Specify the number of PDE's in this system.
npde = 1

# Initialize problem-dependent parameters.
eps = 1.0e-4

# Function defining the PDE system.
def f(t, x, u, ux, uxx, fval):
    fval[0] = eps*uxx[0] - u[0]*ux[0]

    return fval

# Function defining the left spatial boundary condition.
def bndxa(t, u, ux, bval):
    a1 = (0.5 - 4.95*t)*0.5e-1/eps
    a2 = (0.5 - 0.75*t)*0.25/eps
    a3 = 0.1875/eps
    expa1 = 0.0
    expa2 = 0.0
    expa3 = 0.0
    temp = max(a1, a2, a3)

    if ((a1 - temp) >= -35.0):
        expa1 = exp(a1 - temp)
    if ((a2 - temp) >= -35.0):
        expa2 = exp(a2-temp)
    if ((a3-temp) >= -35.0):
        expa3 = exp(a3-temp)

    bval[0] = u[0] - (0.1*expa1+0.5*expa2+expa3) \
                / (expa1+expa2+expa3)

    return bval

# Function defining the right spatial boundary condition.
def bndxb(t, u, ux, bval):
    a1 = (-0.5 - 4.95*t)*0.5e-1/eps
    a2 = (-0.5 - 0.75*t)*0.25/eps
    a3 = -0.3125/eps
    expa1 = 0.0
    expa2 = 0.0
    expa3 = 0.0
    temp = max(a1, a2, a3)

    if ((a1 - temp) >= -35.0):
        expa1 = exp(a1 - temp)
    if ((a2 - temp) >= -35.0):
        expa2 = exp(a2 - temp)
    if ((a3 - temp) >= -35.0):
        expa3 = exp(a3 - temp)

    bval[0] = u[0] - (0.1*expa1+0.5*expa2+expa3) \
                 / (expa1+expa2+expa3)

    return bval

# Function defining the initial conditions.
def uinit(x, u):
    a1 = (-x + 0.5)*0.5e-1/eps
    a2 = (-x + 0.5)*0.25/eps
    a3 = (-x + 0.375)*0.5/eps
    expa1 = 0.0
    expa2 = 0.0
    expa3 = 0.0
    temp = max(a1, a2, a3)

    if ((a1-temp) >= -35.0):
        expa1 = exp(a1-temp)
    if ((a2-temp) >= -35.0):
        expa2 = exp(a2-temp)
    if ((a3-temp) >= -35.0):
        expa3 = exp(a3-temp)

    u[0] = (0.1*expa1+0.5*expa2+expa3)/(expa1+expa2+expa3)

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
initial_mesh = [0, 1]

# Choose output times and points.
tspan = numpy.linspace(0.001, 1, 100)
xspan = numpy.linspace(0, 1, 100)

# Solve this problem.
evaluation = solver.solve(problem_definition, initial_time, initial_mesh,
                           tspan, xspan, atol=1e-6, rtol=1e-6)

u = evaluation.u[0,:,:]

# Plot the solution.
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

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(T, X, u, **styling)

ax.set_xlabel('$x$')
ax.set_ylabel('$t$')
ax.set_zlabel('$u(t,x)$')

plt.savefig('trimesh.png')

