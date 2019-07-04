# Solving the nonlinear Schrodinger system.
# ------------------------------------------------------------------

import bacoli_py
import numpy

# Create Solver object. Here use the Runge-Kutta method for time
# integration and allow a large number of spatial subintervals to be
# used.
solver = bacoli_py.Solver(t_int='r', nint_max=2000)

# The number of PDEs in this system.
npde = 4

# Initialize problem-dependent parameters.
tempt1 = numpy.sqrt(6.0/5.0)
tempt2 = numpy.sqrt(2.0)

# Function defining the PDE system.
def f(t, x, u, ux, uxx, fval):

    fval[0] = -0.5*ux[0] - 0.5*uxx[1] - u[1]      \
        * ((u[0] * u[0] + u[1] * u[1]) + 2.0/3.0  \
        * ((u[2] * u[2] + u[3] * u[3])))

    fval[1] = - 0.5 * ux[1] + 0.5 * uxx[0] + u[0] \
        * ((u[0] * u[0] + u[1] * u[1]) + 2.0/3.0  \
        * ((u[2] * u[2] + u[3] * u[3])))

    fval[2] = 0.5 * ux[2] - 0.5 * uxx[3] - u[3]   \
        * ((u[2] * u[2] + u[3] * u[3]) + 2.0/3.0  \
        * ((u[0] * u[0] + u[1] * u[1])))

    fval[3] = 0.5 * ux[3] + 0.5 * uxx[2] + u[2]   \
        * ((u[2] * u[2] + u[3] * u[3]) + 2.0/3.0  \
        * ((u[0] * u[0] + u[1] * u[1])))

    return fval

# Function defining the left spatial boundary condition.
def bndxa(t, u, ux, bval):
    bval[0] = ux[0]
    bval[1] = ux[1]
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
    u[0] = tempt1/numpy.cosh(tempt2*x)*numpy.cos(0.5*x)
    u[1] = tempt1/numpy.cosh(tempt2*x)*numpy.sin(0.5*x)
    u[2] = tempt1/numpy.cosh(tempt2*x)*numpy.cos(1.5*x)
    u[3] = tempt1/numpy.cosh(tempt2*x)*numpy.sin(1.5*x)

    return u 

# Instantiate problem definition object.
problem_definition = bacoli_py.ProblemDefinition(npde, f=f, 
                                            bndxa=bndxa, 
                                            bndxb=bndxb,
                                            uinit=uinit)

# Set t_0.
initial_time = 0

# Initial spatial mesh.
initial_mesh = numpy.linspace(-30, 90, 101)

# Output points
tspan = numpy.linspace(0.001, 50, 101)
xspan = numpy.linspace(-30, 90, 101)

# Set a high level of error control.
atol = 1.0e-6
rtol = atol

# Solve the nonlienar Schrodinger system.
evaluation = solver.solve(problem_definition, initial_time, initial_mesh,
                           tspan, xspan, atol, rtol)

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
for i in range(npde):
    Z = evaluation.u[i,:,:]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, T, numpy.transpose(Z), **styling)

    ax.set_xlabel('$x$')
    ax.set_ylabel('$t$')
    ax.set_zlabel('$u_{}(t,x)$'.format(str(i+1)))

    # ax.view_init(azim=-210)

    plt.savefig('U{}.png'.format(str(i+1)))

    plt.clf()