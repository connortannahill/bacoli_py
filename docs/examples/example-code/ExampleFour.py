# Solving PDE u_t = u_xx + pi**2 * sin(pi*x) (npde=1), u(0,t)=1,
# u(1,t)=1, u(x,0) = 1.
# --------------------------------------------------------------------
# Here we specify Jacobian matrices for the PDE system, and the 
# boundary conditions. We do this by defining and providing
# the callback functions derivf, bndxa, bndxb. Doing so is known
# to increase the speed of computation drastically, as these matrices
# will not have to be calculated in the software.
#
# Additionally we demonstrate a useful way of outputting the
# numerical solutions from bacoli_py to an output file and later
# use this information to generate a 3D plot.
# --------------------------------------------------------------------
# This example is based off a FORTRAN analogue for the original 
# BACOLI. The original can be found at:
#   http://cs.stmarys.ca/~muir/BACOLI95-3_Source/3-Problems/sincmads.f
# --------------------------------------------------------------------

import bacoli_py as bacoli_py
import numpy
from numpy import pi
from numpy import sin

solver = bacoli_py.Solver()
npde = 1

# Set t0.
initial_time = 0.0

# Function defining the PDE system.
def f(t, x, u, ux, uxx, fval):
    fval[0] = uxx[0] + pi * pi * sin(pi * x)
    return fval

# Functioning defining the analytic Jacobian for the PDE system.
def derivf(t, x, u, ux, uxx, dfdu, dfdux, dfduxx):
    dfdu[0,0] = 0.0
    dfdux[0,0] = 0.0
    dfduxx[0,0] = 1.0

    return dfdu, dfdux, dfduxx

# Function defining the analytic Jacobian for the left boundary
# condition.
def difbxa(t, u, ux, dbdu, dbdux, dbdt):
    dbdu[0,0] = 1.0
    dbdux[0,0] = 0.0
    dbdt[0] = 0.0

    return dbdu, dbdux, dbdt

# Function defining the analytic Jacobian for the right boundary
# condition.
def difbxb(t, u, ux, dbdu, dbdux, dbdt):
    dbdu[0,0] = 1.0
    dbdux[0,0] = 0.0
    dbdt[0] = 0.0

    return dbdu, dbdux, dbdt

# Function defining the left spatial boundary condition.
def bndxa(t, u, ux, bval):
    bval[0] = u[0] - 1.0

    return bval

# Function defining the right spatial boundary condition.
def bndxb(t, u, ux, bval):
    bval[0] = u[0] - 1.0

    return bval

# Function defining the initial conditions.
def uinit(x, u):
    u[0] = 1.0

    return u

# Define the ProblemDefinition object, prodividing the optional
# arguments derivf, difbxa, and difbxb along with the usual parameters.
problem_definition = bacoli_py.ProblemDefinition(npde, f=f, 
                                            bndxa=bndxa, 
                                            bndxb=bndxb,
                                            uinit=uinit,
                                            derivf=derivf,
                                            difbxa=difbxa,
                                            difbxb=difbxb)

initial_time = 0
initial_mesh = numpy.linspace(0, 1, 11)
tspan = numpy.linspace(0.001, 1, 20)
xspan = numpy.linspace(0,1,100)
atol = 1.0e-6
rtol = atol

# Solve this system, passing the optional parameter tstop, the absolute
# end of the temporal domain which can be helpful when performing 
# time-integration.
evaluation = solver.solve(problem_definition, initial_time, initial_mesh,
                              tspan, xspan, atol, rtol, tstop=1.0)

# Get numerical solution from BacoliSolution object.
u = evaluation.u

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

# Output numerical solution to output files.
for i in range(npde):

    fname = "Points{0:d}".format(i+1)

    output_file = open(fname, 'w')

    # Print to file.
    for j in range(len(tspan)):
        for k in range(len(xspan)):
            out_str = '{0: <20.17f} {1: <20.17f} {2: <20.17f}\n'.format(
                xspan[k], tspan[j], u[i,j,k])

            output_file.write(out_str)

# Now we use this output to generate a trimesh.
X, T, U = numpy.loadtxt('Points1', unpack=True)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(X, T, U, **styling)

ax.set_xlabel('$x$')
ax.set_ylabel('$t$')
ax.set_zlabel('$u(t,x)$')

plt.savefig('trimesh.png')