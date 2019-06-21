# Tutorial Example - Modeling the spread of influenza over a spatial domain x
#     over time.
# Model comes from Samzuzzoha et. al 2010.

import bacoli_py
import numpy
from numpy import exp

# Initialize the BacoliPy object with optional parameter nint_max, the maximum
# number of spatial subintervals which will be used in solving this system.
solver = bacoli_py.BacoliPy()

# Initialize problem-dependent parameters.
be = 0.514
gam = 0.20
mu = 5.5e-5
r = 0.0714
ka = 1.0
sig = 0.50
kb = 0.1857
alph = 0.0093
d1 = 0.05
d2 = 0.025
d3 = 0.001
d4 = 0.01

# Function defining the PDE system.
# N = Total population at time t.
#
# S = fval[0] = Number of individuals in host population susceptible to 
#     influenza
#
# E = fval[1] = Number of individuals who have been exposed to influenza and
#     are currently in the latent period of the disease.
#
# I = fval[2] = Number of individuals who are currently infected with influenza.
#
# R = fval[3] = Number of individuals who have recovered from the disease.
def f(t, x, u, ux, uxx, fval):
    N = u[0] + u[1] + u[2] + u[3] 
    fval[0] = -be*u[0]*(u[1] + u[2])/N - mu*u[0] + r*N*(1 - N/ka) + d1*uxx[0]
    fval[1] = be*u[0]*(u[1] + u[2])/N - (mu + sig + kb)*u[1] + d2*uxx[1]
    fval[2] = sig*u[1] - (mu + alph + gam)*u[2] + d3*uxx[2]
    fval[3] = kb*u[1] + gam*u[2] - mu*u[3] + d4*uxx[3]
    return fval

# Function defining the left spatial boundary condition. Homogeneous Neumann
# boundary conditions of form f_{i}(t, u_{i}, u_{i}x) = 0 (i = 1... npde) are
# used at the left spatial boundary.
def bndxa(t, u, ux, bval):
    bval[0] = ux[0]
    bval[1] = ux[1]
    bval[2] = ux[2]
    bval[3] = ux[3]
    return bval

# Function defining the right spatial boundary condition. Homogeneous Neumann
# boundary conditions of form f_{i}(t, u_{i}, u_{i}x) = 0 (i = 1... npde) are
# used at the right spatial boundary.
def bndxb(t, u, ux, bval):
    bval[0] = ux[0]
    bval[1] = ux[1]
    bval[2] = ux[2]
    bval[3] = ux[3]
    return bval

# Function defining the initial conditions.
# Population of susceptible and infected individuals are most concentrated at
# the origin of the spatial domain.
def uinit(x, u):
    u[0] = 0.96*exp(-10*(x**2))
    u[1] = 0.0
    u[2] = 0.04*exp(-100*(x**2))
    u[3] = 0.0
    return u

# Specify the number of PDE's in the system.
npde = 4

# Initialize ProblemDefinition object.
problem_definition = bacoli_py.ProblemDefinition(npde, f=f, 
                                            bndxa=bndxa, 
                                            bndxb=bndxb,
                                            uinit=uinit)

# Set t0.
initial_time = 0.0

# Make initial spatial mesh of 10 uniformly spaced partitions of the spatial 
# domain, with a left boundary of -2 and right boundary of 2.
initial_mesh = numpy.linspace(-2, 2, 11)

# Set the times at which the solution should be output.
tspan = array([0.0, 5.0, 10.0, 15.0, 40.0])

# Set the points at which the solution should be output.
xspan = numpy.linspace(-2,2,100)

# Solve this system for each output time and point.
solution = solver.solve(problem_definition, initial_time, initial_mesh,
                                  tspan, xspan)

# Use matplotlib.pyplot to plot the distribution of the susceptible population 
# at each time.
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

# Plot all all times on same graph.
plt.plot(xspan, solution.u[0,0,:])
plt.plot(xspan, solution.u[0,1,:])
plt.plot(xspan, solution.u[0,2,:])
plt.plot(xspan, solution.u[0,3,:])
plt.plot(xspan, solution.u[0,4,:])

plt.legend(['t=0', 't=5', 't=10', 't=15', 't=40'])

plt.savefig("tutorial_graph.png")
