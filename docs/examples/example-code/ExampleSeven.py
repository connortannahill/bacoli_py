# Solving Two Layer Burgers Equations.
# PDE: u_t = eps*u_xx - u*ux, with initial and boundary conditions
# defined from the exact solution.
# 
# This example uses compiled callback functions to increase speed.
# ------------------------------------------------------------------
# This example is based off a FORTRAN analogue for the original 
# BACOLI. The original can be found at:
#   http://cs.stmarys.ca/~muir/BACOLI95-3_Source/3-Problems/burg2.f
# ------------------------------------------------------------------

"""
This file has two parts:
    a) Generation of an extension module containing compiled Fortran callback
       subroutines for use with bacoli_py.
    b) Application of bacoli_py to solve a problem using these compiled routines.
"""

"""
Part 1

Creating and linking the Fortran callback functions
"""
import numpy.f2py as f2py
import sys

prob_def_f = """
      subroutine f(t, x, u, ux, uxx, fval)
          integer          npde
          parameter       (npde=1)
          double precision t, x, u(npde), ux(npde)
          double precision uxx(npde), fval(npde)

          double precision eps
          parameter       (eps=1d-4)

          fval(1) = eps*uxx(1) - u(1)*ux(1)
      return
      end

      subroutine bndxa(t, u, ux, bval)
          integer          npde
          parameter        (npde=1)
          double precision t, u(npde), ux(npde), bval(npde)
          double precision eps
          parameter       (eps=1d-4)
          double precision a1, a2, a3, expa1, expa2, expa3, temp

          a1 = (0.5d0 - 4.95d0 * t) * 0.5d-1 / eps
          a2 = (0.5d0 - 0.75d0 * t) * 0.25d0 / eps
          a3 = 0.1875d0 / eps
          expa1 = 0.d0
          expa2 = 0.d0
          expa3 = 0.d0
          temp = max(a1, a2, a3)
          if ((a1-temp) .ge. -35.d0) expa1 = exp(a1-temp)
          if ((a2-temp) .ge. -35.d0) expa2 = exp(a2-temp)
          if ((a3-temp) .ge. -35.d0) expa3 = exp(a3-temp)

          bval(1) = u(1) - (0.1d0*expa1+0.5d0*expa2+expa3)
     &                / (expa1+expa2+expa3)
      return
      end

      subroutine bndxb(t, u, ux, bval)
          integer          npde
          parameter       (npde=1)
          double precision t, u(npde), ux(npde), bval(npde)
          double precision eps
          parameter       (eps=1d-4)
          double precision a1, a2, a3, expa1, expa2, expa3, temp

          a1 = (-0.5d0 - 4.95d0 * t) * 0.5d-1 / eps
          a2 = (-0.5d0 - 0.75d0 * t) * 0.25d0 / eps
          a3 = - 0.3125d0 / eps
          expa1 = 0.d0
          expa2 = 0.d0
          expa3 = 0.d0
          temp = max(a1, a2, a3)
          if ((a1-temp) .ge. -35.d0) expa1 = exp(a1-temp)
          if ((a2-temp) .ge. -35.d0) expa2 = exp(a2-temp)
          if ((a3-temp) .ge. -35.d0) expa3 = exp(a3-temp)

          bval(1) = u(1) - (0.1d0*expa1+0.5d0*expa2+expa3)
     &                / (expa1+expa2+expa3)
      return
      end

      subroutine uinit(x, u)
          integer          npde
          parameter        (npde=1)
          double precision x, u(npde)
          double precision eps
          parameter       (eps=1d-4)
          double precision a1, a2, a3, expa1, expa2, expa3, temp

          a1 = (-x + 0.5d0) * 0.5d-1 / eps
          a2 = (-x + 0.5d0) * 0.25d0 / eps
          a3 = (-x + 0.375d0) * 0.5 / eps
          expa1 = 0.d0
          expa2 = 0.d0
          expa3 = 0.d0
          temp = max(a1, a2, a3)
          if ((a1-temp) .ge. -35.d0) expa1 = exp(a1-temp)
          if ((a2-temp) .ge. -35.d0) expa2 = exp(a2-temp)
          if ((a3-temp) .ge. -35.d0) expa3 = exp(a3-temp)

          u(1) = (0.1d0*expa1+0.5d0*expa2+expa3) / (expa1+expa2+expa3)
      return
      end
"""

f2py.compile(prob_def_f.encode('ascii'), modulename='problemdef', verbose=0)

"""
Part 2

Using the compiled callback routines to solve the problem.
"""

import bacoli_py
import numpy
import time
from problemdef import f, bndxa, bndxb, uinit

# Initialize the Solver object.
solver = bacoli_py.Solver()

# Specify the number of PDE's in this system.
npde = 1

# Pack all of these callbacks and the number of PDE's into a 
# ProblemDefinition object.
problem_definition = bacoli_py.ProblemDefinition(npde, f=f._cpointer, 
                                            bndxa=bndxa._cpointer,
                                            bndxb=bndxb._cpointer,
                                            uinit=uinit._cpointer)

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
                           tspan, xspan, vec=False)

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
Z = evaluation.u[0,:]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(T, X, Z, **styling)

ax.set_xlabel('$x$')
ax.set_ylabel('$t$')
ax.set_zlabel('$u(t,x)$')

plt.savefig('trimesh.png')
