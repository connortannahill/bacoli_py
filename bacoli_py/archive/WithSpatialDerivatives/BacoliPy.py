"""
    Class which contains the functionality of the BACOLI_PY python wrapper.
"""

import bacoli95_interface
from ProblemDefinition import ProblemDefinition
from BacoliSolution import BacoliSolution
import numpy
from numpy import array

class BacoliPy:

    """
    A solver object.
    @author: Connor Tannahill
    """

    def __init__(self, nint_max=None, kcol=None, estimator=None, maxord=None,
                 ini_ss=None):
        """A BACOLI solver object. 
        
        Takes as arguments many of the optional parameters controlling
        the underlying functionality of the BACOLI PDE solver.

        Parameters:
            nint_max (int, optional): The maximum number of subintervals allowed 
                along the spatial mesh.

            kcol (int, optional): The number of collocation points per mesh 
                subinterval and determines the order of the piecewise 
                polynomials.

            estimator (int, optional): Determines which of the two spatial error 
                estimation schemes will be used, LOI/LE (estimator = 0) and
                SCI/ST (estimator = 1).

            maxord (int, optional): The maximum order of BDF method used by 
                BACOLI's underlying time integrator DASSL.

            ini_ss (float, optional): Initial stepsize for underlying time
                integrator DASSL.
        """

        # Check for optional arguments, assign unpassed arguments to -1.
        if nint_max is not None:
            self.nint_max = nint_max
        else:
            self.nint_max = -1

        if kcol is not None:
            self.kcol = kcol
        else:
            self.kcol = -1

        if estimator is not None:
            self.estimator = estimator
        else:
            self.estimator = -1

        if maxord is not None:
            self.maxord= maxord
        else:
            self.maxord = -1
        
        if ini_ss is not None:
            self.ini_ss = ini_ss
        else:
            self.ini_ss = -1

        # Convert all arguments into numpy arrays for passing to Fortran.
        try:
            self.nint_max = array(self.nint_max, dtype=numpy.int)
            self.kcol = array(self.kcol, dtype=numpy.int)
            self.estimator = array(self.estimator, dtype=numpy.int)
            self.maxord = array(self.maxord, dtype=numpy.int)
            self.ini_ss = array(self.ini_ss, dtype=numpy.float)
        except ValueError:
            raise ValueError('Could not convert all Bacolipy arguments '
                             + 'into numpy arrays')

    def intialize_problem(self, problem_definition, initial_time,
              initial_mesh, atol=None, rtol=None, dirichlet=None):
        """ Initializes the BACOLI95 object for this problem.
            
            Args:
                problem_definition (:obj:`ProblemDefinition`): ProblemDefinition 
                    object containing callback functions which define the problem 
                    to be solved. Also contains the number of PDE's.

                initial_time (float): Initial output time.

                initial_mesh (:obj:`list` of (float)): The initial spatial mesh. 
                    If user provides array of size 2 then it is assumed to contain 
                    the left and right spatial endpoints and a spatial mesh with 10
                    subintervals is constructed between these two points.

                atol (:obj:`list` of (float), optional): Absolute error tolerance, 
                    optional. Can be scalar or numpy array.

                rtol (:obj:`list` of (float), optional): Relative error tolerance, 
                    optional. Can be scalar or numpy array.

                dirichlet ((int), optional): Specifies if either of the boundary
                    conditions are dirichlet.
        """

        bacoli_obj = bacoli95_interface.bacoli95_interface

        # Add tstart as private instance variable of BacoliPy object.
        self.initial_time = initial_time

        # Check that first argument is a ProblemDefinition Object.
        if not isinstance(problem_definition, ProblemDefinition):
            raise ValueError('Error: First argument must be a ProblemDefinition '
                             + 'object containing necissary callback functions.')
        else:
            self.problem_definition = problem_definition

        # Check for presence of absolute and relative error tolerences and 
        # convert then to numpy arrays for passing to FORTRAN.
        try:
            if atol is not None:
                if isinstance(atol, (float, int)):
                    atol = array(atol, dtype=numpy.float)
                is_atol = True
            else:
                atol = array(-1, dtype=numpy.float)
                is_atol = False

            # Check if rtol is provided.
            if rtol is not None:
                if isinstance(rtol, (float, int)):
                    rtol = array(rtol, dtype=numpy.float)
                is_rtol = True
            else:
                rtol = array(-1, dtype=numpy.float)
                is_rtol = False
        except ValueError:
            raise ValueError('Error: Could not convert absolute and relative '
                             + 'error tolerances to numpy arrays. atol and '
                             + 'rtol must be given as either scalar '
                             + 'quantities, lists, or numpy arrays.')

        # Make sure initial spatial mesh contains at least the left and right 
        # extents of a spatial domain.
        if len(initial_mesh) < 2:
            raise ValueError('Initial spatial mesh must at least contain two '
                             + 'spatial boundary points.')

        # Check for integer arguments, if not present set them to a default value.
        if dirichlet is None:
            dirichlet = -1

        # Convert all arguments to numpy arrays for passing to FOTRAN.
        try:
            dirichlet = array(dirichlet, dtype=numpy.int)
        except ValueError:
            raise ValueError('Could not cast all integer arguments '
                             + 'to numpy arrays.')

        # Initialize Bacoli95 solution object.
        bacoli_obj.bacoli95_init_wrap(npde=problem_definition.npde, 
                        x=initial_mesh, tstart=initial_time, atol=atol, 
                        is_atol=is_atol, rtol=rtol, is_rtol=is_rtol, 
                        kcol=self.kcol, nint_max=self.nint_max, 
                        estimator=self.estimator, dirichlet=dirichlet, 
                        maxord=self.maxord, ini_ss=self.ini_ss)
    
    def solve(self, output_times, output_points, tstop=None,
              nderiv=None, nsteps=None):
        """Method for returning BACOLI95 solution values.

        Makes call to BACOLI, computing an error controlled numerical solution
        for the PDE system.

        Args:
            output_times (:ob:`list` of (int)): Vector or scalar containing 
                all of times the solution will be outputted at. 

            output_points (:ob:`list` of (int)): List containing x values 
                solution is to be output at If not provided the x values used 
                are from the final spatial mesh after BACOLI has performed its 
                adaptive error control.
            
            tstop (float, optional): Indicates the absolute end of the temporal  
                domain. Used by BACOLI's underlying time integrator DASSL.

            nderiv (float, optional): Number of spatial partial derivatives
                whole values to be output are given. Optional

            nsteps (bool, optional): Indicator for if the number of time steps 
                performed by DASSL should be returned. Optional.
        
        Returns:
            bacoli_solution (:obj:`BacoliSolution`): A object containing the
                results which have been computed with BACOLI.
        """

        bacoli_obj = bacoli95_interface.bacoli95_interface

        if output_points is list:
            output_points = numpy.asarray(ouput_points, dtype=numpy.float)

        if tstop is None:
            tstop = -1 

        if nderiv is not None:
            is_nderiv = True
            uout_rank = 3
        else:
            nderiv = -1
            is_nderiv = False
            uout_rank = 2
        
        if nsteps is None:
            nsteps = -1
        
        # Convert all scalar arguments into numpy arrays.
        try:
            tstop = array(tstop, dtype=numpy.float)
            nderiv = array(nderiv, dtype=numpy.int)
            uout_rank = array(uout_rank, dtype=numpy.int)
            nsteps = array(nsteps, dtype=numpy.int)
        except ValueError:
            raise ValueError('Could not cast all integer arguments '
                             + 'to numpy arrays.')
        
        npde = self.problem_definition.npde
        
        # Create list to be returned to user. Values for calls to 
        # bacoli95_vals are appended to this array and final result is returned
        # to caller.
        solution = []
        for i in range(npde):
            solution.append([])

        # Create array used to hold values from calls to bacoli95_vals. May
        # be of rank 2 or 3.
        if uout_rank == 2:
            u_in = numpy.empty(shape=(npde, output_points.size))
            spatial_derivatives = None 
        else:
            # For third shape component, we use nderiv + 1 due to the way in
            # which f2py returns the result of the bacoli95_vals_rank2_wrap
            # call.
            u_in = numpy.empty(shape=(npde, output_points.size, nderiv + 1))
            spatial_derivatives = []
            for i in range(npde):
                spatial_derivatives.append([])

        for tout in output_times:

            # If an output time is equal to the initial time, then call
            # the initial condition
            if tout == initial_time:
                for x in output_times:
                    u = numpy.zeros(npde)
                    u = self.problem_definition.initial_conditions(x, u, npde)

                    if uout_rank == 2:
                        for i in range(npde):
                            solution[i].append(u[i])
                    else:
                        for i in range(npde):
                            solution[i].append(u[i])
                            spatial_derivatives[i] = numpy.zeros(len(xout), nderiv)
                    
            # Call the bacoli95 subroutines interfacing subroutine.
            bacoli_obj.bacoli95_wrap(output_time=tout, 
                        f=self.problem_definition.pde_system,
                        bndxa=self.problem_definition.left_boundary_conditions,
                        bndxb=self.problem_definition.right_boundary_conditions,
                        uinit=self.problem_definition.initial_conditions,
                        derivf=self.problem_definition.derivf,
                        is_derivf=self.problem_definition.is_derivf,
                        difbxa=self.problem_definition.difbxa,
                        is_difbxa=self.problem_definition.is_difbxa,
                        difbxb=self.problem_definition.difbxb,
                        is_difbxb=self.problem_definition.is_difbxb,
                        tstop=tstop, nsteps=nsteps)
            
            # Make call to different bacoli95_vals wrapping functions depending
            # on the desired rank of uout.
            if uout_rank == 2:
                # Temporary array to hold value from this call.
                temp_soln = bacoli_obj.bacoli95_vals_rank2_wrap(uout=u_in,
                                xout=output_points, nderiv=nderiv)
                
                # Build solution array.
                for i in range(npde):
                    solution[i].append(temp_soln[i])
            else:
                # Temporary array to hold value from this call.
                temp_soln = bacoli_obj.bacoli95_vals_rank3_wrap(uout=u_in,
                                xout=output_points, nderiv=nderiv)
                
                # Build solution array.
                for i in range(npde):
                    uout = []
                    for j in range(len(output_points)):
                        uout.append(temp_soln[i][j][0])
                    solution[i].append(uout)
                
                # Build array to contain spatial partial derivatives.
                deriv_arr = []
                for i in range(nderiv):
                    deriv_arr.append([])
              
                for i in range(npde):
                    for j in range(len(output_points)):
                        for k in range(nderiv):
                            diriv_arr[k].append(temp_soln[i][j][k + 1])
                    spatial_derivatives[i].append(deriv_arr)
                
        # Build and return solution object.
        bacoli_info = self.__extract_bacoli_info(bacoli_obj)

        bacoli_solution = BacoliSolution(output_points, output_times, solution,
                                         bacoli_info, spatial_derivatives)
        return bacoli_solution

    def __extract_bacoli_info(self, bacoli_obj):
        """Method to extract BACOLI information.

        Extracts all the information forming a BACOLI95_SOL object from
        our wrapping FORTRAN module.

        Args:
            bacoli_obj (:obj:`BACOLI95 object`): Object interfacing with the
            FORTRAN95 wrapping module.

        Returns:
            bacoli_info (:obj:`dict`): Dictionary containing information about
                the computation performed by BACOLI.
        """
        bacoli_obj = bacoli95_interface.bacoli95_interface

        # Initialize BACOLI info.
        bacoli_obj.unpack_sol_parameters()

        bacoli_info = {
            "rtol" : bacoli_obj.rtol,
            "atol" : bacoli_obj.atol,
            "final_mesh" : bacoli_obj.x,
            "bspline_coefficients" : bacoli_obj.y,
            "npde" : bacoli_obj.npde,
            "nint_max" : bacoli_obj.nint_max,
            "kcol" : bacoli_obj.kcol,
            "estimator" : bacoli_obj.estimator,
            "maxord" : bacoli_obj.maxord,
            "ini_ss" : bacoli_obj.ini_ss,
            "t0" : bacoli_obj.t0,
            "nint" : bacoli_obj.nint,
            "mflag" : bacoli_obj.mflag,
            "idid" : bacoli_obj.idid,
            "num_remeshings" : bacoli_obj.num_remeshings,
            "num_ini_remeshings" : bacoli_obj.num_ini_remeshings,
            "num_cold_restarts" : bacoli_obj.num_cold_restarts,
            "num_accepted_time_steps" : bacoli_obj.num_accepted_time_steps,
            "prev_bdf_order" : bacoli_obj.prev_bdf_order,
            "prev_time_step_size" : bacoli_obj.prev_time_step_size
        }

        return bacoli_info