import bacoli_interface
from bacoli_py.ProblemDefinition import ProblemDefinition
from bacoli_py.Evaluation import Evaluation
import numpy as np
from numpy import array
import numbers

class Solver:

    """PDE solver wrapping the BACOLI and BACOLRI software packages.

    Computes an error controlled numerical solution to a given system
    of parabolic 1D partial differential equations. 
    """

    def __init__(self, nint_max=500, kcol=4, t_int='b', s_est='loi', maxord=None,
                 ini_ss=None):
        """A solver object, may be passed optional arguments which determine
           the functionality of the underlying solver. 

        Parameters
        ----------
        nint_max : int
            The maximum number of subintervals allowed for the spatial mesh.
        kcol : int
            The number of collocation points per mesh subinterval which
            determines the order of the piecewise polynomials basis used
            in spatial discretization.
        t_int : string
            Determines which of the two adaptive error control time integration
            schemes be used. t_int = 'b' corresponds to an adaptive order BDF.
            method, t_int = 'r' corresponds to a 5th order Runge-Kutta Method.
            Default value of 'b'.
        s_est : string 
            Determines which of the two spatial error estimation schemes will
            be used, LOI/LE (s_est = 'loi') and SCI/ST (s_est = 'sci'). Default
            value of 'loi'.
        maxord : int
            The maximum order of BDF method to be used in time integration.
            Only used if t_int = 'b'. 1 <= maxord <= 5
        ini_ss : float 
            Initial stepsize for time integration. If not
            provided will be chosen automatically.
        
        Raises
        ------
        ValueError
            If kcol is outside of acceptable range.
        ValueError
            If est is not either 'loi/le' or 'sci/st'
        ValueError
            If maxord is outside of acceptable range.
        ValueError
            If arguments can not be cast into numpy arrays.
        ValueError
            If nint_max <= 0.
        TypeError
            If any arguments are of incorrect type.
        """

        self.bacoli_obj = bacoli_interface.bacoli_interface

        if isinstance(nint_max, int):
            if nint_max <= 0:
                raise ValueError("nint_max must be > 0.")
        else:
            raise TypeError("nint_max must be an integer.")

        if isinstance(kcol, int):
            if kcol <= 2 or kcol >= 11:
                raise ValueError("kcol inside range 2 < kcol < 11.")
        else:
            raise TypeError("kcol must be an integer.")

        if s_est != 'loi' and s_est != 'sci':
            raise ValueError("est may only have values 'loi/le' "
                           + "for LOI/LE scheme or 'sci/st' for SCI/ST "
                           + "scheme.")
        elif s_est == 'loi':
            s_est = 0
        else:
            s_est = 1

        if isinstance(maxord, int):
            if maxord < 1 or maxord > 5:
                raise ValueError("maxord must be in range 0 < maxord < "
                            + "6.")
            else:
                self.is_maxord = True

                if (t_int == 'r'):
                    print("Note: when using Solver with t_int ='r', setting " \
                        + "maxord has no effect.")
        elif maxord == None:
            maxord = 0
            self.is_maxord = False
        else:
            raise TypeError("maxord must be an integer.")

        # If ini_ss not provided, set to default value and handle it on the 
        # Fortran side.
        if ini_ss != None:
            if isinstance(ini_ss, numbers.Number):
                ini_ss    = float(ini_ss)
                self.is_ini_ss = True
            else:
                raise TypeError("ini_ss must be a float.")
        else:
            ini_ss = -1
            self.is_ini_ss = False

        if t_int != 'b' and t_int != 'r':
            raise ValueError("t_int must either be 'b', or 'r'.")
        else:
            # Convert t_int to the form used by the Fortran software.
            if t_int == 'b':
                t_int = 0
            else:
                t_int = 1

        # Convert all arguments into numpy arrays for passing to Fortran.
        try:
            self.nint_max = array(nint_max, dtype=np.int)
            self.kcol = array(kcol, dtype=np.int)
            self.t_int = array(t_int, dtype=np.int)
            self.s_est = array(s_est, dtype=np.int)
            self.maxord = array(maxord, dtype=np.int)
            self.ini_ss = array(ini_ss, dtype=np.float64)
        except ValueError:
            print('Could not convert all Solver arguments into numpy arrays')
            raise

    def solve(self, problem_definition, initial_time, initial_mesh, tspan,
              xspan, atol=1e-4, rtol=1e-4, dirichlet=False, tstop=None,
              compiled_callbacks=False, deriv=False):
        """Solves a system of Partial Differential Equations numerically.

        Parameters
        ----------
        problem_definition : :class:`ProblemDefinition`
            Object containing callback functions which define the problem to be solved.
        initial_time : float 
            Initial point on the time domain.
        initial_mesh : castable to floating point ndarray
            The initial spatial mesh. If it has length 2, then these are assumed to be the
            boundaries of the spatial domain and an initial mesh is automatically generated
            which is adapted to the behaviour of the initial conditions.
        tspan : castable to floating point ndarray
            Vector or scalar containing times at which the solution will be output. 
        xspan : castable to floating point ndarray
            Vector or scalar containing points at which the solution will be output. 
        atol : castable to floating point ndarray
            Absolute error tolerance. Can be either a scalar or a numpy array.
        rtol : castable to floating point ndarray
            Relative error tolerance. Can be either a scalar or a numpy array.
        dirichlet : bool
            Specifies if both of the boundary conditions are dirichlet.
        tstop : float
            Indicates the absolute end of the temporal domain. Used by BACOLI's
            underlying time integrator DASSL. Only used if t_int = 'b'.
        compiled_callbacks   : bool
            Indicates whether compiled callbacks are given in problem_definition.
        deriv : bool
            Indicates that the returned Solution object should contain the first
            spatial derivative at each point.
            
        Returns
        -------
        bacoli_solution : :class:`Evaluation`
            A object containing the results which have been computed with BACOLI.
        
        Raises
        ------
        ValueError
            If any arguments have incorrect values.
        TypeError
            If any arguments are of incorrect type.
        """

        # Check that first argument is a ProblemDefinition object.
        if not isinstance(problem_definition, ProblemDefinition):
            raise TypeError('First argument must be a ProblemDefinition '
                           + 'object containing necissary callback functions.')
        
        # Validate initial time.
        if not isinstance(initial_time, numbers.Number):
            raise TypeError('initial_time must be a scalar quantity.')



        # Make sure initial spatial mesh contains at least the left and right 
        # extents of a spatial domain and if it is of correcy type.
        if not isinstance(initial_mesh, np.ndarray):
            try:
                initial_mesh = np.asarray(initial_mesh, dtype=np.float64)
            except ValueError:
                print('Could not convert initial_mesh into numpy array.')
                raise

        # Validate vectorization flag
        if not isinstance(compiled_callbacks, bool):
            raise ValueError('compiled_callbacks must have type bool.')
                
        # Check that initial_mesh has >= 2 elements
        if len(initial_mesh) < 2:
            raise ValueError('initial_mesh must contain >= 2 elements.')

        # Check for monotone mesh.
        if not all(x<y for x, y in zip(initial_mesh, initial_mesh[1:])):
            raise ValueError('initial_mesh must contain strictly increasing '
                + 'values.')

        # Validate tspan.
        if not isinstance(tspan, np.ndarray):
            try:
                tspan = np.asarray(tspan, dtype=np.float64)
            except ValueError:
                print('Could not convert tspan into a numpy array.')
                raise

        if tspan.size != 1:
            if not all(x<y for x, y in zip(tspan, tspan[1:])):
                raise ValueError('tspan must contain strictly increasing '
                    + 'values.')

        # Validate xspan.
        if not isinstance(xspan, np.ndarray):
            try:
                xspan = np.asarray(xspan, dtype=np.float64)
            except ValueError:
                print('Could not convert xspan into a numpy array.')
                raise

        if not all(x<y for x, y in zip(xspan, xspan[1:])):
            raise ValueError('xspan must contain strictly increasing '
                + 'values.')

        npde = problem_definition.npde

        # Check if atol is provided and in an acceptable form. If it is not
        # provided, will be set to default value when the call to FORTRAN is
        # made.
        if isinstance(atol, numbers.Number):
            atol = np.asarray(atol, dtype=np.float64)
        elif not isinstance(atol, np.ndarray):
            try:
                atol = np.asarray(atol, dtype=np.float64) 
            except ValueError:
                print('Could not convert atol into numpy array.')
                raise

            if len(atol) != npde:
                raise ValueError('atol must either be a scalar or a numpy array '
                            + 'of size npde.')
        else:
            if len(atol) != npde:
                raise ValueError('atol must either be a scalar or a numpy array '
                            + 'of size npde.')

        # Check if rtol is provided and in an acceptable form. If it is not
        # provided, will be set to default value when the call to FORTRAN is
        # made.
        if isinstance(rtol, numbers.Number):
            rtol = np.asarray(rtol, dtype=np.float64)
        elif not isinstance(rtol, np.ndarray):
            try:
                rtol = np.asarray(rtol, dtype=np.float64) 
            except ValueError:
                print('Could not convert rtol into numpy array.')
                raise

            if len(rtol) != npde:
                raise ValueError('rtol must either be a scalar or a numpy array '
                            + 'of size npde.')
        else:
            if len(rtol) != npde:
                raise ValueError('rtol must either be a scalar or a numpy array '
                              + 'of size npde.')

        # Validate dirichlet and convert to acceptable form.
        if dirichlet == False:
            dirichlet = 0 
        elif dirichlet == True:
            dirichlet = 1
        else:
            raise TypeError('dirichlet must be a boolean value.')

        is_tstop = None
        # Validate tstop.
        if tstop == None:
            tstop = -1 
            is_tstop = False
        else:
            # Check if tstop is a number by attempting to convert it to float.
            if not isinstance(tstop, numbers.Number):
                raise TypeError('tstop must be a number.')
            else:
                is_tstop = True
                if (self.t_int == 1):
                    print("Note: when using Solver with t_int ='r', setting " \
                        + "tstop has no effect.")
        
        # Initialize Bacoli95 solver object.
        idid = self.bacoli_obj.initialize(npde=problem_definition.npde,
            nint_max=self.nint_max, kcol=self.kcol, s_est=self.s_est,
            t_est=self.t_int, maxord=self.maxord, is_maxord=self.is_maxord,
            atol=atol, rtol=rtol, t0=initial_time, x=initial_mesh,
            xspan=xspan, dirichlet=dirichlet, ini_ss=self.ini_ss,
            is_ini_ss=self.is_ini_ss, tstop=tstop,
            is_tstop=is_tstop)

        # Check that initialization was successful.
        if idid < 0:
            raise RuntimeError(__get_error_message(self.bacoli_obj.idid))

        # Memory used to contain array slices in convernient format when
        # vectorization is to be used
        if not compiled_callbacks:
            u_sliced    = np.ndarray(shape=(npde), dtype=np.ndarray)
            ux_sliced   = np.ndarray(shape=(npde), dtype=np.ndarray)
            uxx_sliced  = np.ndarray(shape=(npde), dtype=np.ndarray)
            fval_sliced = np.ndarray(shape=(npde), dtype=np.ndarray)

            def __fvec(t, x, u, ux, uxx, fval, npde, vnpts):
                """Function for vectorized BACOLI"""

                for i in range(int(npde)):
                    l = i*vnpts
                    r = (i+1)*vnpts

                    u_sliced[i]    = u[l:r]
                    ux_sliced[i]   = ux[l:r]
                    uxx_sliced[i]  = uxx[l:r]
                    fval_sliced[i] = fval[l:r]


                # Make vectorized call with f.
                problem_definition.f(t, x, u_sliced, ux_sliced, 
                    uxx_sliced, fval_sliced)

                for i in range(int(npde)):
                    fval[i*vnpts:(i+1)*vnpts] = fval_sliced[i]

                # Return fval, array will have been modified with this call.
                return fval
        else:
            def __fvec(t, x, u, ux, uxx, fval, npde, vnpts):
                pass

        # Create arrays to be returned to user. Values for calls to 
        # bacoli95_vals are appended to this array and final result is returned
        # to caller.
        solution = np.empty(shape=(npde, tspan.size, xspan.size))

        # Array containing derivative information
        solution_deriv = np.empty(shape=(npde, tspan.size, xspan.size))

        # Create array used to hold values from calls to bacoli95_vals.
        u_in = np.empty(shape=(npde, xspan.size))

        # Set nderiv based on whether the derivative valuse was requested.
        nderiv = 1 if deriv else 0

        # Main loop
        for i in range(tspan.size):
            # If an output time is equal to the initial time, then call
            # the function specifying initial conditions. 
            if tspan[i] == initial_time:
                raise ValueError('tspan can not include the initial point in time.')

            # Call the bacoli95 subroutines interfacing subroutine. Advances the
            # solution to the next value in tspan.
            idid = self.bacoli_obj.solve(tout=tspan[i],
                f=problem_definition.f,
                fvec=__fvec,
                bndxa=problem_definition.bndxa,
                bndxb=problem_definition.bndxb,
                uinit=problem_definition.uinit,
                derivf=problem_definition.derivf,
                is_derivf=problem_definition.is_derivf,
                difbxa=problem_definition.difbxa,
                is_difbxa=problem_definition.is_difbxa,
                difbxb=problem_definition.difbxb,
                is_difbxb=problem_definition.is_difbxb,
                vec=(not compiled_callbacks))

            if idid < 0:
                raise RuntimeError("An error was raised during this computation.\n ")

            # Evaluate the solution at this point.
            solution[:,i,:], solution_deriv[:,i,:] \
                = self.bacoli_obj.vals(npde=npde, xspan_size=len(xspan), nderiv=nderiv)

            if idid < 0:
                raise RuntimeError("An error was raised while evaluating the solution.\n ")

        # Free the memory on the Fortran side.
        self.bacoli_obj.sol_teardown()

        # Build and return solution object.
        if deriv:
            return Evaluation(tspan, xspan, solution, solution_deriv)
        else:
            return Evaluation(tspan, xspan, solution)

    @staticmethod
    def __get_error_message(idid):
        """Returns error message corresponding to error indicator set by BACOLI.
        """

        return {
            '-1000' : 'Memory allocation error occured.',
       
        }.get(idid, 'Undefined Error: Please report to developers.')
