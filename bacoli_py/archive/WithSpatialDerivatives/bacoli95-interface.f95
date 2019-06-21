! A simple interface for bacoli95.f95 to facilitate its use with a python
! wrapper. Calls bacoli95 subroutines based on user defined specifications
! in order to solve a system of 1D Parabolic PDE's. Stores solution variables
! to allow easy retrieval using the python driver.
! @author Connor Tannahill

module bacoli95_interface
    use bacoli95_mod
    implicit none

    ! Globally scoped variable of bacoli95 sol type to allow easy
    ! manipulation of variables and extraction of data to python wrapper.
    type(bacoli95_sol)   :: sol

    ! Globally scoped variable of bacoli95 flags type to send to subroutines
    ! with optional parameters. Indicates to these subroutines whether or not
    ! optional parameters sent should be considered.
    type(bacoli95_flags) :: flags

    ! Globally scoped variables for easy access with a python script. Necissary
    ! due to f2py's incompatability with Fortran 90 user defined types. The 
    ! values for these are unpacked from the bacoli95_sol type when calls to
    ! unpack_sol_parameters are made. Includes spline solutions as well, 
    ! unpacked with unpack_spline_parameters. 

    ! SOL PARAMETERS:
    ! -----------------------------------------------------------------------
    integer           :: npde
    !                        npde is the number of partial differential
    !                        equations in the system being solved.
    !                        Set by bacoli95_init.
    integer           :: nint_max
    !                        A maximum number of spatial mesh
    !                        subintervals is used to (indirectly) limit
    !                        the amount of work and storage employed by
    !                        BACOLI. This is that maximum.
    !                        Set by bacoli95_init.
    integer          :: kcol
    !                        BACOLI uses collocation methods. kcol is
    !                        the number of collocation points per mesh
    !                        subinterval and determines the order of
    !                        piecewise polynomials used.
    !                        Set by bacoli95_init.
    integer           :: estimator
    !                        BACOLI has two options for spatial error
    !                        estimation schemes, LOI/LE (estimator=0)
    !                        and SCI/ST (estimator=1).
    !                        Set by bacoli95_init.
    integer           :: maxord
    !                        The underlying time integrator, DASSL, uses
    !                        multistep BDF methods. maxord may be used
    !                        to limit the maximum order of BDF method
    !                        used, which can be useful since higher
    !                        order methods have smaller stability
    !                        regions.
    !                        Set by bacoli95_init.
    real(kind=8)          :: ini_ss
    !                        The user may optionally specify an initial
    !                        stepsize for the underlying time integrator
    !                        to attempt. This is that initial stepsize.
    !                        Set by bacoli95_init.
    real(kind=8), allocatable :: atol(:)
    !                        The absolute error tolerance, either
    !                        of length 1 or length npde. BACOLI
    !                        and DASSL work to (independently)
    !                        control errors associated with
    !                        spatial discretization and time
    !                        stepping using these tolerances.
    real(kind=8), allocatable :: rtol(:)
    !                        The relative error tolerance, either
    !                        of length 1 or length npde. BACOLI
    !                        and DASSL work to (independently)
    !                        control errors associated with
    !                        spatial discretization and time
    !                        stepping using these tolerances.
    !
    ! variable problem parameters (input-output)
    real(kind=8)          :: t0
    !                        BACOLI constructs B-spline interpolants
    !                        of solution approximations in space for
    !                        specific values of the time variable, t.
    !                        t0 is the starting value of t on input,
    !                        and on output (after a call to bacoli95)
    !                        it represents the value of t associated
    !                        with the current B-spline interpolant.
    !                        Set by bacoli95_init and bacoli95.
    real(kind=8), allocatable :: x(:)
    !                        On input, x is the initial spatial mesh on
    !                        which to attempt time stepping. As BACOLI
    !                        adapts the mesh as part of its error
    !                        control, this is modified on output.
    !                        Set by bacoli95_init and bacoli95.
    integer           :: nint
    !                        nint is the number of subintervals in the
    !                        spatial mesh, x. As BACOLI adapts the
    !                        mesh as part of its error control, this
    !                        number is variable.
    !                        Set by bacoli95_init and bacoli95.
    integer           :: mflag(12)
    !                        A vector of control flags for BACOLI.
    !                        These flags are handled by this wrapper
    !                        module. Consult BACOLI (bacoli.f)
    !                        documentation for details.
    !                        Set by bacoli95_init and bacoli95.
    !
    ! variable problem parameters (output)
    integer           :: idid
    !                        idid is the status / error code which
    !                        indicates if solution information is
    !                        available or the nature of an error.
    !                        It is used independently by bacoli95_init
    !                        and bacoli95/BACOLI to report errors in
    !                        initializations and integration,
    !                        respectively.
    !                        idid = 0 is the initial value set by a
    !                        successful call to bacoli95_init.
    !                        idid > 0 indicates that time integration
    !                        to the current value of t0 was successful.
    !                        -1000 < idid < 0 indicates an error arose
    !                        within the call to BACOLI, and a message
    !                        was printed to stdout (unit=6).
    !                        idid <= -1000 indicates that an error
    !                        occurred during the bacoli95_init call,
    !                        and a message was printed to the unit
    !                        eunit. eunit is a module variable,
    !                        initially 6 (for stdout).
    !                        For more details regarding idid, refer
    !                        to BACOLI (bacoli.f) documentation.
    !                        Set by bacoli95_init and bacoli95.
    real(kind=8), allocatable :: y(:)
    !                        BACOLI approximates solutions using
    !                        B-splines, and y is the vector of
    !                        coefficients which, when combined with
    !                        a B-spline basis, provides solution
    !                        information. To extract solution
    !                        information, use bacoli95_vals or an
    !                        alternative described along with the
    !                        bacoli95_splines type (see below).
    !                        Set by bacoli95.
    ! miscellaneous counters and values (output by bacoli95)
    integer           :: num_remeshings
    !                        Number of spatial remeshings.
    integer           :: num_ini_remeshings
    !                        Number of spatial remeshings on the
    !                        initial step.
    integer           :: num_cold_restarts
    !                        Number of cold restarts of DASSL.
    !                        When BACOLI adapts the spatial mesh,
    !                        the DAE system passed to DASSL is
    !                        different from the previous one and
    !                        a restart of some sort is required.
    !                        BACOLI attempts to do warm restarts
    !                        several times before resorting to
    !                        a cold restart, so these are a sign
    !                        of difficulties in the spatial error
    !                        estimation.
    integer           :: num_accepted_time_steps
    !                        Number of time steps taken by DASSL
    !                        whose spatial error estimate was
    !                        found to be within tolerances.
    integer           :: min_len_ipar
    !                        The number of locations within ipar
    !                        currently being used.
    integer           :: min_len_rpar
    !                        The number of locations within rpar
    !                        currently being used.
    integer           :: prev_bdf_order
    !                        The order of the BDF method employed
    !                        by DASSL on the most recent time step.
    real(kind=8)      :: prev_time_step_size
    !                        The size of the most recent time step
    !                        taken by DASSL.

contains

    ! Set logical flags for bacoli95_init arguments.
    ! Explanations of passed parameters are given with sol parameters.
    ! Default values used for "unpassed" optional parameters:
    ! -- tstart    = -1
    ! -- atol      = [-1] 
    ! -- is_atol   = False
    ! -- rtol      = [-1]
    ! -- is_rtol   = False
    ! -- kcol      = -1
    ! -- nint_max  = -1
    ! -- estimator = -1
    ! -- dirichlet = -1
    ! -- maxord    = -1
    ! -- ini_ss    = -1
    ! Default values required since f2py always provides all parameters.
    subroutine bacoli95_init_wrap(npde,                      &
                                  x, x_size,                 &
                                  tstart,                    &
                                  atol, atol_size, is_atol,  &
                                  rtol, rtol_size, is_rtol,  &
                                  kcol, nint_max, estimator, &
                                  dirichlet, maxord, ini_ss)
        ! Parameters for f2py allocation of x, atol and rtol.
        integer                              :: x_size, atol_size, rtol_size
        ! Parameters to indicate the presence of atol and rtol.
        logical, intent(in)                  :: is_atol, is_rtol
    !f2py integer intent (hide), depend (x) :: x_size = len(x)
        real(kind=8), intent(in)             :: x(x_size)
        real(kind=8), intent(in)             :: tstart
    !f2py integer intent (hide), depend (atol) :: atol_size = len(atol)
        real(kind=8), intent(in)             :: atol(atol_size)
    !f2py integer intent (hide), depend (rtol) :: rtol_size = len(rtol)
        real(kind=8), intent(in)             :: rtol(rtol_size)
        integer,  intent(in)                 :: kcol
        integer,  intent(in)                 :: nint_max
        integer,  intent(in)                 :: estimator
        integer,  intent(in)                 :: dirichlet
        integer,  intent(in)                 :: maxord
        integer,  intent(in)                 :: npde
        real(kind=8), intent(in)             :: ini_ss

        ! Initialize the flags structured type.
        if (tstart /= -1.0d0) then
            flags%is_tstart = .true.
        end if
        if (is_atol) then
            flags%is_atol = .true.
        end if
        if (is_rtol) then
            flags%is_rtol = .true.
        end if
        if (kcol /= -1) then
            flags%is_kcol = .true.
        end if
        if (nint_max /= -1) then
            flags%is_nint_max = .true.
        end if
        if (estimator /= -1) then
            flags%is_estimator = .true.
        end if
        if (dirichlet /= -1) then
            flags%is_dirichlet = .true.
        end if
        if (maxord /= -1) then
            flags%is_maxord = .true.
        end if
        if (ini_ss /= -1) then
            flags%is_ini_ss = .true.
        end if
        
        call bacoli95_init(sol, flags, npde, x, tstart, atol, rtol, kcol, &
                nint_max, estimator, dirichlet, maxord, ini_ss)
    end subroutine bacoli95_init_wrap

    ! Wrapping subroutine to call bacoli95 subroutine.
    ! ------------------------------------------------
    ! Default values for "unpassed" optional parameters:
    ! -- tstop   = -1
    ! -- nsteps = -1
    ! -- dirivf, difbxa, difbxb = Dummy callback functions.
    ! Values for flags if dummy functions have not been sent:
    ! -- is_derivf = False
    ! -- is_difbxa = False
    ! -- is_difbxb = False
    ! Default values required since f2py always provides all parameters.
    subroutine bacoli95_wrap(output_time, f, bndxa, bndxb, uinit, derivf,     &
                             is_derivf, difbxa, is_difbxa, difbxb, is_difbxb, &
                             tstop, nsteps)

        real(kind=8), intent(in)              :: tstop
        real(kind=8), intent(in)              :: output_time
        integer, intent(in)                   :: nsteps

        external                             f
        external                             bndxa
        external                             bndxb
        external                             uinit
        external                             derivf
        external                             difbxa
        external                             difbxb
        optional  derivf, difbxa, difbxb
        ! Booleans to indicate if optional callback functions have been 
        ! provided.
        logical :: is_derivf, is_difbxa, is_difbxb

        if (is_derivf) then
            flags%is_derivf = .true.
        end if
        if (is_difbxa) then
            flags%is_difbxa = .true.
        end if
        if (is_difbxb) then
            flags%is_difbxb = .true.
        end if
        if (nsteps /= -1) then
            flags%is_nsteps = .true.
        end if
        if (tstop /= -1) then 
            flags%is_tstop = .true.
        end if
        
        call bacoli95(sol, flags, output_time, f, bndxa, bndxb, uinit, &
                derivf, difbxa, difbxb, tstop, nsteps)
    end subroutine bacoli95_wrap

    
    ! Wrapping subroutine to make call to bacoli95_vals. Used in case
    ! where output array is of rank 2.
    ! Default value for "unpassed" optional parameter:
    ! -- nderiv - -1
    ! Cases for different uout ranks:
    !     rank = 2: nderiv has not been provided.
    !     rank = 3: nderiv has been provided.
    subroutine bacoli95_vals_rank2_wrap(uout, xout, xout_size, nderiv)

        real(kind=8), intent(out)  :: uout(:,:)
        ! Parameter to allow f2py allocation of xout.
        integer                                 :: xout_size
    !f2py integer intent(hide), depend (xout) :: xout_size = len(xout)
        real(kind=8), intent(in)                :: xout(xout_size)
        ! Optional parameter flags.
        integer,  intent(in)                    :: nderiv

        ! TODO, remove these loop variables, they are for testing.
        integer                                 :: i, j

        if (nderiv /= -1) then
            flags%is_nderiv = .true.
        end if

        call bacoli95_vals(sol, flags, xout, uout, nderiv)
    end subroutine bacoli95_vals_rank2_wrap

    ! Wrapping subroutine to make call to bacoli95_vals. Used in case
    ! where output array is of rank 2.
    ! Default value for "unpassed" optional parameter:
    ! -- nderiv - -1
    ! Cases for different uout ranks:
    !     rank = 2: nderiv has not been provided.
    !     rank = 3: nderiv has been provided.
    subroutine bacoli95_vals_rank3_wrap(uout, xout, xout_size, nderiv)

        real(kind=8), intent(out)   :: uout(:,:,:)
        ! Parameter to allow f2py allocation of xout.
        integer                                  :: xout_size
    !f2py integer intent(hide), depend (xout) :: xout_size = len(xout)
        real(kind=8), intent(in)                 :: xout(xout_size)
        ! Optional parameter flags.
        integer,  intent(in)                     :: nderiv

        if (nderiv /= -1) then
            flags%is_nderiv = .true.
        end if

        call bacoli95_vals(sol, flags, xout, uout, nderiv)
    end subroutine bacoli95_vals_rank3_wrap 

    ! Subroutine to unpack all parts of the solution type. All variables
    ! from this structured type are unpacked into analogous local variables
    ! to allow easy extraction by a python wrapper.
    ! ---------------------------------------------------------------------
    ! Future work must be done on this subroutine to allow extraction of
    ! the final spatial mesh, x and solution coeffiencients y.
    subroutine unpack_sol_parameters()

        ! Allocate output arrays.
        allocate(rtol(size(sol%rtol)), atol(size(sol%rtol)), x(size(sol%x)),  &
                 y(size(sol%y)))

        rtol = sol%rtol
        atol = sol%atol
        x = sol%x
        y = sol%y
        npde = sol%npde
        nint_max = sol%nint_max
        kcol = sol%kcol
        estimator = sol%estimator
        maxord = sol%maxord
        ini_ss = sol%ini_ss
        t0 = sol%t0
        nint = sol%nint
        mflag = sol%mflag
        idid = sol%idid
        num_remeshings = sol%num_remeshings
        num_ini_remeshings = sol%num_ini_remeshings
        num_cold_restarts = sol%num_cold_restarts
        num_accepted_time_steps = sol%num_accepted_time_steps
        prev_bdf_order = sol%prev_bdf_order
        prev_time_step_size = sol%prev_time_step_size
    end subroutine unpack_sol_parameters

end module bacoli95_interface