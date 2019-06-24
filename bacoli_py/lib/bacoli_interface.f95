module type_interface
    type, public :: interface_sol 
        integer                   :: npde

        integer                   :: kcol

        integer                   :: t_est

        integer                   :: s_est

        double precision, pointer :: atol(:)

        double precision, pointer :: rtol(:)

        double precision          :: t0

        double precision, pointer :: x(:)

        integer                   :: nint

        integer                   :: nint_max

        integer                   :: mflag(12)

        double precision, pointer :: y(:)

        double precision, pointer :: rpar(:)

        integer,  pointer         :: ipar(:)

        double complex,  pointer  :: cpar(:)

        double precision, pointer :: xspan(:)

        integer                   :: idid

        double precision, pointer :: output_mem_u(:)

        double precision, pointer :: output_mem_ux(:)

        double precision, pointer :: output_mem_utemp(:,:,:)

        logical                   :: gen_mesh

    end type
end module type_interface

module bacoli_interface
    use type_interface

    implicit none
    ! Interface with wraps the BACOLI and BACOLRI PDE solvers for use
    ! by a Python caller.

    ! Sol object to encapsulate BACOL(R)I's solution information,
    ! preserving it through multiple calls.
    type(interface_sol) :: sol

    private

    public :: initialize
    ! Initializes the sol object with all input parameters needed for calls to
    ! BACOL(R)I.
    !
    ! subroutine initialize(npde, nint_max, kcol, s_est, t_est, maxord, &
    !     is_maxord, atol, atol_size, rtol, rtol_size, t0, x, x_size,   &
    !     xspan, xspan_size, dirichlet, ini_ss, is_ini_ss, tstop,       &
    !     is_tstop)
    !
    ! Note: All arguments suffixed with '_size' need not be passed explicitly
    !       in Python call to the subroutine.
    ! Args:
    ! npde:      The number of PDE's to be solved.
    ! nint_max:  The maximum number of spatial subintervals allowed in this
    !            computation.
    ! kcol:      The number of collaction points used per spatial subinterval.
    ! s_est:     The spatial error control scheme to be used in this computation.
    !            s_est = 0: Use LOI/LE Local Extrapolation
    !                       Error control mode.
    !            s_est = 1: Use SCI/ST Superconvergent
    !                       Interpolation Error Control Mode.
    ! t_est:     The temporal error control scheme to be used in this computation.
    !            t_est = 0: Use variable order BDF method with
    !                       adaptive error control provided
    !                       by DASSL.
    !            s_est = 1: Use 5th order Runge-Kutta method
    !                       with adaptive error control provided
    !                       by RADAU5.
    ! maxord:    Maximum order of the BDF method to be used in time integration if using
    !            DASSL for time stepping.
    ! is_maxord: Flag indicating if maxord is specified.
    ! atol:      Absolute error tolerances to be used in this compuation.
    ! rtol:      Relative error tolerances to be used in this compuation.
    ! t0:        The time at which to begin the compuation.
    ! x:         The initial spatial mesh.
    ! xspan:     The x-values at which the solution will be output.
    ! dirichlet: Indicates if both boundary conditions are dirichlet.
    ! ini_ss:    Initial step size to be used in time integration.
    ! is_ini_ss: Indicates if ini_ss is to considered. Else is chosen by
    !            time stepping software.
    ! t_stop:    Indicates the absolute end of the temporal domain.
    ! is_t_stop: Indicates if t_stop is provided.

    public :: solve
    ! Subroutine which solves the initialized PDE system for a given time
    ! tout.
    ! subroutine solve(tout, f, bndxa, bndxb, uinit, derivf, is_derivf,  &
    !     difbxa, is_difbxa, difbxb, is_difbxb)
    !
    ! Args:
    ! tout:      The time that the solution should be computed for.
    ! f:         The PDE system to be solved.
    ! bndxa:     The left boundary conditions.
    ! bndxb:     The right boundary conditions.
    ! uinit:     The initial conditions.
    ! derivf:    Analytic jacobian of the PDE system.
    ! is_derivf: Indicates if derivf provided.
    ! difbxa:    Analytic jacobian of the left boundary condition.
    ! is_difbxa: Indicates if difbxa provided.
    ! difbxb:    Analytic jacobian of the right boundary condition.
    ! is_difbxa: Indicates if difbxb provided.

    public :: vals
    ! Subroutine to extract solution information from BACOL(R)I.
    ! subroutine vals(u)
    !
    ! Args:
    ! u:          Array which will hold solution information.
    ! ux:         Array which will hold solution derivative information.
    ! npde:       Number of PDE's in this system.
    ! xspan_size: Number of output points on the x-axis.
    ! nderiv:     0 if no derivative informationr requested. Otherwise 1.

    ! public :: deriv_vals
    ! ! Subroutine to extract solution derivative information from
    ! ! BACOL(R)I. subroutine vals(u)
    ! !
    ! ! Args:
    ! ! ux:         Array which will hold derivative information.
    ! ! npde:       Number of PDE's in this system.
    ! ! xspan_size: Number of output points on the x-axis.


    public :: sol_teardown
    ! Subroutine to deallocate memory and free pointers which were
    ! set in call to initialize. This should be done following the
    ! conclusion of all computations.
    ! subroutine sol_teardown()

    contains

        subroutine initialize(npde, nint_max, kcol, s_est, t_est, maxord,     &
            is_maxord, atol, atol_size, rtol, rtol_size, t0, x, x_size,       &
            xspan, xspan_size, dirichlet, ini_ss, is_ini_ss, tstop, is_tstop, &
            idid)
            implicit none

            ! f2py provided arguments
            integer, intent(in)  :: atol_size, rtol_size, x_size, xspan_size
            integer, intent(in)  :: npde, nint_max, kcol, s_est, t_est, maxord, &
                dirichlet
            integer, intent(out) :: idid
            logical, intent(in)  :: is_maxord, is_ini_ss, is_tstop
            double precision, intent(in) :: t0, ini_ss, tstop
        !f2py integer intent (hide), depend (atol) :: atol_size = len(atol)
            double precision, intent(in) :: atol(atol_size)
        !f2py integer intent (hide), depend (rtol) :: rtol_size = len(atol)
            double precision, intent(in) :: rtol(rtol_size)
        !f2py integer intent (hide), depend (x) :: x_size = len(x)
            double precision, intent(in) :: x(x_size)
        !f2pt integer intent (hide), depend (xspan) :: xspan_size = len(xspan)
            double precision, intent(in) :: xspan(xspan_size)

            ! Initialize constant values in sol structured type
            sol%idid     = 0
            sol%npde     = npde
            sol%kcol     = kcol
            sol%s_est    = s_est
            sol%t_est    = t_est
            sol%nint     = size(x) - 1
            sol%nint_max = nint_max
            sol%mflag    = 0
            sol%t0       = t0

            nullify(sol%x, sol%y, sol%atol, sol%rtol, sol%rpar, sol%ipar, &
                sol%cpar)

            call allocate_arrays(x, xspan, rtol, atol, npde, kcol, nint_max)
            ! If there was an error in array allocation, return error flag to
            ! caller.
            if (sol%idid < 0) then
                goto 101
            end if

            ! Check for scalar error tolerances and set appropriate flags
            if (any(sol%atol /= sol%atol(1))) then
                sol%mflag(2) = 1
            end if
            if (any(sol%rtol /= sol%rtol(1))) then
                sol%mflag(2) = 1
            end if

            ! Set flag for which spatial error estimation scheme to use
            if (t_est == 0) then
                sol%mflag(8) = s_est
            else
                sol%mflag(5) = s_est
            end if

            ! If BACOLI being used, handle case where maximum BDF order is specified.
            if (t_est == 0) then
                if (is_maxord) then
                    sol%mflag(7) = 1
                    sol%ipar(15) = maxord
                end if
            end if

            ! Check for and handle case where initial stepsize if speciifed by user.
            if (is_ini_ss) then
                if (t_est == 0) then
                    sol%mflag(6) = 1
                    sol%rpar(2)  = ini_ss
                else
                    sol%mflag(4) = 1
                    sol%rpar(2)  = ini_ss
                end if
            end if

            ! Check for and handle case where tstop is provided by user.
            if (t_est == 0) then
                if (is_tstop) then
                    sol%mflag(3) = 1
                    sol%rpar(1)  = tstop
                end if
            end if

            ! Check if both boundary conditions are dirichet.
            if (dirichlet /= 0) then
                if (t_est == 0) then
                    sol%mflag(5) = 1
                else
                    sol%mflag(3) = 1
                end if
            end if
        101 continue

            idid = sol%idid
        end subroutine initialize
!-----------------------------------------------------------------------

        subroutine solve(tout, f, fvec, bndxa, bndxb, uinit, derivf, is_derivf, &
            difbxa, is_difbxa, difbxb, is_difbxb, idid)
            implicit none

            double precision, intent(in) :: tout        
            integer, intent(out)         :: idid
            logical, intent(in)          :: is_derivf, is_difbxa, is_difbxb
            logical                      :: use_fd_both, use_fd_bconds, use_fd_pde
            external                        f
            ! BACOLI --> BACOLIVEC
            external                        fvec
            external                        bndxa
            external                        bndxb
            external                        uinit
            external                        derivf
            external                        difbxa
            external                        difbxb
            external                        bacoli
            external                        bacolri
            external                        meshgen
            double precision              :: xa, xb
            double precision, allocatable :: wm(:), xi(:)
            integer                       :: ier, i
 
            ! Set flags about whether or not to use finite differences.
            use_fd_both = (.not. is_derivf) .and. ((.not. is_difbxa) &
                .or. (.not. is_difbxb))
            use_fd_bconds = is_derivf .and. ((.not. is_difbxa)     &
                .or. (.not. is_difbxb))
            use_fd_pde = (.not. is_derivf) .and. is_difbxa .and.     &
                is_difbxb

            if (sol%mflag(1) == 0 .and. sol%gen_mesh) then
                write(*,*) 'in sol using generated mesh'
                ! Extract endpoints
                xa = sol%x(1)
                xb = sol%x(2)

                ! User-supplied mesh case
                allocate(wm((sol%nint+1)+3*sol%npde), stat=ier)

                allocate(xi(sol%nint+1), stat=ier)

                do i = 1, sol%nint+1
                    xi(i) = xa + (xb - xa)*dble(i-1)/dble(sol%nint)
                end do

                call meshgen(xa, xb, sol%x, xi, sol%nint, uinit, sol%npde, wm, ier)

                ! write(*,*) 'nint ', sol%nint
                ! write(*,*) 'xa ', xa
                ! write(*,*) 'xb ', xb

                ! Free memory
                deallocate(wm)
                deallocate(xi)
            end if


            ! write(*,*) 'continuation?'
            ! write(*,*) sol%t0

            if (sol%t_est == 0) then
                if (use_fd_both) then
                    sol%mflag(9) = 0
                    call bacoli(sol%t0, tout, sol%atol, sol%rtol, sol%npde,   &
                        sol%kcol, sol%nint_max, sol%nint, sol%x, sol%mflag,   &
                        sol%rpar, size(sol%rpar), sol%ipar, size(sol%ipar),   &
                        sol%y, sol%idid, fvec, dummy_derivf, bndxa, dummy_difbx, &
                        bndxb, dummy_difbx, uinit, uinitvec)
                else if (use_fd_bconds) then
                    sol%mflag(9) = 1
                    call bacoli(sol%t0, tout, sol%atol, sol%rtol, sol%npde,   &
                        sol%kcol, sol%nint_max, sol%nint, sol%x, sol%mflag,   &
                        sol%rpar, size(sol%rpar), sol%ipar, size(sol%ipar),   &
                        sol%y, sol%idid, fvec, derivf, bndxa, dummy_difbx,       &
                        bndxb, dummy_difbx, uinit, uinitvec)
                else if (use_fd_pde) then
                    sol%mflag(9) = 2
                    call bacoli(sol%t0, tout, sol%atol, sol%rtol, sol%npde,   &
                        sol%kcol, sol%nint_max, sol%nint, sol%x, sol%mflag,   &
                        sol%rpar, size(sol%rpar), sol%ipar, size(sol%ipar),   &
                        sol%y, sol%idid, fvec, dummy_derivf, bndxa, difbxa,      &
                        bndxb, difbxb, uinit, uinitvec)
                else
                    sol%mflag(9) = 3
                    call bacoli(sol%t0, tout, sol%atol, sol%rtol, sol%npde,   &
                        sol%kcol, sol%nint_max, sol%nint, sol%x, sol%mflag,   &
                        sol%rpar, size(sol%rpar), sol%ipar, size(sol%ipar),   &
                        sol%y, sol%idid, fvec, derivf, bndxa, difbxa, bndxb,     &
                        difbxb, uinit, uinitvec)
                end if
            else
                if (use_fd_both) then
                    sol%mflag(6) = 0
                    call bacolri(sol%t0, tout, sol%atol, sol%rtol, sol%npde,  &
                        sol%kcol, sol%nint_max, sol%nint, sol%x, sol%mflag,   &
                        sol%rpar, size(sol%rpar), sol%ipar, size(sol%ipar),   &
                        sol%cpar, size(sol%cpar), sol%y, sol%idid, f, fvec,   &
                        dummy_derivf, bndxa, dummy_difbx, bndxb, dummy_difbx, &
                        uinit, uinitvec)
                else if (use_fd_bconds) then
                    sol%mflag(6) = 1
                    call bacolri(sol%t0, tout, sol%atol, sol%rtol, sol%npde,  &
                        sol%kcol, sol%nint_max, sol%nint, sol%x, sol%mflag,   &
                        sol%rpar, size(sol%rpar), sol%ipar, size(sol%ipar),   &
                        sol%cpar, size(sol%cpar), sol%y, sol%idid, f, fvec,   &
                        derivf, bndxa, dummy_difbx, bndxb, dummy_difbx,       &
                        uinit, uinitvec)
                else if (use_fd_pde) then
                    sol%mflag(6) = 2
                    call bacolri(sol%t0, tout, sol%atol, sol%rtol, sol%npde,  &
                        sol%kcol, sol%nint_max, sol%nint, sol%x, sol%mflag,   &
                        sol%rpar, size(sol%rpar), sol%ipar, size(sol%ipar),   &
                        sol%cpar, size(sol%cpar), sol%y, sol%idid, f, fvec,   &
                        dummy_derivf, bndxa, difbxa, bndxb, difbxb, uinit, uinitvec)
                else
                    sol%mflag(6) = 3
                    call bacolri(sol%t0, tout, sol%atol, sol%rtol, sol%npde,  &
                        sol%kcol, sol%nint_max, sol%nint, sol%x, sol%mflag,   &
                        sol%rpar, size(sol%rpar), sol%ipar, size(sol%ipar),   &
                        sol%cpar, size(sol%cpar), sol%y, sol%idid, f, fvec,   &
                        derivf, bndxa, difbxa, bndxb, difbxb, uinit, uinitvec)
                end if
            end if

            ! If not first call, set problem continuation flag.
            if (sol%mflag(1) /= 1) then
                sol%mflag(1) = 1
            else
                sol%idid = 1
            end if

            idid = sol%idid
        end subroutine solve

!-----------------------------------------------------------------------

        subroutine vals(u, ux, npde, xspan_size, nderiv)
            implicit none
            integer, intent(in)           :: npde, xspan_size
            integer, intent(in)           :: nderiv
            double precision, intent(out) :: u(npde, xspan_size)
            double precision, intent(out) :: ux(npde, xspan_size)
            integer :: i, j

            external                         values
            external                         valuesri


            ! Get values at given points in xspan
            if (sol%t_est == 0) then
                call values(sol%kcol, sol%xspan, sol%nint, sol%x, sol%npde, &
                    size(sol%xspan), nderiv, sol%output_mem_utemp, sol%y, sol%output_mem_ux)
            else
                call valuesri(sol%kcol, sol%xspan, sol%nint, sol%x, sol%npde, &
                    size(sol%xspan), nderiv, sol%output_mem_utemp, sol%y, sol%output_mem_ux)
            end if

            ! Extract derivative values.
            do i = 1, npde
                do j = 1, xspan_size
                    u(i,j)  = sol%output_mem_utemp(i,j,1)
                    ux(i,j) = sol%output_mem_utemp(i,j,2)
                end do
            end do
        end subroutine vals

!-----------------------------------------------------------------------

        subroutine sol_teardown()
            implicit none
            if (associated(sol%x))    then
                deallocate(sol%x) ;    nullify(sol%x)
            end if
            if (associated(sol%y))    then
                deallocate(sol%y) ;    nullify(sol%y)
            end if
            if (associated(sol%atol)) then
                deallocate(sol%atol) ; nullify(sol%atol)
            end if
            if (associated(sol%rtol)) then
                deallocate(sol%rtol) ; nullify(sol%rtol)
            end if
            if (associated(sol%rpar)) then
                deallocate(sol%rpar) ; nullify(sol%rpar)
            end if
            if (associated(sol%ipar)) then
                deallocate(sol%ipar) ; nullify(sol%ipar)
            end if
            if (associated(sol%xspan)) then
                deallocate(sol%xspan) ; nullify(sol%xspan)
            end if

            if (associated(sol%output_mem_u)) then
                deallocate(sol%output_mem_u) ; nullify(sol%output_mem_u)
            end if
            if (associated(sol%output_mem_ux)) then
                deallocate(sol%output_mem_ux) ; nullify(sol%output_mem_ux)
            end if
            if (associated(sol%output_mem_utemp)) then
                deallocate(sol%output_mem_utemp) ; nullify(sol%output_mem_utemp)
            end if

            if (sol%t_est == 1) then
                if (associated(sol%cpar)) then
                    deallocate(sol%cpar) ; nullify(sol%cpar)
                end if
            end if
        end subroutine

!-----------------------------------------------------------------------

        subroutine allocate_arrays(x, xspan, rtol, atol, npde, kcol, nint_max)
            implicit none
            integer                      :: cpar, rpar, ipar, ly, ier
            double precision, intent(in) :: x(:), xspan(:), rtol(:), atol(:)
            integer, intent(in)          :: npde, kcol, nint_max
            integer                      :: lrp, lcp, lip
            integer                      :: lenwrk_u, lenwrk_ux

            if (sol%t_est == 0) then
                ! Set work array sizes for BACOLI.
                lrp = 113 + 59*npde + 27*nint_max + 13*npde*npde       &
                   + 9*kcol + 24*kcol*nint_max + 6*nint_max*kcol*kcol &
                   + 27*npde*nint_max*kcol + 7*nint_max*npde          &
                   + 2*npde*npde*nint_max*kcol*kcol                   &
                   + 4*npde*npde*kcol*nint_max
                ! BACOLI --> BACOLIVEC
                ! lrp = 113 + 65*npde + 27*nint_max + 13*npde*npde       &
                !     + 9*kcol + 24*kcol*nint_max + 6*nint_max*kcol*kcol &
                !     + 33*npde*nint_max*kcol + 7*nint_max*npde          &
                !     + 2*npde*npde*nint_max*kcol*kcol                   &
                !     + 4*npde*npde*kcol*nint_max


                ly = npde*(kcol*nint_max+2)

                lip = 115 + ly

                ! The LOI/LE scheme requires slightly less storage.
                if (sol%s_est == 0) then
                    lrp = lrp - 15*nint_max + 3*kcol - 8*kcol*nint_max &
                        + kcol*kcol - nint_max*kcol*kcol               &
                        - 3*nint_max*npde + 2
                end if

                ! BACOLI --> BACOLIVEC
                lrp = lrp + 6*npde*(kcol*nint_max+1)

                ! Allocate the work arrays.
                allocate(sol%ipar(lip), sol%rpar(lrp), sol%y(ly), stat=ier)
                sol%ipar = 0
                sol%rpar = 0
                sol%y    = 0
            else
                ! Set work array sizes for BACOLRI.
                ly = npde*(kcol*nint_max+2)

                lcp = npde*(2+kcol*nint_max)+npde*npde &
                    *(4+nint_max*kcol*(kcol+2))

                lip = 102+3*npde*(nint_max*kcol+4)

                lrp = 55+8*npde*npde*nint_max*kcol+39*npde &
                    +9*kcol+28*nint_max+24*nint_max*kcol   &
                    +4*npde*npde*nint_max*kcol*kcol        &
                    +17*npde*nint_max*kcol+21*npde*npde    &
                    +6*nint_max*kcol*kcol                  &
                    +7*npde*nint_max

                ! The LOI/LE scheme requires slightly less storage.
                if (sol%s_est == 0) then
                    lrp = lrp+2+3*kcol+kcol*kcol-18*nint_max &
                        -nint_max*kcol*kcol                  &
                        -9*nint_max*kcol-3*npde*nint_max          
                end if

                ! BACOLRI --> BACOLRIVEC
                lrp = lrp + 6*npde*(kcol*nint_max+1)


                ! Allocate the work arrays.
                allocate(sol%ipar(lip), sol%rpar(lrp), sol%y(ly),  &
                    sol%cpar(lcp), stat=ier)
                if (ier /= 0) then
                    sol%idid = -1000 ;  goto 100
                end if

                sol%ipar = 0
                sol%rpar = 0
                sol%y    = 0
                sol%cpar = 0
            end if

            ! Allocate spatial mesh, initialize with initial mesh given by
            ! user.
            allocate(sol%x(nint_max+1), stat=ier)
            if (ier /= 0) then
                sol%idid = -1000 ;  goto 100
            end if
            !sol%x = 0
            sol%x(1:size(x)) = x
            sol%x(size(x)+1:size(x)) = 0

            if (size(x) > 2) then
                ! Using the user mesh in solve()
                write(*,*) 'using user mesh'
                sol%gen_mesh = .false.

                ! nint is based on the user input
                sol%nint = size(x) - 1
            else
                write(*,*) 'using generated mesh'
                ! Will be generating the mesh in solve()
                sol%gen_mesh = .true.

                ! Set number of mesh subintervals as max(10% of nint_max, 10).
                sol%nint = max(ceiling(0.1*sol%nint_max), 10)
            end if



            ! Allocate error tolerance arrays
            allocate(sol%atol(npde), stat=ier)
            if (ier /= 0) then
                sol%idid = -1000;    goto 100
            end if
            sol%atol = atol(1)

            allocate(sol%rtol(npde), stat=ier)
            if (ier /= 0) then
                sol%idid = -1000;    goto 100
            end if
            sol%rtol = rtol(1)

            ! Allocate the array of x values at which solution will be output
            allocate(sol%xspan(size(xspan)), stat=ier)
            if (ier /= 0) then
                sol%idid = -1000;    goto 100
            end if
            sol%xspan = xspan

            ! Allocate memory for the values routine
            lenwrk_u = (kcol+2) + kcol*(nint_max+1) + 4
            allocate(sol%output_mem_u(lenwrk_u), stat=ier)
            if (ier /= 0) then
                sol%idid = -1000;    goto 100
            end if

            lenwrk_ux = (kcol+2)*2 + kcol*(nint_max+1) + 4
            allocate(sol%output_mem_ux(lenwrk_ux), stat=ier)
            if (ier /= 0) then
                sol%idid = -1000;    goto 100
            end if

            ! Allocate the temp array holding the solution and derivative values.
            allocate(sol%output_mem_utemp(sol%npde, size(xspan), 2), stat=ier)
            if (ier /= 0) then
                sol%idid = -1000;    goto 100
            end if


        100 continue
            ! If error, simply return to caller which will check which error
            ! flag was set.
            return
        end subroutine

!-----------------------------------------------------------------------

        subroutine dummy_derivf(t, x, u, ux, uxx, dfdu, dfdux, &
              dfduxx, npde)
           implicit none
           ! Dummy - used if user does not provide derivf.
           integer,  intent(in)  :: npde
           double precision, intent(in)  :: t, x, u(npde)
           double precision, intent(in)  :: ux(npde), uxx(npde)
           double precision, intent(out) :: dfdu(npde,npde)
           double precision, intent(out) :: dfdux(npde,npde)
           double precision, intent(out) :: dfduxx(npde,npde)

           dfdu = 0 ; dfdux = 0 ; dfduxx = 0
        end subroutine

!-----------------------------------------------------------------------
        subroutine uinitvec(x, u, vnpts, npde)
            integer :: vnpts, npde
            integer :: x(vnpts)
            double precision :: u(vnpts*npde)
        end subroutine
!-----------------------------------------------------------------------

        subroutine dummy_difbx(t, u, ux, dbdu, dbdux, dbdt, npde)
           implicit none
           ! Dummy - used if user does not provide difbxa, difbxb.
           integer,  intent(in)  :: npde
           double precision, intent(in)  :: t, u(npde), ux(npde)
           double precision, intent(out) :: dbdu(npde,npde)
           double precision, intent(out) :: dbdux(npde,npde)
           double precision, intent(out) :: dbdt(npde)

           dbdu = 0 ; dbdux = 0 ; dbdt = 0
        end subroutine

end module bacoli_interface