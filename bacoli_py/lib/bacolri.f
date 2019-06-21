      subroutine bacolri(t0, tout, atol, rtol, npde, kcol, nintmx, nint,
     &                   x, mflag, rpar, lrp, ipar, lip, cpar, lcp, y,
     &                   idid, f, fvec, derivf, bndxa, difbxa, bndxb,
     &                   difbxb, uinit, uinitvec)

c-----------------------------------------------------------------------
c Purpose:
c       The purpose of BACOLRI is to solve NPDE dimensional systems of
c       second order parabolic partial differential equations (PDEs)
c       in one space variable of the form:
c
c            dU
c            -- (t,x) = f ( t, x, U(t,x), U_x (t,x), U_{xx} (t,x) ) ,
c            dt
c
c       where xa < x < xb and t > t0, with initial conditions at
c       time t = t0 given by:
c
c                          u(t0,x) = u_0(x),
c
c       for x_a <= x <= x_b, subject to separated boundary conditions
c       given by:
c
c                   b_{xa} ( t, U(t,x_a), U_x (t,x_a) ) = 0,
c
c                   b_{xb} ( t, U(t,x_b), U_x (t,x_b) ) = 0,
c
c       for t > t0 and x = x_a, x = x_b, respectively.
c
c       Guide to the above notation:
c          dU
c          -- (t,x)     denotes the first partial derivative of U(t,x)
c          dt           with respect to the time variable t.
c
c          U_x (t,x)    denotes the first partial derivative of U(t,x)
c                       with respect to space variable x.
c
c          U_{xx} (t,x) denotes the second partial derivative of U(t,x)
c                       with respect to space variable x.
c
c       Also, the above functions are NPDE dimensional vector functions.
c
c       BACOLRI is a method of lines algorithm which uses bspline
c       collocation to discretize the spatial domain [x_a,x_b].
c       The output is a vector of bspline coefficients which
c       can be used to calculate the approximate solution U(t,x) and
c       its spatial derivatives at (tout,x) where x_a <= x <= x_b
c       and t0 < tout.
c
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c Setup of BACOLRI:
c       BACOLRI requires that the user specifies the system of PDEs and
c       the related initial and boundary conditions as also sets
c       input parameters (which define the bspline space and the
c       requested error tolerances) and allocates work storage.
c
c       The calling sequence of BACOLRI is:
c
c       call bacolri(t0, tout, atol, rtol, npde, kcol, nint, nintmx, x,
c    &             mflag, rpar, lrp, ipar, lip, cpar, lcp, y, idid,
c    &             f, derivf, bndxa, difbxa, bndxb, difbxb, uinit)
c
c       which will generate the vector y of bspline coefficients for the
c       approximation of U at t = tout, upon successful completion.
c       Generally, the call to BACOLRI will be followed by a
c       call to VALUES to calculate the solution at a set of points:
c
c       call values(kcol, xsol, nint, x, npde, npts, nderiv,
c    &             usol, y, work)
c
c       The details of the parameters to VALUES are documented within
c       the source code for that routine. The input parameters for
c       BACOLRI are dealt with in detail below, but a quick summary is:
c
c       [t0, tout] is the time domain of the problem.
c       atol is the absolute error tolerance.
c       rtol is the relative error tolerance.
c       npde is the number of components in the PDE system.
c       kcol, nint define the bspline space.
c       nintmx is the maximum number of subintervals allowed.
c       mflag(1:12) is used to control the operation of BACOLRI.
c       x is the spatial mesh.
c       rpar(lrp) is a floating point work array.
c       ipar(lip) is an integer work array.
c       cpar(lcp) is a complex work array.
c       y is the set of coefficients in the bspline approximation to U
c       at t = tout.
c       idid is an exit status flag provided by BACOLRI.
c       The user must check idid to determine what further action needs
c       to be taken.

c       PDE system definition subroutines are to be set up as follows:
c
c       f(t, x, u, ux, uxx, fval, npde)
c       npde is an integer, the rest are double precision.
c       t, x are scalars, u, ux, uxx are vectors of length npde.
c       This subroutine defines the right hand side of the PDE system,
c       and ut = f(t, x, u, ux, uxx) should be returned in fval.
c
c       derivf(t, x, u, ux, uxx, dfdu, dfdux, dfduxx, npde)
c       npde is an integer, the rest are double precision.
c       t, x are scalars, u, ux, uxx are vectors of length npde.
c       dfdu, dfdux, dfduxx are Jacobians of f evaluated at (t, x).
c
c       bndxa(t, u, ux, bval, npde)
c       npde is an integer, the rest are double precision.
c       t is scalar, u, ux, bval are vectors of length npde.
c       This subroutine defines the left boundary conditions,
c       where x = x_a, by b(t, u, ux) = 0.
c       Return the residual b(t, u, ux) in bval.
c
c       difbxa(t, u, ux, dbdu, dbdux, dbdt, npde)
c       npde is an integer, the rest are double precision.
c       t is scalar, u, ux are vectors of length npde.
c       This subroutine defines the differentiated left boundary,
c       where x = x_a.
c       Return the values of the Jacobians of bndxa evaluated at t in
c       dbdu, dbdux and dbdt.
c
c       bndxb and difbxb are the same as bndxa and difbxa, but for the
c       right boundary conditions where x = x_b.
c
c       uinit(x, u, npde)
c       npde is an integer, x is a double precision scalar and
c       u is a double precision vector of length npde.
c       This subroutine defines the initial condition of the PDE system
c       by returning the value of u0(x) = u(t0, x) in u.
c       This initial condition should be C1-continuous.
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
c
        double precision        t0
c       On input, t0 < tout is the initial time. On output, t0 is the
c       current time, t0 <= tout.
c
        double precision        tout
c       tout is the desired final output time. After a successful
c       return from BACOLRI, the time stepping may be resumed by
c       changing tout so that t0 < tout and setting mflag(1) = 1
c       to indicate a continuation of the previous problem.
c
        double precision        atol(npde)
c       atol is the absolute error tolerance requested and
c       is a scalar quantity if mflag(2) = 0.
c
c       If the PDE components vary in importance, then vector error
c       tolerances may be used by setting mflag(2) = 1. In this
c       case, the dimension of atol must be npde. The user will
c       define atol(1), atol(2), ..., atol(npde) appropriately.
c       Note that a change from scalar to vector tolerances (or vice
c       versa) constitutes a new problem, and BACOLRI will have to
c       be reinitialized.
c
        double precision        rtol(npde)
c       rtol is the relative error tolerance request and is a scalar
c       quantity if mflag(2) = 0.
c
c       If the PDE components vary in importance, then vector error
c       tolerances may be used by setting mflag(2) = 1. In this
c       case, the dimension of rtol must be npde. The user will define
c       rtol(1), rtol(2), ..., rtol(npde) appropriately.
c       Note that a change from scalar to vector tolerances (or vice
c       versa) constitutes a new problem, and BACOLRI will have to
c       be reinitialized.
c
        integer                 npde
c       npde is the number of components in the system of PDEs.
c       npde > 0.
c
        integer                 kcol
c       kcol is the number of collocation points to be used in each
c       subinterval. 1 < kcol <= mxkcol.
c
c       The order of the bsplines used will be (kcol+2).
c
        integer                 nint
c       at input, nint is the number of subintervals defined by the
c       spatial mesh x at the initial time t0. at output, nint is
c       the number of subintervals at tout. nint >= 1.
c
        integer                 nintmx
c       the maximum number of subintervals that the user requires.
c
        double precision        x(nintmx+1)
c       x is the spatial mesh which divides the interval [x_a,x_b]
c       as: x_a = x(1) < x(2) < x(3) < ... < x(nint+1) = x_b.
c       At input, x(1:nint+1) stores the mesh points at the initial
c       time t0.  at output, x(1:nint+1) stores the mesh points at tout.
c
        integer                 mflag(12)
c       This vector determines the interaction of BACOLRI with RADAU5
c       and which error estimation scheme BACOLRI will employ.
c
c       How to set mflag(1):
c
c       On the initial call to BACOLRI with a new problem, set
c       mflag(1) = 0, which indicates that BACOLRI and RADAU5 should
c       perform the initialization steps that are required by each code,
c       respectively.
c
c       In order to continue time stepping in the current problem after
c       a successful return from BACOLRI, set mflag(1) = 1,
c       idid = 1, and ensure that t0 < tout.
c
c       How to set mflag(2):
c
c       If scalar absolute and relative error tolerances (atol and rtol)
c       are desired, then set mflag(2) = 0.
c
c       For vector absolute and relative error tolerances, set
c       mflag(2) = 1, define atol(1), ..., atol(npde), and
c       rtol(1), ..., rtol(npde), as described above, ensuring that
c       the dimension of each of atol and rtol is at least npde.
c
c       How to set mflag(3):
c
c       If both boundary conditions are dirichlet, set mflag(3) = 1;
c       else, set mflag(3) = 0.
c
c       How to set mflag(4):
c
c       If the user wants to specify an initial stepsize, set
c       mflag(4) = 1, and define rpar(2) = the initial stepsize;
c       else, set mflag(4) = 0;
c
c       How to set mflag(5):
c
c       If the user wants to run BACOLI with LOI error estimation,
c       set mflag(5) = 0;
c       for the SCI error estimation scheme,
c       set mflag(5) = 1.
c
c       How to set mflag(6) :

c       If the user wants BACOLI to use finite difference
c       approximations in place of the user-provided analytic
c       partial derivative subroutines derivf, difbxa and difbxb,
c           set mflag(6) = 0;
c       If the user implements derivf, but not difbxa or difbxb,
c           set mflag(6) = 1;
c       If the user implements difbxa and difbxb, but not derivf,
c           set mflag(6) = 2;
c       If the user implements all of derivf, difbxa and difbxb,
c           set mflag(6) = 3.
c
c       mflag(7:12): reserved for future use.
c
        integer                 lrp
c       lrp is the size of the rpar storage array and must satisfy:
c       lrp >= 55+8*npde*npde*nintmx*kcol+39*npde
c            + 9*kcol+28*nintmx+24*nintmx*kcol
c            + 4*npde*npde*nintmx*kcol*kcol
c            + 17*npde*nintmx*kcol+21*npde*npde
c            + 6*nintmx*kcol*kcol
c            + 7*npde*nintmx
c     -      + 2+3*kcol+kcol*kcol-18*nintmx
c     -      - nintmx*kcol*kcol
c     -      - 9*nintmx*kcol-3*npde*nintmx
c       BACOLRI --> BACOLRIVEC
c       Modified memory requirements.
c       lrp >= 55+8*npde*npde*nintmx*kcol+39*npde
c            + 9*kcol+28*nintmx+24*nintmx*kcol
c            + 4*npde*npde*nintmx*kcol*kcol
c            + 17*npde*nintmx*kcol+21*npde*npde
c            + 6*nintmx*kcol*kcol
c            + 7*npde*nintmx
c                  change:
c            --------------------
c            + 6*npde*(kcol*nintmx+1)
c            --------------------
c     -      + 2+3*kcol+kcol*kcol-18*nintmx
c     -      - nintmx*kcol*kcol
c     -      - 9*nintmx*kcol-3*npde*nintmx
c       where the last three lines account for the (slightly) reduced
c       memory requirements of the loi error estimator.
c
        integer                 lip
c       lip is the size of the ipar integer work array and must satisfy:
c       lip>=102+3*npde*(nintmx*kcol+4)
c
        integer                 lcp
c       lcp is the size of the cpar complex work array and must satisfy:
c       lcp>= npde*(kcol*nintmx+2)
c           + npde*npde*(4+nintmx*kcol*(kcol+2))

        external                f
        external                fvec
        external                derivf
        external                bndxa
        external                difbxa
        external                bndxb
        external                difbxb
        external                uinit
        external                uinitvec
c                               See "Setup of BACOLI" above.

c
c       Work Storage:
        double precision        rpar(lrp)
c       rpar is a floating point work array of size lrp.
c
        integer                 ipar(lip)
c       ipar is an integer work array of size lip.
c
        double complex          cpar(lcp)
c       cpar is a complex work array of size lcp.
c
c       Output:
        double precision        y(npde*(kcol*nintmx+2))
c       On successful return from BACOLRI, y(1:npde*(kcol*nint+2)) is
c       the vector of bspline coefficients at the current time.
c
        integer                 idid
c       idid is the BACOLRI exit status flag which is based on the exit
c       status from RADAU5 plus some additional status codes based on
c       error checking performed by BACOLRI on initialization. Positive
c       values of idid indicate a successful return. Negative values of
c       idid indicate an error which may or may not be fatal. The exact
c       descriptions of idid return values will be discussed below.
c
C                   IDID= 1  COMPUTATION SUCCESSFUL,
C                   IDID=-1  INPUT IS NOT CONSISTENT,
C                   IDID=-2  LARGER NMAX IS NEEDED,
C                   IDID=-3  STEP SIZE BECOMES TOO SMALL,
C                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR.
C
c-----------------------------------------------------------------------

c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        integer                 mxkcol
        parameter              (mxkcol = 10)
c                               mxkcol is the maximum number of
c                               collocation points per subinterval.
c
        integer                 maxrsh
        parameter              (maxrsh = 40)
c                               maxrsh is the maximum number of
c                               remesh times at one time step,
c                               i.e., ipar(icount) must be less than or
c                               equal to maxrsh
c
        integer                 mxrhin
        parameter              (mxrhin = 500)
c                               mxrhin is the maximum number of
c                               remesh times at the initial step,
c                               i.e., icount must less than or
c                               equal to mxrhin
c
        double precision        point1
        parameter              (point1 = 0.1D0)
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c Local variables:
c
        integer                 kerr
c                               kerr is what kcol would be for the
c                               solution whose error is being estimated.
c                               This helps explain how the LOI error
c                               estimate is working.
c
        integer                 neq
c                               neq=npde*ncpts is the number of
c                               bspline coefficients (or DAEs) when
c                               using radau_{kcol}.
c
        integer                 leniw
c                               leniw = 20 + neq is the length of the
c                               integer work array required by RADAU5.
c
        integer                 lenpd
c                               lenpd is the size of the Almost Block
c                               Diagonal (ABD) Jacobian required by
c                               radau_{kcol}.
c                               lenpd=npde*npde*(2*nconti
c                                      +kcol*(kcol+nconti)*nint)
c
        integer                 lenrw
c                               lenrw = 20+12*neq+3*lenpd
c                               is the total size of the floating point
c                               work array required by radau_{kcol}.
c
        integer                 lenin
c                               lenin is the size of the floating
c                               point work array used by INIY
c                               when using radau_{kcol}.
c                               lenin>=lenpd+2*neq+npde*4+2*npde*npde
c                               BACOLRI --> BACOLRIVEC
c                               lenin>=lenpd+2*neq+npde*4+2*npde*npde
c                                      + 3*npde*nint*kcol
c
        integer                 lenri
c                               lenri is the size of the floating
c                               point work array used by REINIT when
c                               using radau_{kcol}.
c
        integer                 lenrj
c                               lenrj is the size of the floating
c                               point work array used by RES and JAC.
c                               lenrj>=6*npde+5*npde*npde.
c
        integer                 lencof
c                               Length of the work storage array for
c                               Hermite-Birkhoff coefficients and
c                               values of B-spline basis functions used
c                               inside errest's interpolation
c                               subroutines. They only need to be
c                               recalculated when the mesh changes.
c                               The SCI estimate requires more storage
c                               than the LOI estimate here.
c
        integer                 lenerr
c                               lenerr is the size of the floating point
c                               work array used by ERREST.
c                               lenerr>=2*npde*necpts+npde*nint.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the total
c                               number of collocation points when using
c                               radau_{kcol}.
c
        integer                 necpts
c                               necpts=(kcol+3)*nint is the total number
c                               of collocation points used for
c                               error estimate.
c
        integer                 icflag
c                               This is the status flag from the almost
c                               block diagnonal factorization routine,
c                               CRDCMP.
c                               icflag =  0, indicates non-singularity.
c                               icflag = -1, indicates singularity.
c                               icflag =  1, indicates invalid input.
c
        double precision        torign
c                               torign is the initial time, i.e. = t0
c                               at the beginning.
c
        integer                 ninold
c                               ninold is the number of subintervals
c                               before the current remeshing.
c
        integer                 ninpre
c                               ninpre is the number of subintervals
c                               when ipar(icount) = 0 before remeshing.
c
        integer                 neqpre
c                               neqpre is the number of bspline
c                               coefficients when ipar(icount)=0 before
c                               remeshing when using radau_{kcol+1}.
c
        integer                 irold
c                               irold is the value of ipar(ixold) before
c                               remeshing.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
c
c-----------------------------------------------------------------------
c Direct pointers into the RPAR floating point work array:
        integer                 iiniss
        parameter              (iiniss =  2)
c                               rpar(iiniss) = the initial stepsize when
c                               mflag(4) = 1.
c
        integer                 ierrat
        parameter              (ierrat =  3)
c                               rpar(ierrat) = the value of the largest
c                               component of rpar(ipar(iercom)).
c
        integer                 it0
        parameter              (it0    =  4)
c                               rpar(it0)    = t0 at the last accepted
c                               time step.
c
        integer                 irpstr
        parameter              (irpstr = 11)
c                               rpar(1:irpstr-1) are reserved to store
c                               floating point scalar quantities.
c
c-----------------------------------------------------------------------
c Direct pointers into the IPAR integer work array:
        integer                 inpde
        parameter              (inpde  =  1)
c                               ipar(inpde) = npde
c
        integer                 ikcol
        parameter              (ikcol  =  2)
c                               ipar(ikcol) = kcol.
c
        integer                 inint
        parameter              (inint  =  3)
c                               ipar(inint) = nint.
c
        integer                 incpt
        parameter              (incpt =  4)
c                               ipar(incpt) = ncpts.
c
        integer                 ineq
        parameter              (ineq  =  5)
c                               ipar(ineq) = neq.
c
        integer                 iipstp
        parameter              (iipstp =  6)
c                               ipar(iipstp) = the minimum size of ipar.
c
        integer                 irpstp
        parameter              (irpstp =  7)
c                               ipar(irpstp) = the minimum size of rpar.
c
        integer                 iest
        parameter              (iest   =  8)
c                               ipar(iest) = mflag(5), for LOI or SCI.
c
        integer                 irshin
        parameter              (irshin =  9)
c                               ipar(irshin) is the number of remeshing
c                               times at the initial step.
c
        integer                 isteps
        parameter              (isteps = 10)
c                               ipar(isteps) is the number of time steps
c                               on the current problem.
c
        integer                 irmesh
        parameter              (irmesh = 11)
c                               ipar(irmesh) is the number of remeshing
c                               times after BACOLRI starts the initial
c                               step.
c
        integer                 istblc
        parameter              (istblc = 13)
c                               ipar(istblc) is the number of steps
c                               BACOLRI has taken before the latest cold
c                               start.
c
        integer                 irshfg
        parameter              (irshfg = 14)
c                               ipar(irshfg) is a flag for redefining
c                               all the pointers.
c                               ipar(irshfg) = 0, the initial step or
c                                                 any step not needing
c                                                 remesh;
c                                            = 1, a step needing remesh.
c
        integer                 icount
        parameter              (icount = 15)
c                               ipar(icount) is the number of remeshing
c                               times at the current step.
c
        integer                 istart
        parameter              (istart = 16)
c                               ipar(istart) is a flag to begin the
c                               code.
c                               ipar(istart) = 0, the initial step;
c                                            = 1, not the initial step.
c
        integer                 imflg6
        parameter              (imflg6 = 17)
c                               ipar(imflg6) stores the value of
c                               mflag(6) in order for it to be passed
c                               down to caljac.
c
        integer                 iradi
        parameter              (iradi  = 61)
c                               ipar(iradi) stores, before remeshing,
c                               the first 20 elements of the integer
c                               point work array in RADAU5.
c
        integer                 iiwork
        parameter              (iiwork = 81)
c                               ipar(iiwork) is the integer work array
c                               for RADAU5.
c
        integer                 ipivot
        parameter              (ipivot = 101)
c                               ipar(ipivot-1+i), i = 1, neq, contains
c                               the pivoting information from the
c                               factorization of the temporary matrix
c                               for RADAU5.
c
c-----------------------------------------------------------------------
c Indirect pointers into the RPAR floating point work array:
        integer                 ih
        parameter              (ih     = 21)
c                               rpar(ipar(ih)) stores the mesh step
c                               size sequence.
c
        integer                 ixcol
        parameter              (ixcol = 22)
c                               rpar(ipar(ixcol)) stores the
c                               collocation points when using
c                               radau_{kcol}.
c
        integer                 ixbs
        parameter              (ixbs  = 23)
c                               rpar(ipar(ixbs)) stores the breakpoint
c                               sequence when using radau_{kcol}.
c
        integer                 iy
        parameter              (iy    = 24)
c                               rpar(ipar(iy)) stores the vector of
c                               solution components to the DAE system
c                               when using radau_{kcol}.
c
        integer                 iabtp
        parameter              (iabtp = 26)
c                               rpar(ipar(iabtp)) stores the top block
c                               of the ABD collocation matrices when
c                               using radau_{kcol}.
c
        integer                 iabbk
        parameter              (iabbk = 27)
c                               rpar(ipar(iabbk)) stores the nint
c                               blocks in the middle of the ABD
c                               collocation matrices when using
c                               radau_{kcol}.
c
        integer                 iabbt
        parameter              (iabbt = 28)
c                               rpar(ipar(iabbt)) stores the bottom
c                               block of the ABD collocation matrices
c                               when using radau_{kcol}.
c
        integer                 irwork
        parameter              (irwork = 29)
c                               rpar(ipar(irwork)) stores the floating
c                               point work array for RADAU5. And it
c                               is also used to be a work storage for
c                               the subroutine INIY to get
c                               the initial condition.
c
        integer                 iwkrj
        parameter              (iwkrj  = 30)
c                               rpar(ipar(iwkrj)) stores an additional
c                               work array required by RES and JAC.
c
        integer                 ibasi
        parameter              (ibasi = 31)
c                               rpar(ipar(ibasi)) stores the basis
c                               function values at the collocation
c                               points when using radau_{kcol}.
c                               rpar(ipar(ibasi)) contains
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts). A(k,j,i) stores
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        integer                 iatol
        parameter              (iatol  = 32)
c                               rpar(ipar(iatol)) = atol.
c
        integer                 irtol
        parameter              (irtol  = 33)
c                               rpar(ipar(irtol)) = rtol.
c
        integer                 iexcol
        parameter              (iexcol = 34)
c                               rpar(ipar(iexcol)) stores the
c                               collocation points which are used for
c                               error estimate.
c
        integer                 iewts
        parameter              (iewts  = 35)
c                               rpar(ipar(iewts)) stores the gaussian
c                               weights which are used for error
c                               estimate.
c
        integer                 iebas
        parameter              (iebas = 36)
c                               rpar(ipar(iebas)) stores the values
c                               of the nonzero basis functions at
c                               rpar(ipar(iexcol)) when using
c                               radau_{kcol}.
c
        integer                 iecoef
        parameter              (iecoef = 37)
c           BACOLR --> BACOLRI: unsure if the documentation of this
c                               variable is correct
c                               iercoef holds the values of the nonzero
c                               B-spline basis functions and the
c                               Hermite-Birkhoff coefficients used in
c                               the selected interpolation subroutine.
c                               They are reusable provide the spatial
c                               mesh does not change.
c
        integer                 iercom
        parameter              (iercom = 38)
c                               rpar(ipar(iercom)) stores the error
c                               estimate for each component.
c
        integer                 ierint
        parameter              (ierint = 39)
c                               rpar(ipar(ierint)) stores the error
c                               estimate at each subinterval.
c
        integer                 iework
        parameter              (iework = 40)
c                               rpar(ipar(iework)) stores the floating
c                               point work array for errest.
c
        integer                 ix
        parameter              (ix     = 50)
c                               rpar(ipar(ix)) is the copy of x in rpar
c                               for use by errest.
c
        integer                 ixold
        parameter              (ixold  = 51)
c                               rpar(ipar(ixold)) stores the mesh point
c                               sequence when ipar(icount) = 0 before
c                               remeshing.
c
        integer                 iypre
        parameter              (iypre  = 53)
c                               rpar(ipar(iypre)) stores the values of
c                               rpar(ipar(iy)) at the previous step.
c                               It is required for a hot restart after
c                               remeshing.
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               colpnt
c                               RADAU5
c                               iniy
c                               meshsq
c                               reinit
c                               remesh
        external                radfcn
        external                radjac
        external                radmas
        external                solout
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c                               dcopy
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, Oct 18, 2006.
c
c-----------------------------------------------------------------------

c     Check validity of the used elements of the mflag vector.
      do 1 i = 1, 5
         if ((mflag(i) .lt. 0) .or. (mflag(i) .gt. 1)) goto 710
   1  continue
      if ((mflag(6) .lt. 0) .or. (mflag(6) .gt. 3)) goto 711

      ipar(irshfg) = 0

c     The SCI scheme estimates the error for the computed collocation
c     solution, while the LOI scheme estimates the error for a
c     collocationsolution of one order lower. (In LOI mode, the code
c     computes an estimate of the error that would be obtained were the
c     code to compute a collocation solution that was one order lower.
c     This is similar to local extrapolation for initial value
c     solutions.)
c     kerr keeps track of this for calls to meshsq and remesh.
      if (mflag(5) .eq. 0) then
c        LOI scheme
         kerr = kcol - 1
      elseif (mflag(5) .eq. 1) then
c        SCI scheme
         kerr = kcol
      endif

c     Check for continuation of a previous problem.
      if (mflag(1) .eq. 1) then
         ipar(istart) = 1

         neq   = ipar(ineq)
         lenpd = npde*npde*(nconti+nconti+kcol*(kcol+nconti)*nint)
         leniw = 20 + neq*3
         lenrw = 20 + 12*neq + 3*lenpd

         goto 200
      else

c        Check if the user specifies an initial stepsize
         if (mflag(4) .eq. 1) then
            if (((tout-t0)*rpar(iiniss)) .lt. zero) goto 720
            if (rpar(iiniss) .eq. zero) goto 725
         endif

         do 10 i = 1, iradi
            ipar(i) = 0
   10    continue
         do 20 i = irpstr, lrp
            rpar(i) = zero
   20    continue
         torign = t0
      endif

c-----------------------------------------------------------------------
c     On the initial call or after remeshing, check for valid input and
c     initialize the workspace.
c-----------------------------------------------------------------------

  100 continue

c     Check validity of npde, kcol, and nint.
      if (npde .le. 0) goto 730
      if ((kcol .le. 2) .or. (kcol .gt. mxkcol)) goto 740
      if ((nint .le. 0) .or. (nint .gt. nintmx)) goto 750

c     Check for a monotone mesh.
      do 110 i = 1, nint
         if (x(i) .ge. x(i+1)) goto 760
  110 continue

c-----------------------------------------------------------------------
c     Calculate the extra storage requirements of res and jac.
      lenrj = (6 + 5 * npde) * npde

c     Calculate the number of collocation points when using
c     radau_{kcol}.
      ncpts = nint * kcol + nconti

c     Calculate the number of DAEs when using radau_{kcol}.
      neq = npde * ncpts

c     Size of the ABD iteration matrix when using radau_{kcol}.
      lenpd = npde*npde*(nconti+nconti+kcol*(kcol+nconti)*nint)

c     Calculate the extra storage requirements of iniy when
c     using radau_{kcol}.
C     BACOLRI --> BACOLRIVEC
C       lenin = lenpd + 2 * neq + 2 * npde * (2 + npde)
      lenin = lenpd + 2 * neq + 2 * npde * (2 + npde)
     &              + 3 * npde * nint * kcol

c-----------------------------------------------------------------------
c     Total size of the RADAU5 floating point work array.
      lenrw = 20 + 12 * neq + 3 * lenpd

c     Total size of the RADAU5 integer work array.
      leniw = 20 + neq * 3

c-----------------------------------------------------------------------
c     Calculate the number of quadrature points used for error
c     estimate and the extra storage requirements of errest.
c     The SCI scheme needs slightly more storage here.
      if (mflag(5) .eq. 0) then
         necpts = (kcol + 2) * nint
         lenerr = (2 * necpts + nint) * npde
     &          + ((kcol - 3) * nint + (nint + 1) * 2) * npde
         lencof = (kcol + 1) * (kcol + 2)
     &          + (kcol + nconti) * (kcol - 3) * nint
     &          + (kcol + nconti) * 2 * (nint + 1)
      elseif (mflag(5) .eq. 1) then
         necpts = (kcol + 3) * nint
         lenerr = (2 * necpts + nint) * npde
     &          + ((kcol - 2) * nint + (nint + 1) * 2) * npde
         lencof = (kcol + 4) * (kcol + 3) * nint
     &          + (kcol + nconti) * (kcol - 2) * nint
     &          + (kcol + nconti) * 2 * (nint + 1)
      endif

c-----------------------------------------------------------------------
c     Save the input parameters in ipar, the integer communication
c     storage array.
      ipar(inpde)  = npde
      ipar(ikcol)  = kcol
      ipar(inint)  = nint
      ipar(incpt) = ncpts
      ipar(ineq)  = neq
      ipar(iest)   = mflag(5)

c     mflag(6) is loaded into ipar so that it may be communicated to
c     caljac. This value can in theory change between calls to BACOLRI.
      ipar(imflg6) = mflag(6)

c-----------------------------------------------------------------------
c     Calculate the offsets into rpar, the floating point storage array.
c-----------------------------------------------------------------------
      ipar(iatol)  = irpstr
      ipar(irtol)  = ipar(iatol)  + npde

      ipar(ih)     = ipar(irtol)  + npde

      ipar(iy)    = ipar(ih)     + nint

      ipar(ixcol) = ipar(iy)    + neq
      ipar(ixbs)  = ipar(ixcol) + ncpts
      ipar(iabtp) = ipar(ixbs)  + ncpts + kcol + nconti
      ipar(iabbk) = ipar(iabtp) + npde * npde * nconti
      ipar(iabbt) = ipar(iabbk) + npde * npde * nint * kcol
     &                              * (kcol + nconti)
      ipar(ibasi) = ipar(iabbt) + npde * npde * nconti

      ipar(irwork) = ipar(ibasi) + (kcol + nconti) * 3 * ncpts
C     BACOLRI --> BACOLRIVEC
C       ipar(iwkrj)  = ipar(irwork) + lenrw
      ipar(iwkrj)  = ipar(irwork) + lenrw + 3*npde*(kcol*nint+1)

C     BACOLRI --> BACOLRIVEC
C       ipar(iexcol) = ipar(iwkrj)  + lenrj
      ipar(iexcol) = ipar(iwkrj)  + lenrj + 3*npde*(kcol*nint+1)
      ipar(iewts)  = ipar(iexcol) + necpts
      ipar(ierint) = ipar(iewts)  + necpts
      ipar(iercom) = ipar(ierint) + nint
      ipar(iebas) = ipar(iercom) + npde
      ipar(iecoef) = ipar(iebas) + (kcol + nconti) * necpts
      ipar(iework) = ipar(iecoef) + lencof

      ipar(ix)     = ipar(iework) + lenerr
      ipar(ixold)  = ipar(ix)     + nintmx + 1

      ipar(iypre)  = ipar(ixold)  + nintmx + 1

c     The offset is different between the initial call and remeshing.
      if ((ipar(irshfg) .ne. 0) .and. (ipar(istart) .eq. 1)) then
         ipar(irpstp) = ipar(iypre) + neqpre - 1
      else
         ipar(irpstp) = ipar(iypre) + neq - 1
      endif

c     Check for a sufficiently large rpar floating point work array.
      if (lrp .lt. ipar(irpstp)) goto 770

c     Calculate the offsets into the integer storage array.
      ipar(iipstp) = ipivot + 3 * neq - 1

c     Check for a sufficiently large ipar integer work array.
      if (lip .lt. ipar(iipstp)) goto 780

c     Check whether it is initial call or for remeshing.
      if ((ipar(irshfg) .ne. 0) .and. (ipar(istart) .ne. 0)) goto 300

c-----------------------------------------------------------------------
c     Save the atol and rtol in rpar, the real communication
c     storage array.
      if (mflag(2) .eq. 0) then
         rpar(ipar(iatol)) = atol(1)
         rpar(ipar(irtol)) = rtol(1)
      else
         do 120 i = 1, npde
            rpar(ipar(iatol)-1+i) = atol(i)
            rpar(ipar(irtol)-1+i) = rtol(i)
  120    continue
      endif

c-----------------------------------------------------------------------
c     Perform initializations for using radau_{kcol}.
c-----------------------------------------------------------------------
c     Check whether it is initial call or for remeshing.
      if (ipar(icount) .eq. 0) then

c        Set the initial stepsize if applicable.
         if (mflag(4) .eq. 0) then
            if (mflag(2) .eq. 0) then
               rpar(iiniss) = max(atol(1),rtol(1))
            else
               rpar(iiniss) = zero
               do 130 i = 1, npde
                  rpar(iiniss) = max(rpar(iiniss), atol(i))
  130          continue
               do 140 i = 1, npde
                  rpar(iiniss) = max(rpar(iiniss), rtol(i))
  140          continue
            endif
            if (((tout-t0)*rpar(iiniss)) .lt. zero)
     &         rpar(iiniss) = - rpar(iiniss)
         endif
      else
         ipar(irshin) = ipar(irshin) + 1
      endif
C      write(*,*) 'calling meshq'
      call meshsqri(kerr, nint, x, rpar(ipar(irwork)), rpar(ipar(ih)),
     &            rpar(ipar(iexcol)), rpar(ipar(iewts)))
c      write(*,*) 'calling colpnt'
      call colpntri(kcol, nint, ncpts, x, rpar(ipar(ih)),
     &            rpar(ipar(irwork)), rpar(ipar(ixcol)),
     &            rpar(ipar(ixbs)))

      icflag = 0

C     BACOLRI --> BACOLRIVEC
C       call iniyri(t0, npde, kcol, nint, neq, ncpts, mflag(3),
C      &          rpar(ipar(ixcol)), rpar(ipar(ixbs)),
C      &          rpar(ipar(iabbk)), rpar(ipar(ibasi)), rpar(ipar(iy)),
C      &          ipar(ipivot), rpar(ipar(irwork)), lenin, icflag,
C      &          mflag(6), f, bndxa, difbxa, bndxb, difbxb, uinit)
      call iniyri(t0, npde, kcol, nint, neq, ncpts, mflag(3),
     &          rpar(ipar(ixcol)), rpar(ipar(ixbs)),
     &          rpar(ipar(iabbk)), rpar(ipar(ibasi)), rpar(ipar(iy)),
     &          ipar(ipivot), rpar(ipar(irwork)), lenin, icflag,
     &          mflag(6), bndxa, difbxa, bndxb, difbxb, uinit,
     &          uinitvec)

      if (icflag .ne. 0) then
         idid = -66
         goto 600
      endif

      ipar(irshfg) = 0


c     Copy rpar(ipar(iy)) to rpar(ipar(iypre)).
c      write(*,*) 'calling dcopy'
      call dcopyri(neq, rpar(ipar(iy)), 1, rpar(ipar(iypre)), 1)

      do 150 i = 1, 20
         ipar(iiwork-1+i) = 0
  150 continue

      do 160 i = 1, 20
         rpar(ipar(irwork)-1+i) = zero
  160 continue

      goto 400

c-----------------------------------------------------------------------
c     This is not the first call for the problem, and integration is to
c     continue.
c-----------------------------------------------------------------------
  200 continue

c     Examine idid to determine if RADAU5 can be called again.
      if (idid .ne. 1) goto 790

      goto 400

c-----------------------------------------------------------------------
c     Initialization after remeshing.
c-----------------------------------------------------------------------
  300 continue

      ipar(irmesh) = ipar(irmesh) + 1

      do 310 i = 1, 20
         ipar(iiwork-1+i) = 0
  310 continue
      do 320 i = 1, 20
         rpar(ipar(irwork)-1+i) = zero
  320 continue

      rpar(ipar(irwork)-1+3) = point1

c     Move the last two sections of rpar (indexed by ipar(ixold)
c     and ipar(iypre)) to their new positions after a remeshing
c     in such a way as to avoid overwriting what you just copied.
c     The new and old ranges may overlap.
      if (nint .lt. ninold) then
c        rpar is shrinking; copy from the front
c         write(*,*) 'calling dcopy'
         call dcopyri(nintmx+1+neqpre, rpar(irold),
     &              1, rpar(ipar(ixold)), 1)
      elseif (nint .gt. ninold) then
c        rpar is growing; copy in reverse
c         write(*,*) 'calling dcopy'
         call dcopyri(nintmx+1+neqpre, rpar(irold),
     &              -1, rpar(ipar(ixold)), -1)
      endif

      lenri  = lenpd + kcol + nconti + kcol * (ninpre + 1)
     &         + 2 * nconti

c      write(*,*) 'calling meshq'
      call meshsqri(kerr, nint, x, rpar(ipar(irwork)+20),
     &            rpar(ipar(ih)),
     &            rpar(ipar(iexcol)), rpar(ipar(iewts)))

c      write(*,*) 'calling reinit'
      call reinitri(npde, kcol, kcol, nint, ninpre, ncpts, neq,
     &            neqpre, x, rpar(ipar(ixold)),
     &            rpar(ipar(iypre)), rpar(ipar(irwork)+20), lenri,
     &            ipar(ipivot), rpar(ipar(ih)), rpar(ipar(ixbs)),
     &            rpar(ipar(ixcol)), rpar(ipar(ibasi)),
     &            rpar(ipar(iy)), rpar(ipar(iabbk)), icflag)

      if (icflag .ne. 0) then
         idid = -66
         goto 600
      endif

c-----------------------------------------------------------------------
c     Time integration loop for RADAU5.
c-----------------------------------------------------------------------

  400 continue

c     Update the copy of x in rpar for use by errest.
      do 401 i = 1, nint + 1
         rpar(ipar(ix)-1+i) = x(i)
  401 continue

c      write(*,*) 'calling radau5'
      call radau5(neq, radfcn, t0, rpar(ipar(iy)), tout, rpar(iiniss),
     &            rtol, atol, mflag(2), radjac, radmas, solout,
     &            rpar(ipar(irwork)), lenrw, ipar(iiwork), leniw, rpar,
     &            ipar, cpar, idid, f, fvec, derivf, bndxa, difbxa,
     &            bndxb,
     &            difbxb, uinit)
c      write(*,*) 'leaving radau5'

c-----------------------------------------------------------------------
c     Check for a successful time step and decide whether to continue
c     integration or to perform a remeshing.
c-----------------------------------------------------------------------

      if (idid .le. 0) goto 600

      if (idid .eq. 2) then

c        The current step is rejected. A greater number or resmeshings
c        is to be allowed at the first step, when the user-provided
c        mesh may be entirely arbitrary.
         if (ipar(istart) .eq. 1 .and. ipar(icount) .eq. maxrsh .or.
     &   ipar(istart) .eq. 0 .and. ipar(icount) .eq. mxrhin) goto 610

c        For the first remeshing at the current step, save nintpre and
c        neqpre at the last successful step.
         if (ipar(icount) .eq. 0) then
            ninpre = nint
            neqpre = neq
         endif

         do 410 i = 14, 20
            ipar(iradi-1+i) = ipar(iiwork-1+i)
  410    continue
         ninold = nint
         irold = ipar(ixold)

c        Update xold.
         if (ipar(icount) .eq. 0) then
            do 420 i = 1, ninpre + 1
               rpar(ipar(ixold)-1+i) = x(i)
  420       continue
         endif
c         write(*,*) 'calling remesh'
         call remeshri(ipar(istart), ipar(icount), nintmx,
     &               ninold, rpar(ierrat), rpar(ipar(ierint)),
     &               ipar(irshfg), nint, kerr, x, rpar(ipar(iework)))

         if (ipar(istart) .eq. 1) then

c           This is not the initial step.
            t0 = rpar(it0)

            ipar(istblc) = ipar(istblc) + ipar(iradi-1+17) - 1

         else

c           This is the initial step.
            t0 = torign

         endif

         goto 100

      else

c        The current step is the last step.
         goto 500

      endif

c-----------------------------------------------------------------------
c     Successful return section.
c-----------------------------------------------------------------------
  500 continue

c     Retrieve the value of mflag(1).
      mflag(1) = 1

c     Retrieve the output vector y from the rpar communication array.
      do 510 i = 1, neq
         y(i) = rpar(ipar(iy)-1+i)
  510 continue

c     Retrieve information on the time stepping from the ipar array.
      ipar(isteps) = ipar(istblc) + ipar(iiwork-1+17)

      return

c-----------------------------------------------------------------------
c     Unsuccessful return section.
c-----------------------------------------------------------------------
  600 continue
      write(6,9999) 'ERROR: BACOLRI runtime error in time stepping.'
      write(6,9999) '       An error code and message should have'
      write(6,9999) '       been issued by RADAU5.'
      return
  610 continue
      if (ipar(istart) .eq. 1) then
         write(6,9998) 'ERROR: BACOLI has remeshed ', maxrsh
     &                 , ' times at', ' t0 =', rpar(it0)
      else
         write(6,9998) 'ERROR: BACOLI has remeshed ', maxrsh
     &                 , ' times at', ' t0 =', torign
      endif
      idid = -41
      return

c-----------------------------------------------------------------------
c     The following section is the return point for invalid input.
c-----------------------------------------------------------------------

  710 continue
      write(6,9999) 'ERROR: BACOLRI input violation.'
      write(6,9999) 'Require:  0 <= mflag(i) <= 1, i = 1, 2, ..., 6.'
      idid = -51
      return
  711 continue
      write(6,9999) 'ERROR: BACOLRI input violation.'
      write(6,9999) ' Require:  0 <= mflag(6) <= 3'
      idid = -51
      return
  720 continue
      write(6,9999) 'ERROR: BACOLRI input violation.'
      write(6,9999) 'Require:  if mflag(4) = 1, tout must be in front'
      write(6,9999) 'of t0.'
      idid = -53
      return
  725 continue
      write(6,9999) 'ERROR: BACOLRI input violation.'
      write(6,9999) 'Require:  if mflag(4) = 1, rpar(2) must be the'
      write(6,9999) 'initial stepsize, thus nonzero.'
      idid = -54
      return
  730 continue
      write(6,9999) 'ERROR: BACOLRI input violation.'
      write(6,9999) 'Require: npde > 0.'
      idid = -55
      return
  740 continue
      write(6,9999) 'ERROR: BACOLRI input violation.'
      write(6,9999) 'Require: 2 < kcol <=', mxkcol, '.'
      idid = -56
      return
  750 continue
      if (istart .eq. 1) then
          write(6,9998) 'ERROR: nint >', nintmx, ' at'
     &                  , ' t0 =' ,rpar(it0)
      else
          write(6,9998) 'ERROR: nint >', nintmx, ' at'
     &                  , ' t0 =', torign
      end if
      idid = -57
      return
  760 continue
      write(6,9999) 'ERROR: BACOLRI input violation.'
      write(6,9999) 'Require: x(1) < x(2) < ... < x(nint+1).'
      idid = -58
      return
  770 continue
      write(6,9999) 'ERROR: BACOLRI input violation.'
      write(6,9999) 'Require: lrp >= ', ipar(irpstp), '.'
      idid = -59
      return
  780 continue
      write(6,9999) 'ERROR: BACOLRI input violation.'
      write(6,9999) 'Require: lip >= ', ipar(iipstp), '.'
      idid = -60
      return
  790 continue
      write(6,9999) 'ERROR: BACOLRI input violation.'
      write(6,9999) 'IDID .ne. 1, on a continuation call of BACOLRI'
      write(6,9999) 'If IDID > 1, set idid = 1 and tout (t0 < tout)'
      write(6,9999) 'If IDID < -1, the code cannot be continued due to'
      write(6,9999) '              a previous error.'
      idid = -61
      return

c-----------------------------------------------------------------------
 9998 format(a,i4,a,a,e12.5)
 9999 format(a,i8,a,i4,a,i4,a,i4,a,i4)
c-----------------------------------------------------------------------
      end
      subroutine calfcn(npde, kcol, nint, ncpts, neq, xcol, fbasis,
     &                      t, y, work, fr, f, fvec, bndxa, bndxb)

c-----------------------------------------------------------------------
c Purpose:
c       This subroutine is called by radfcn. It provides a lower-level
c       interface to calculate the right side of the DAEs.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, May 21, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        integer                 neq
c                               neq=npde*ncpts is the number of bsplines
c                               coefficients (or DAEs).

        external                f
        external                fvec
        external                bndxa
        external                bndxb

c
        double precision        xcol(ncpts)
c                               xcol stores the collocation
c                               points when using kcol collocation
c                               points at each subinterval.
c
        double precision        fbasis((kcol+nconti)*3*ncpts)
c                               fbasis stores the basis function values
c                               at the collocation points. It acts like
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts). A(k,j,i) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        double precision        t
c                               T is the current time.
c
        double precision        y(neq)
c                               y is the vector of bspline
c                               coefficients at the current time.
c
c       Work storage:
c     BACOLRI --> BACOLRIVEC
C         double precision        work(4*npde+2*npde*npde)
        double precision        work(3*npde*nint*kcol)
c                               work is a floating point work array
c                               of size 4*npde+2*npde*npde.
c
c       Output:
        double precision        fr(neq)
c                               fr is the vector at the right side of
c                               the DAEs.
c
c-----------------------------------------------------------------------
c       Loop indices:
        integer                 i
        integer                 j
c       BACOLRI --> BACOLRIVEC
        integer                 k
c
c       Indices:
        integer                 ii
        integer                 jj
        integer                 mm
        integer                 kk
c
c-----------------------------------------------------------------------
c       Pointers into the floating point work array work:
        integer                 iu
c                               work(iu) stores the approximation to
c                               u(t,x).
c
        integer                 iux
c                               work(iux) stores the approximation to
c                               the first spatial derivative of u(t,x).
c
        integer                 iuxx
c                               work(iuxx) stores the approximation to
c                               the second spatial derivative of u(t,x).
c
c       BACOLRI --> BACOLRIVEC
c       ----------------------
        integer                 swapiu
        integer                 swapiux
        integer                 swapiuxx
c                               swap memory for vectorized calls.
        integer                 vnpts
c                               number of collocation points which the
c                               vectorized f user routine must evaluate
        integer                 voffset
c                               Offset into work arrays for vectorized
c                               BACOLRI
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bndxa
c                               bndxb
c                               f
c                               eval
c-----------------------------------------------------------------------

c     BACOLRI --> BACOLRIVEC
      vnpts = nint * kcol

c     Set pointers into the temporary floating point work array.
      iu     = 1
      iux    = iu     + npde
      iuxx   = iux    + npde
c     BACOLRI --> BACOLRIVEC
      swapiu   = iuxx    + npde
      swapiux  = swapiu  + npde*vnpts
      swapiuxx = swapiux + npde*vnpts

c-----------------------------------------------------------------------
c     Loop over the nint blocks of collocation equations.

      do 20 i = 1, nint

c        ii is the value of ileft for the current collocation point.
         ii = kcol + nconti + (i - 1) * kcol

         do 10 j = 1, kcol

c           jj is the pointer of collocation point.
            jj = (i - 1) * kcol + j + 1

c           mm is the pointer of fr.
            mm = (jj - 1) * npde + 1

c           BACOLRI --> BACOLRIVEC
c           set offsets into vectorized work arrays
            voffset = j - 1 + kcol*(i-1)

c           kk is the pointer of the basis function values at
c           the current collocation point.
            kk =(jj-1)*(kcol+nconti)*3+1

c           Generate the approximate solution and its spatial
c           derivatives at the current collocation point.
            call evalri(npde,kcol,ii,jj,ncpts,work(iu),work(iux),
     &                work(iuxx),fbasis(kk),y)

c           Evaluate the function defining the PDE at the current
c           collocation point, storing the result in fr.
c           BACOLRI --> BACOLRIVEC
C             call f(t, xcol(jj), work(iu), work(iux),
C      &              work(iuxx), fr(mm), npde)
            do 181 k = 1, npde
                work(swapiu+(k-1)*vnpts+voffset) = work(iu+(k-1))
                work(swapiux+(k-1)*vnpts+voffset) = work(iux+(k-1))
                work(swapiuxx+(k-1)*vnpts+voffset) 
     &                  = work(iuxx+(k-1))
  181       continue

   10    continue
   20 continue

c     BACOLRI --> BACOLRIVEC
      call fvec(t, xcol(2), work(swapiu), work(swapiux),
     &       work(swapiuxx), fr(npde+1), vnpts, npde)

c     Copy all of the yprime values into the swap memory (swapiu
c     arbitrarily chosen).
      call dcopy(nint*kcol*npde, fr(npde+1), 1, work(swapiu), 1)

c     Put the points all in back in their place.
      do 201 i = 1, vnpts
          do 211 j = 1, npde
              fr(npde+npde*(i-1)+j)=work(swapiu+(j-1)*vnpts+i-1)
  211     continue
  201 continue
c-----------------------------------------------------------------------
c     Calculate (fr(i), i=1, npde), which depend on the left
c     boundary point.
      call evalri(npde, kcol, kcol+2, 1, ncpts, work(iu), work(iux),
     &          work(iuxx), fbasis(1), y)
      call bndxa(t, work(iu), work(iux), fr(1), npde)

c-----------------------------------------------------------------------
c     Calculate (fr(i), i=neq-npde+1, neq), which depend on the right
c     boundary point.
      call evalri(npde, kcol, ncpts, ncpts, ncpts, work(iu), work(iux),
     &          work(iuxx), fbasis(1+(ncpts-1)*(kcol+nconti)*3), y)
      call bndxb(t, work(iu), work(iux), fr(neq-npde+1), npde)

      return
      end
      subroutine caljacri(npde, kcol, nint, ncpts, neq, xcol, fbasis,
     &                  abdtop, abdbot, t, y, work, dfdy, ifgfdj,
     &                  f, derivf, bndxa, difbxa, bndxb, difbxb)

c-----------------------------------------------------------------------
c Purpose:
c       This subroutine is called by radjac. It provides a lower-level
c       interface to generate the Jacobian matrix at the right side
c       of the DAEs.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, May 22, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcolis the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        integer                 neq
c                               neq=npde*ncpts is the number of bspline
c                               coefficients (or DAEs).
        integer                 ifgfdj
c                               Are finite difference approximations
c                               being used to approximate Jacobian 
c                               matrices? Set to 0 to approximate 
c                               derivf, difbxa and difbxb.
c                               Set to 1 to approxiate only difbxa and 
c                               difbxb.
c                               Set to 2 to approximate only derivf.
c                               Set to 3 if all of derivf, difbxa and
c                               difbxb are provided.

        external                f
        external                derivf
        external                bndxa
        external                difbxa
        external                bndxb
        external                difbxb

c
        double precision        xcol(ncpts)
c                               xcol stores the collocation
c                               points when using kcol collocation
c                               points at each subinterval.
c
        double precision        abdtop(npde*npde*nconti)
c                               abdtop stores the top block of the ABD
c                               matrices.
c
        double precision        abdbot(npde*npde*nconti)
c                               abdbot stores the bottom block of the
c                               ABD matrices.
c
        double precision        fbasis((kcol+nconti)*3*ncpts)
c                               fbasis stores the basis function values
c                               at the collocation points. It acts like
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts). A(k,j,i) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        double precision        t
c                               t is the current time.
c
        double precision        y(neq)
c                               y is the vector of bspline
c                               coefficients at the current time.
c
c       Work storage:
        double precision        work(2*npde+5*npde*npde)
c                               work is a floating point work array
c                               of size 4*npde+7*npde*npde.
c
c       Output:
        double precision        dfdy(npde*npde*(2*nconti
     *                               +nint*kcol*(kcol+nconti)))
c                               dfdy is the ABD Jacobian matrix at the
c                               right side of the DAEs.
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 ifytop
c                               ifytop is the pointer into dfdy where
c                               the top block of the ABD Jacobian is
c                               stored.
c
        integer                 ifyblk
c                               ifyblk is the pointer into dfdy where
c                               the nint blocks in the middle of the ABD
c                               Jacobian are stored.
c
        integer                 ifybot
c                               ifybot is the pointer into dfdy where
c                               the bottom block of the ABD Jacobian is
c                               stored.
c
        integer                 nsiztb
c                               nsiztb is the size of the top block
c                               as same as the bottom block of the ABD
c                               Jacobian.
c
        integer                 nsizbk
c                               nsizbk is the size of a subblock in
c                               the middle of ABD Jacobian.
c
c-----------------------------------------------------------------------
c       Loop indices:
        integer                 i
        integer                 j
        integer                 k
        integer                 m
        integer                 n
c
        integer                 ii
        integer                 ij
        integer                 jj
        integer                 kk
        integer                 nn
        integer                 mm
        integer                 jk
        integer                 jk2
        integer                 jk3
        integer                 mn
        integer                 mn2
        integer                 mn3
c
c-----------------------------------------------------------------------
c       Pointers into the floating point work array work:
        integer                 iu
c                               work(iu) stores the approximation to
c                               u(t,x).
c
        integer                 iux
c                               work(iux) stores the approximation to
c                               the first spatial derivative of u(t,x).
c
        integer                 iuxx
c                               work(iuxx) stores the approximation to
c                               the second spatial derivative of u(t,x).
c
        integer                 idfdu
c                               work(idfdu) stores the Jacobian of f
c                               with respect to u.
c
        integer                 idfdux
c                               work(idfdux) stores the Jacobian of f
c                               with respect to u_x.
c
        integer                 idfuxx
c                               work(idfuxx) stores the Jacobian of f
c                               with respect to u_xx.
c
        integer                 idbdu
c                               work(idbdu-1+i), i=1, npde*npde,
c                               contains dbdu(npde,npde). That is,
c                               dbdu(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the unknown function u.
c
        integer                 idbdux
c                               work(idbdux-1+i), i=1, npde*npde,
c                               contains dbdux(npde,npde), That is,
c                               dbdux(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the spatial derivative of the
c                               unknown function u.
c
        integer                 idbdt
c                               work(idbdt-1+i), i=1, npde, contains
c                               the partial derivative of the i-th
c                               component of the vector b with respect
c                               to time t.
c
        integer                 ifdwrk
c                               work(ifdwrk-1+i), i=1, 2*npde, is used
c                               as work storage fdderivf and fdbndx.
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               derivf
c                               fdderivf
c                               difbxa
c                               difbxb
c                               fdbndx
c                               eval
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c                               dcopy
c
c-----------------------------------------------------------------------
 
c     Set pointers into the temporary floating point work array.
      iu     = 1
      iux    = iu     + npde
      iuxx   = iux    + npde
      idfdu  = iuxx   + npde
      idfdux = idfdu  + npde * npde
      idfuxx = idfdux + npde * npde
      idbdu  = idfuxx + npde * npde
      idbdux = idbdu  + npde * npde
      idbdt  = idbdux + npde * npde
      ifdwrk = idbdt  + npde

c     Set the indices into dfdy which define the ABD Jacobian.
      ifytop = 1
      ifyblk = ifytop + nconti * npde * npde
      ifybot = ifyblk + nint * npde * npde * kcol * (kcol + nconti)

c-----------------------------------------------------------------------
c     Calculate the size of top (or bottom) block and the size of a
c     subblock in the middle.
      nsiztb = npde * npde * nconti
      nsizbk = npde * npde * kcol * (kcol + nconti)

c     Initialize fytop, fyblk and fybot to zero.
      do 10 i = 1, nsiztb
         dfdy(ifytop-1+i) = zero
         dfdy(ifybot-1+i) = zero
   10 continue
      do 20 i = 1, nint * nsizbk
         dfdy(ifyblk-1+i) = zero
   20 continue

c-----------------------------------------------------------------------
c     Loop over the nint blocks of collocation equations and
c     caluculate the portion of dfdy which depends on them.

      do 70 i = 1, nint

c        ii+1 is the pointer to the first element at the i-th subblock
c        of the jacobian matrix, i = 1, nint.
         ii = ifyblk - 1 + (i - 1) * nsizbk

c        ij is the value of ileft for the current collocation point.
         ij = kcol + nconti + (i - 1) * kcol

         do 60 j = 1, kcol

c           jj+1 is the pointer to the first element corresponding to
c           the j-th collocation point in the i-th interval.
            jj = ii + (j - 1) * npde

c           mm is the index of the current collocation point.
            mm = (i - 1) * kcol + j + 1

c           Generate the approximate solution and its spatial
c           derivatives at the current collocation point.
            call evalri(npde,kcol,ij,mm,ncpts,work(iu),work(iux),
     &                work(iuxx),fbasis(1+(mm-1)*(kcol+nconti)*3),y)

c           Generate dfdu, dfdux, and dfdux at the current
c           collocation point (the j-th point of the i-th
c           subinterval).
            if (ifgfdj .eq. 0 .or. ifgfdj .eq. 2) then
                call fdderivf(t, xcol(1+(i-1)*kcol+j), work(iu),
     &                        work(iux), work(iuxx), work(idfdu),
     &                        work(idfdux), work(idfuxx), npde,
     &                        f, work(ifdwrk))
            else
                call derivf(t, xcol(1+(i-1)*kcol+j), work(iu),
     &                      work(iux), work(iuxx), work(idfdu),
     &                      work(idfdux), work(idfuxx), npde)
            end if

            do 50 k = 1, kcol + nconti

c              kk+1 is the pointer to the first element of a npde by
c              npde submatrix, which is corresponding to the j-th
c              collocation point in the i-th interval, and the k-th
c              nonzero basis function.
               kk = jj + (k-1) * npde * npde * kcol

c              jk is the pointer to the k-th nonzero function at the
c              mm-th collocation point in the basis function,
c              fbasis(1).
               jk = (mm - 1) * (kcol + nconti) * 3 + k

c              jk2 is the pointer to the first derivative for the
c              above basis function.
               jk2 = jk + kcol + nconti

c              jk3 is the pointer to the second derivative for the
c              above basis function.
               jk3 = jk2 + kcol + nconti

               do 40 m = 1, npde
                  do 30 n = 1, npde

c                    nn is the pointer to the (n, m) element of the
c                    npde by npde submatrix.
                     nn = kk + (m-1)*npde*kcol + n

c                    mn is the pointer to the (n, m) element of dfdu.
                     mn = idfdu - 1 + (m - 1) * npde + n

c                    mn2 is the pointer to the (n, m) element of dfdux.
                     mn2 = mn + npde * npde

c                    mn3 is the pointer to the (n, m) element of dfduxx.
                     mn3 = mn2 + npde * npde

c                    now set up the value in pd at the place nn.
                     dfdy(nn) = work(mn) * fbasis(jk)
     &                          + work(mn2) * fbasis(jk2)
     &                          + work(mn3) * fbasis(jk3)

   30             continue
   40          continue
   50       continue
   60    continue
   70 continue

c-----------------------------------------------------------------------
c     Update the values at the left boundary.
      call evalri(npde, kcol, kcol+2, 1, ncpts, work(iu), work(iux),
     &          work(iuxx), fbasis(1), y)
      if (ifgfdj .lt. 2) then
          call fdbndxri(t, work(iu), work(iux), work(idbdu),
     &                work(idbdux), work(idbdt), npde,
     &                bndxa, work(ifdwrk))
      else
          call difbxa(t, work(iu), work(iux), work(idbdu),
     &                work(idbdux), work(idbdt), npde)
      end if

c     Update the top block of the collocation matrix dG/dY'.
      do 90 j = 1, npde
         do 80 i = 1, npde
            ii = (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            abdtop(jj) =
     &            fbasis(2+kcol+nconti) * work(idbdux-1+mm)
            abdtop(ii) =
     &            work(idbdu-1+mm) - abdtop(jj)
   80    continue
   90 continue

c-----------------------------------------------------------------------
c     Update the values at the right boundary.
      call evalri(npde, kcol, ncpts, ncpts, ncpts, work(iu), work(iux),
     &          work(iuxx), fbasis(1+(ncpts-1)*(kcol+nconti)*3), y)
      if (ifgfdj .lt. 2) then
          call fdbndxri(t, work(iu), work(iux), work(idbdu),
     &                work(idbdux), work(idbdt), npde,
     &                bndxb, work(ifdwrk))
      else
          call difbxb(t, work(iu), work(iux), work(idbdu),
     &                work(idbdux), work(idbdt), npde)
      end if

c     Update the bottom block of the collocation matrix.
      do 110 j = 1, npde
         do 100 i = 1, npde
            ii = (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            abdbot(ii) =
     &            fbasis(1+kcol+kcol+nconti+(ncpts-1)*(kcol+nconti)*3)
     &            * work(idbdux-1+mm)
            abdbot(jj) =
     &            work(idbdu-1+mm) - abdbot(ii)
  100    continue
  110 continue

c-----------------------------------------------------------------------
c     Copy abdtop and abdbot to the corresponding parts of dfdy.

      call dcopyri(nsiztb, abdtop, 1, dfdy(ifytop), 1)
      call dcopyri(nsiztb, abdbot, 1, dfdy(ifybot), 1)

c-----------------------------------------------------------------------
      return
      end
      subroutine calmas(npde, kcol, nint, abdblk, am)
c-----------------------------------------------------------------------
c Purpose:
c       This subroutine is called by radmas. It provides a lower-level
c       interface to generate the mass-matrix.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, May 20, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcolis the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        double precision        abdblk(npde*npde*nint*kcol
     *                                 *(kcol+nconti))
c                               abdblk stores the nint
c                               blocks in the middle of the ABD
c                               collocation matrices when using kcol
c                               collocation points at each subinterval.
c
c       Output:
        double precision        am(npde*npde*(2*nconti
     *                             +nint*kcol*(kcol+nconti)))
c                               am is the ABD mass-matrix.
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 iamtop
c                               iamtop is the pointer into am where the
c                               top block of the ABD matrix is stored.
c
        integer                 iamblk
c                               iamblk is the pointer into am where the
c                               nint blocks in the middle of the ABD
c                               matrix are stored.
c
        integer                 iambot
c                               iambot is the pointer into am where the
c                               bottom block of the ABD matrix is
c                               stored.
c
        integer                 nels
c                               nels is the size of a subblock in
c                               the middle of ABD Jacobian.
c
c-----------------------------------------------------------------------
c       Loop indices:
        integer                 i
c
c-----------------------------------------------------------------------
      nels = npde*npde*kcol*(kcol+nconti)

c     Set the pointers into the mass-matrix.
      iamtop = 1
      iamblk = iamtop + npde*npde*nconti
      iambot = iamblk + nint*nels

c     Initialize the top and bottom of the mass-matrix to zero.
      do 10 i = 1, npde*npde*nconti
         am(i) = zero
   10 continue
      do 20 i = iambot, iambot+npde*npde*nconti-1
         am(i) = zero
   20 continue

c-----------------------------------------------------------------------
c     The nint blocks at the middle of the mass-matrix are set up.
      do 30 i = 1, nint*nels
         am(i+npde*npde*nconti) = abdblk(i)
   30 continue

      return
      end
      subroutine colpntri(kcol, nint, ncpts, x, h, work, xcol, xbs)

c-----------------------------------------------------------------------
c Purpose:
c       This routine generates the piecewise polynomial space breakpoint
c       sequence, and calculates the collocation point sequence.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, April 3, 2001.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        integer                 mxkcol
        parameter              (mxkcol = 10)
c                               mxkcol is the maximum number of
c                               collocation points per subinterval.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        double precision        x(nint+1)
c                               x is the spatial mesh which divides the
c                               interval [x_a, x_b] as: x_a = x(1) <
c                               x(2) < x(3) < ... < x(nint+1) = x_b.
c
        double precision        h(nint)
c                               h is the mesh step size sequence.
c
c       Work Storage:
        double precision        work(kcol*kcol)
c                               work is a floating point work storage
c                               array of size lw.
c
c       Output:
        double precision        xcol(ncpts)
c                               The sequence of collocation points on
c                               the interval [a,b].
c
        double precision        xbs(ncpts+kcol+nconti)
c                               The breakpoint sequence.
c                               xbs(i)=x(1), i=1, kcol+nconti;
c                               xbs((i-1)*kcol+nconti+j)=x(i),
c                                    i=2, nint;  j=1, kcol
c                               xbs(ncpts+i)=x(nint+1), i=1,kcol+nconti.
c
c-----------------------------------------------------------------------
c Local Variables:
        double precision        rho(mxkcol+1)
c                               rho stores the Gaussian points.
c
        double precision        wts(mxkcol+1)
c                               wts stores the Gaussian weights.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 ii
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               gauleg
c
c-----------------------------------------------------------------------
c     Generate the piecewise polynomial space breakpoint sequence.
      do 10 i = 1, kcol + nconti
         xbs(i) = x(1)
         xbs(i + ncpts) = x(nint + 1)
   10 continue
      do 30 i = 2, nint
         ii = (i - 2) * kcol + kcol + nconti
         do 20 j = 1, kcol
            xbs(ii + j) = x(i)
   20    continue
   30 continue

c-----------------------------------------------------------------------
c     Compute the Gaussian points.
      call gaulegri(kcol, kcol*kcol, rho, wts, work, 2)

c     Define the collocation point sequence.
      xcol(1) = x(1)
      do 50 i = 1, nint
         ii = (i - 1) * kcol + 1
         do 40 j = 1, kcol
            xcol(ii + j) = x(i) + h(i) * rho(j)
   40    continue
   50 continue
      xcol(ncpts) = x(nint + 1)

      return
      end
      subroutine errestri(kcol, nint, npde, neq, npts, icount,
     &                  xsol, wts, xbs, y, istart, mflag2, atol,
     &                  rtol, lenwk, work, errbas, ercoef, lencof,
     &                  errrat, errint, errcom, ieflag, x, h, est)

c-----------------------------------------------------------------------
c Purpose:
c       This routine computes the error estimate for each subinterval
c       and for each component of the PDE system, and decides whether a
c       remeshing is necessary or not. If a remeshing is deemed
c       necessary, the distribution of mesh points is determined by this
c       error estimate.
c
c-----------------------------------------------------------------------
c
c Last modified by Jack Pew, August 8, 2011.
c
c Explanation of Modifications: This version of errest has replaced the
c evaluation of a second computed solution with interpolated values from
c the (now only) computed solution. According to a choice of est, it will
c either use local extrapolation (the LOI scheme) or superconvergent
c interpolation (the SCI scheme) as the replacement.
c
c------------------------------------------------------------------------
        implicit none
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        integer                 nintsm
        parameter              (nintsm = 15)
c                               when the current step is the first step
c                               after remeshing, we require
c                               if nint <= nintsm
c                                    errrat < saffa2
c                               else
c                                    saffa1 < errrat < saffa2.
c                               endif
c
        double precision        zero
        parameter              (zero = 0.0d0)
c
        double precision        one
        parameter              (one = 1.0d0)
c
        double precision        two
        parameter              (two = 2.0d0)
c
        double precision        saffa1
        parameter              (saffa1 = 0.1d0)
c
        double precision        saffa2
        parameter              (saffa2 = 0.4d0)
c                               These safety factors are heuristic
c                               values that bound the acceptable global
c                               scaled error estimate after a remeshing.
c                               They are used in Spatial Error Test II.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval in the
c                               computed solution.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 neq
c                               neq=npde*(nint*kcol+nconti) is the
c                               number of bspline coefficients (or
c                               DAEs) when using dassl_kcol.
c
        integer                 npts
c                               npts is the number of points in the
c                               x vector, which is equal to
c                               nint*quad.
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        double precision        xsol(npts)
c                               xsol is the npts Gauss-Legendre
c                               points at which the solution are
c                               to be calculated for the L2-norm
c                               of the error.
c
        double precision        wts(npts)
c                               wts is the npts Gauss-Legendre
c                               weights at the corresponding xsol.
c
        double precision        x(nint+1)
c                               The current mesh.
c
        double precision        h(nint)
c                               The mesh step size sequence.
c
        double precision        xbs((kcol+1)*nint+nconti+nconti)
c                               xbs is the breakpoint sequence when
c                               using dassl_kcol.
c
        double precision        y(neq)
c                               y is the vector of bspline
c                               coefficients when using dassl_kcol.
c
        integer                 istart
c                               istart is a flag to begin the code.
c                               istart = 0, it is the initial step;
c                                      = 1, it is not the initial step.
c
        integer                 mflag2
c                               mflag2 = 0, scalar atol and rtol.;
c                               mflag2 = 1, vector atol and rtol.
c
        integer                 est
c                               est=0 => use LOI error estimate
c                               est=1 => use SCI error estimate

        double precision        atol(npde)
c                               atol is the absolute error tolerance
c                               request and is a scalar quantity if
c                               mflag2 = 0.
c
        double precision        rtol(npde)
c                               rtol is the relative error tolerance
c                               request and is a scalar quantity if
c                               mflag2 = 0.
c
        integer                 lenwk
c                               lenwk is the size of the work storage
c                               array and must satisfy:
c                               lenwk >= 2*npde*nint*(kcol+2)
c                                        +npde*nint
c                                        +(kcol-3)*nint*npde
c                                        +2*(nint+1)*npde
c                               if using the LOI scheme, or:
c                               lenwk >= 2*npde*nint*(kcol+3)
c                                        +npde*nint
c                                        +(kcol-2)*nint*npde
c                                        +2*(nint+1)*npde
c                               if using the SCI scheme.
c
        integer                 lencof
c                               lencof is the size of the coefficient
c                               storage array and must satisfy:
c                               lencof >= (kcol+1)*(kcol+2)
c                                         +(kcol+nconti)*(kcol-3)*nint
c                                         +(kcol+nconti)*2*(nint+1)
c                               if using the LOI scheme, or:
c                               lencof >= (kcol+4)*(kcol+3)*nint
c                                         +(kcol+nconti)*(kcol-2)*nint
c                                         +(kcol+nconti)*2*(nint+1)
c                               if using the SCI scheme.
c
c       Work Storage:
        double precision        work(lenwk)
c                               work is a floating point work storage
c                               array of size lenwk.
c
c       output:
        double precision        errbas((kcol+nconti)*npts)
c                               errbas holds the values of the nonzero
c                               basis functions at xsol when using
c                               dassl_kcol.
c                               They are reusable between remeshings.
c
        double precision        ercoef(lencof)
c                               ercoef holds the values of the nonzero
c                               B-spline basis functions and the
c                               Hermite-Birkhoff coefficients used in
c                               the selected interpolation subroutine.
c                               They are reusable provide the spatial
c                               mesh does not change.
c
        double precision        errrat
c                               errrat is the value of the largest
c                               component of errcom.
c
        double precision        errint(nint)
c                               errint is the error estimate at each
c                               subinterval.
c
        double precision        errcom(npde)
c                               errcom is the error estimate for
c                               each component of pdes at the whole
c                               range, i.e. from x_a to x_b.
c
        integer                 ieflag
c                               ieflag is a status flag for remesh.
c                               ieflag = 0, indicates no need remeshing.
c                               ieflag = 1, indicates need remeshing.
c
C-----------------------------------------------------------------------
c Local Variables:
        double precision        errsum
c                               errsum is the sum of errint.
c
        double precision        errmax
c                               errmax is the maximum value of
c                               errint(i), i = 1, nint.
c
        double precision        aerr
c                               aerr is the average value of errint(i),
c                               i = 1, nint.
c
        double precision        disind
c                               disind is equal to errmax/aerr, and it
c                               indicates the error distribution over
c                               the mesh.
c
        integer                 quad
c                               quad is the number of Gaussian
c                               quadrature points used in the L2-norm.
c
        double precision        power
c                               The exponent applied to each errint(i)
c                               component. This is one/dble(kcol+2) for
c                               the SCI and one/dble(kcol+1) for the LOI
c                               scheme.
c
        double precision        rL, rR
c                               Mesh subinterval size ratios
c                              (for scaling the SCI estimate.)
c
c       Pointers into the floating point work arrays:
        integer                 iusol1
c                               work(iusol1) stores the evaluations of
c                               the collocation solution at the npts
c                               Gaussian quadrature points.
c
        integer                 iusol2
c                               work(iusol2) stores the evaluations of
c                               interpolant at the npts Gaussian
c                               quadrature points.
c                               Under the LOI scheme, they are one
c                               order of accuracy lower than those in
c                               work(iusol1), and under the SCI scheme
c                               they are one order higher.
c
        integer                 ierrci
c                               work(ierrci) stores the error estimate
c                               for each subinterval for each component.
c
        integer                 iliu
c                               work(iliu) stores evaluations of the
c                               computed solution at the superconvergent
c                               points internal to each subinterval.
c
        integer                 ilium
c                               work(ilium) stores evaluations of the
c                               computed solution at the mesh points.
c
        integer                 ih
c                               ercoef(ih) stores the H coefficients
c                               for the Hermite-Birkhoff interpolant.
c
        integer                 ihd
c                               ercoef(ihd) stores the Hbar coefficients
c                               for the Hermite-Birkhoff interpolant.
c
        integer                 ig
c                               ercoef(ig) stores the G coefficients
c                               for the Hermite-Birkhoff interpolant.
c
        integer                 ibassc
c                               ercoef(ibassc) stores the evaluations
c                               of the B-spline basis functions at
c                               the superconvergent points internal
c                               to each subinterval.
c
        integer                 ibasm
c                               ercoef(ibasm) stores the evaluations of
c                               the B-spline basis functions at the
c                               mesh points. It and ercoef(ibassc) are
c                               analogs to the errbas vector, which
c                               stores basis function evaluations at
c                               quadrature points for the subroutine
c                               errval.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i, j, m, ij, im, mm
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               errval
c                               scint
c                               lowint
c
c-----------------------------------------------------------------------
c     quad is the number of quadrature points used per subinterval
      if (est .eq. 0) then
         quad = kcol + 2
      elseif (est .eq. 1) then
         quad = kcol + 3
      endif

c     Set pointers into two work arrays known to errest, work and ercoef
c     The SCI scheme requires more storage here, as its H-B coefficients
c     are different from one subinterval to the next due to the
c     dependence on adjacent subinterval sizes.
c     Also, the SCI uses a greater number of quadrature points.
      iusol1 = 1
      iusol2 = iusol1 + npde * nint * quad
      ierrci = iusol2 + npde * nint * quad
      iliu   = ierrci + npde * nint
      if (est .eq. 0) then
         ilium  = iliu   + npde * nint * (kcol - 3)
         ih     = 1
         ihd    = ih     + 2 * quad
         ig     = ihd    + 2 * quad
         ibassc = ig     + (kcol - 3) * quad
         ibasm  = ibassc + (kcol + nconti) * (kcol - 3) * nint
      elseif (est .eq. 1) then
         ilium  = iliu   + npde * nint * (kcol - 2)
         ih     = 1
         ihd    = ih     + 2 * quad * nint
         ig     = ihd    + 2 * quad * nint
         ibassc = ig     + kcol * quad * nint
         ibasm  = ibassc + (kcol + nconti) * (kcol - 2) * nint
      endif

c-----------------------------------------------------------------------
c     Generate the values of the collocation solution and the
c     interpolant at the Gaussian quadrature points, stored in xsol, and
c     save them in work(iusol1) and work(iusol2), respectively.
c     The LOI and the SCI schemes use quadrature rules of different
c     degrees.
c     errval evaluates the computed solution, while lowint and scint
c     evaluate the interpolant.

      call errvalri(kcol, nint, npde, neq, quad, istart, icount,
     &            xbs, xsol, y, errbas, work(iusol1))

      if (est .eq. 0) then
         call lowintri(kcol, nint, npde, istart, icount, xbs, x, y,
     &         ercoef(ih), ercoef(ihd), ercoef(ig), ercoef(ibassc),
     &         ercoef(ibasm), work(iliu), work(ilium), work(iusol2))
      elseif (est .eq. 1) then
         call scintri(kcol, nint, npde, istart, icount, xbs, x, y,
     &         ercoef(ih), ercoef(ihd), ercoef(ig), ercoef(ibassc),
     &         ercoef(ibasm), work(iliu), work(ilium), work(iusol2))
      endif

c-----------------------------------------------------------------------
c     Initialization task.
      do 10 i = 1, nint
         errint(i) = zero
   10 continue

      do 20 i = 1, npde
         errcom(i) = zero
   20 continue

      do 30 i = 1, npde * nint
         work(ierrci - 1 + i) = zero
   30 continue

c-----------------------------------------------------------------------
c     Calculate the error estimate at each subinterval for each
c     component of PDEs.

      if (mflag2 .eq. 0) then
c        Use scalar error tolerance.
         do 60 m = 1, npde
            do 50 i = 1, nint
               do 40 j = 1, quad
                  ij = (i - 1) * quad + j
                  mm = npde * (ij - 1) + m
                  im = ierrci - 1 + (m - 1) * nint + i
                  work(im) = work(im)
     &                     + ((work(iusol1-1+mm) - work(iusol2-1+mm))
     &                     / (atol(1) + rtol(1)*abs(work(iusol1-1+mm))))
     &                     **2 * wts(ij)
   40          continue
   50       continue
   60    continue
      else
c        Use vector error tolerance (to weight PDEs in the system.)
         do 90 m = 1, npde
            do 80 i = 1, nint
               do 70 j = 1, quad
                  ij = (i - 1) * quad + j
                  mm = npde * (ij - 1) + m
                  im = ierrci - 1 + (m - 1) * nint + i
                  work(im) = work(im)
     &                     + ((work(iusol1-1+mm) - work(iusol2-1+mm))
     &                     / (atol(m) + rtol(m)*abs(work(iusol1-1+mm))))
     &                     **2 * wts(ij)
   70          continue
   80       continue
   90    continue
      endif

c-----------------------------------------------------------------------
c     Scale work(ierrci) to `correct' for overestimations by the
c     SCI due to mesh ratios. We divide by the product of the ratios
c     instead of the only the larger one because for internal
c     subintervals, the underestimation of errors in layer regions is an
c     issue. (This is also the case in the original BACOLI.)
      if (EST .eq. 1) then
         do j = 1, npde
c           i = 1
            rL = zero
            rR = h(2) / h(1)
            ij = ierrci + (j - 1) * nint
            work(ij) = work(ij) * rR * rR
            do i = 2, nint-1
               rL = h(i-1) / h(i)
               rR = h(i+1) / h(i)
               ij = ij + 1
               work(ij) = work(ij) * rL * rR
            end do
c           i = nint
            rL = h(nint-1) / h(nint)
            rR = zero
            ij = ij + 1
            work(ij) = work(ij) * rL * rL
         end do
      end if

c     Calculate errint and errcom.
      do 110 j = 1, npde
         do 100 i = 1, nint
            ij = ierrci - 1 + (j - 1) * nint + i
            errint(i) = errint(i) + work(ij)
            errcom(j) = errcom(j) + work(ij)
  100    continue
  110 continue

c     When using the LOI scheme, error is actually estimated for a
c     solution that is one order lower than that which is computed.
      if (est .eq. 0) power = one/dble(kcol+1)
      if (est .eq. 1) power = one/dble(kcol+2)

c     Take the square root and update errint and errcom.
      do 120 i = 1, nint
         errint(i) = sqrt(errint(i))
         errint(i) = errint(i) ** power
  120 continue

      do 130 i = 1, npde
         errcom(i) = sqrt(errcom(i))
  130 continue

c-----------------------------------------------------------------------
c     Decide whether remeshing is needed.
      ieflag = 0

c     update errrat (the max of errcom.)
      errrat = zero
      do 140 i = 1, npde
         if (errcom(i) .gt. errrat) then
            errrat = errcom(i)
         endif
  140 continue

c     Calculate errsum to be the sum of the errint. Find the maximum
c     errint(i) and save it in errmax.
      errsum = errint(1)
      errmax = errint(1)
      do 150 i = 2, nint
         if (errmax .lt. errint(i)) errmax = errint(i)
         errsum = errint(i) + errsum
  150 continue

c     Let aerr be the mean value of errint(i).
      aerr = errsum/dble(nint)

c     Calculate disind (which is a measure of error distribution.)
      disind = errmax/aerr

c     This is Spatial Error Test I.
      if (disind .gt. two) then
c        Mesh not well distributed. Remeshing needed.
         ieflag = 1
      else
c        Passed Test I, now do Test II.
         if ((istart .eq. 1) .and. (icount .eq. 0)) then
            if (errrat .ge. one) ieflag = 1
         else
            if (nint .gt. nintsm) then
               if ((errrat .ge. saffa2) .or. (errrat .le. saffa1))
     &            ieflag = 1
            else
               if (errrat .ge. saffa2) ieflag = 1
            endif
         endif
      endif

      return
      end
      subroutine errvalri(kcol, nint, npde, neq, nptse, istart, icount,
     &                  xbs, xsol, y, errbas, usol)

c-----------------------------------------------------------------------
c Purpose:
c       This routine computes the values of the (kcol+nconti) nonzero
c       bspline basis function at each Gaussian point of xsol.
c       Then determine the solution usol, which is used for error
c       estimate, at xsol.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 29, 2001.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 neq
c                               neq=npde*(kcol*nint+2) is the number of
c                               bspline coefficients.
c
        integer                 nptse
c                               nptse is the number of Gaussian points
c                               in each subinterval for the error
c                               estimate.
c
        integer                 istart
c                               istart is a flag to begin the code.
c                               istart = 0, it is the initial step;
c                                      = 1, it is not the initial step.
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        double precision        xbs((kcol+1)*nint+nconti+nconti)
c                               The breakpoint sequence.
c
        double precision        xsol(nptse*nint)
c                               xsol is a set of spatial points at which
c                               the solution are to be calculated for
c                               error estimate.
c
        double precision        y(neq)
c                               y is the vector of bspline coefficients.
c
c       output:
        double precision        errbas(kcol+nconti, nptse*nint)
c                               errbas is the values of the nonzero
c                               basis functions at xsol.
c
        double precision        usol(npde, nptse*nint)
c                               uval is the solution at xsol.
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 ileft
c                               breakpoint information.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 k
        integer                 m
        integer                 jj
        integer                 mm
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bsplvd
c
c-----------------------------------------------------------------------

c     check whether errbas is necessary to be calculated.
      if ((istart .eq. 1) .and. (icount .eq. 0)) goto 30

c     calculate errbas.
      do 20 i = 1, nint
         ileft = kcol + nconti + (i - 1) * kcol
         do 10 j = 1, nptse
            jj = (i - 1) * nptse + j
            call bsplvdri(xbs, kcol+nconti, xsol(jj), ileft,
     &          errbas(1,jj),
     &                  1)
   10    continue
   20 continue

   30 continue

c     compute the values of usol at xsol.
      do 70 i = 1, nint
         do 60 j = 1, nptse
            jj = (i - 1) * nptse + j
            do 50 k = 1, npde
               usol(k,jj) = zero
               do 40 m = 1, kcol + nconti
                  mm = npde * (m + (i - 1) * kcol - 1) + k
                  usol(k,jj) = usol(k,jj) + y(mm) * errbas(m,jj)
   40          continue
   50       continue
   60    continue
   70 continue

      return
      end
      subroutine evalri(npde,kcol,ileft,icpt,ncpts,uval,uxval,uxxval,
     &                fbasis,y)

c-----------------------------------------------------------------------
c Purpose:
c       This routine evaluates u(k), ux(k), and uxx(k), k=1 to npde,
c       at the icpt-th collocation point.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, Feb. 11, 2001.
c
c-----------------------------------------------------------------------
c Constants
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 ileft
c                               breakpoint information.
c
        integer                 icpt
c                               the index of the collocation point.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        double precision        fbasis((kcol+nconti)*3)
c                               Basis function values at the icpt-th
c                               collocation point.
c                               fbasis(k+(j-1)*(kcol+nconti)) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti).
c
        double precision        y(ncpts*npde)
c                               y is the vector of bspline coefficients.
c
c       Output:
        double precision        uval(npde)
c                               uval gives the approximation to
c                               u(t,x).
c
        double precision        uxval(npde)
c                               uxval gives the approximation to
c                               the first spatial derivative of u(t,x).
c
        double precision        uxxval(npde)
c                               uxxval gives the approximation to
c                               the second spatial derivative of u(t,x).
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 j
        integer                 m
c-----------------------------------------------------------------------
      do 10 j = 1, npde
         uval(j)   = zero
         uxval(j)  = zero
         uxxval(j) = zero
   10 continue
      if (icpt .ne. 1 .and. icpt .ne. ncpts) then
         do 30 j = 1, npde
            do 20 m = 1, kcol + nconti
               uval(j)   = uval(j) + fbasis(m) *
     &                     y((ileft-kcol-3+m) * npde + j)
               uxval(j)  = uxval(j) + fbasis(m+kcol+nconti) *
     &                     y((ileft-kcol-3+m) * npde + j)
               uxxval(j) = uxxval(j) + fbasis(m+2*(kcol+nconti)) *
     &                     y((ileft-kcol-3+m) * npde + j)
   20          continue
   30       continue
      else
         if (icpt .eq. 1) then
            do 40 j = 1, npde
               uval(j)   = uval(j)   + fbasis(1) * y(j)
               uxval(j)  = uxval(j)  + fbasis(1+kcol+nconti) * y(j)
     &                     + fbasis(2+kcol+nconti) * y(npde + j)
               uxxval(j) = uxxval(j) + fbasis(1+2*(kcol+nconti)) * y(j)
     &                     + fbasis(2+2*(kcol+nconti)) * y(npde + j)
     &                     + fbasis(3+2*(kcol+nconti)) * y(2*npde + j)
   40       continue
         else
            do 50 j = 1, npde
               uval(j)   = uval(j)   + fbasis(kcol+nconti)
     &                     * y((ncpts - 1) * npde + j)
               uxval(j)  = uxval(j)  + fbasis((kcol+nconti)*2)
     &                     * y((ncpts - 1) * npde + j)
     &                     + fbasis((kcol+nconti)*2-1)
     &                     * y((ncpts - 2) * npde + j)
               uxxval(j) = uxxval(j) + fbasis((kcol+nconti)*3)
     &                     * y((ncpts - 1) * npde + j)
     &                     + fbasis((kcol+nconti)*3-1)
     &                     * y((ncpts - 2) * npde + j)
     &                     + fbasis((kcol+nconti)*3-2)
     &                     * y((ncpts - 3) * npde + j)
   50       continue
         endif
      endif
      return
      end
      subroutine iniyri(t0, npde, kcol, nint, neq, ncpts, ifglin, xcol,
     &                xbs, abdblk, fbasis, y, ipivot, work, lw, icflag,
     &                ifgfdj, bndxa, difbxa, bndxb, difbxb, uinit,
     &                uinitvec)

c-----------------------------------------------------------------------
c Purpose:
c       This routine performs the initialization tasks required by
c       inital including:
c
c               calculating the Bspline basis functions,
c               constructing abdblk of the collocation matrices and
c               determining y(t0).
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, November 8, 2001.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
        double precision        negone
        parameter              (negone = -1.0D0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        double precision        t0
c                               t0 is the initial time.
c
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 neq
c                               neq=npde*(kcol*nint+nconti) is
c                               the number of bsplines
c                               coefficients (or DAEs).
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        integer                 ifglin
c                               ifglin is a flag for the boundary
c                               conditions.
c                               ifglin = 1, indicate both derichlet
c                                           boundary conditions;
c                                      = 0, else.
c
        double precision        xcol(ncpts)
c                               The sequence of collocation points on
c                               the interval [x_a, x_b].
c
        double precision        xbs(ncpts+kcol+nconti)
c                               The breakpoint sequence.
c                               xbs(i)=x(1), i=1, kcol+nconti;
c                               xbs((i-1)*kcol+nconti+j)=x(i),
c                                    i=2, nint;  j=1, kcol
c                               xbs(ncpts+i)=x(nint+1), i=1,kcol+nconti.
c
        integer                 lw
c                               lw is the size of the work storage
c                               array and must satisfy:
c                               lw >= 2*npde*npde*nconti+
c                                     npde*npde*kcol*(kcol+nconti)*nint
c                                     +2*neq+2*npde+2*npde*npde
        integer                 ifgfdj
c                               Are finite difference approximations
c                               being used to approximate Jacobian
c                               matrices?
c                               Set to 0 to approximate derivf, difbxa 
c                               and difbxb.
c                               Set to 1 to approxiate only difbxa and 
c                               difbxb.
c                               Set to 2 to approximate only derivf.
c                               Set to 3 if all of derivf, difbxa and 
c                               difbxb are provided.
c
C     BACOLRI --> BACOLRIVEC
C         external                f
C         external                bndxa
C         external                difbxa
C         external                bndxb
C         external                difbxb
C         external                uinit
        external                bndxa
        external                difbxa
        external                bndxb
        external                difbxb
        external                uinit
        external                uinitvec
c                               refer to the preamble of BACOLI.
c
c       Work Storage:
        integer                 ipivot(neq)
c                               pivoting information from the
c                               factorization of the temporary matrix.
c
        double precision        work(lw)
c                               work is a floating point work storage
c                               array of size lw.
c
c       Output:
        double precision        abdblk(npde*npde*nint*kcol
     &                                 *(kcol+nconti))
c                               The nint blocks in the middle of
c                               the matrix A.
c
        double precision        fbasis(kcol+nconti, 3, ncpts)
c                               Basis function values at the collocation
c                               points. fbasis(k,j,i) contains the
c                               values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        double precision        y(neq)
c                               y = y(t0) is the initial vector of
c                               bspline coefficients.
c
        integer                 icflag
c                               This is the status flag from COLROW
c                               which is called by crdcmp.
c                               icflag =  0, indicates non-singularity.
c                               icflag = -1, indicates singularity.
c                               icflag =  1, indicates invalid input.
c-----------------------------------------------------------------------
c Local Variables:
        integer                 ileft
c                               breakpoint information.
c
        integer                 nels
c                               the number of elements in one
c                               collocation block of work.
c
c       Pointers into the floating point work array:
        integer                 iabdtp
c                               work(iabdtp) contains a copy of the top
c                               block which is required since crdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbk
c                               work(iabdbk) contains a copy of abdblk
c                               which is required since crdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbt
c                               work(iabdbt) contains a copy of the
c                               bottom block which is required since
c                               crdcmp overwrites the input collocation
c                               matrix.
c
        integer                 idelta
c                               work(idelta) contains the residual which
c                               indicates how well y satisfies to the
c                               boundary condition and the initial
c                               condition at the internal collocation
c                               points.
c
        integer                 ivcol
c                               work(ivcol) contains the values of u
c                               at the internal collocation points.
c
        integer                 iu
c                               work(iu) stores the approximation to
c                               u(t,x).
c
        integer                 iux
c                               work(iux) stores the approximation to
c                               the first spatial derivative of u(t,x).
c
        integer                 iuxx
c                               work(iuxx) stores the approximation to
c                               the second spatial derivative of u(t,x).
c
        integer                 idbdu
c                               work(idbdu-1+i), i=1, npde*npde,
c                               contains dbdu(npde,npde). That is,
c                               dbdu(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the unknown function u.
c
        integer                 idbdux
c                               work(idbdux-1+i), i=1, npde*npde,
c                               contains dbdux(npde,npde). That is,
c                               dbdux(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the spatial derivative of the
c                               unknown function u.
c
        integer                 idbdt
c                               work(idbdt-1+i), i=1, npde, contains
c                               the partial derivative of the i-th
c                               component of the vector b with respect
c                               to time t.
c
        integer                 ifdwrk
c                               work(ifdwrk-1+i), i=1, 2*npde, is used
c                               as work storage fdderivf and fdbndx.
        double precision        x(neq)
c                               BACOLRI --> BACOLRIPY
c                               Array for copying RHS of linear system
c                               to be solved with LAMSOL. This is done
c                               since LAMSOL uses different arrays for
c                               input and output. This will be used for
c                               input.
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 l
        integer                 m
        integer                 ii
        integer                 jj
        integer                 ll
        integer                 mm
c
c     BACOLRI --> BACOLRIVEC
        integer                 vecoff
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bndxa
c                               bndxb
c                               fdbndx
c                               bsplvd
c                               difbxa
c                               difbxb
c                               eval
c                               uinit
c                               crdcmp
c                               crslve
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c       double precision:
c                               dcopy
c                               dscal
c
c-----------------------------------------------------------------------
      nels = npde*npde*kcol*(kcol+nconti)

c     Set the pointers into the floating point work array.
      iabdtp = 1
      iabdbk = iabdtp + npde*npde*nconti
      iabdbt = iabdbk + nint*nels
      idelta = iabdbt + npde*npde*nconti
      ivcol  = idelta + neq
      iu     = ivcol  + neq-2*npde
      iux    = iu     + npde
      iuxx   = iux    + npde
      idbdu  = iuxx   + npde
      idbdux = idbdu  + npde*npde
      idbdt  = idbdux + npde*npde
      ifdwrk = idbdt  + npde
c     BACOLRI --> BACOLRIVEC
      vecoff = ifdwrk + 2*npde

c-----------------------------------------------------------------------
c     Initialize abdblk, the top block and the bottom block to zero.
      do 20 i = 1, npde*npde*nconti
         work(iabdtp-1+i) = zero
         work(iabdbt-1+i) = zero
   20 continue
      do 30 i = 1, nint*nels
         abdblk(i) = zero
   30 continue

c-----------------------------------------------------------------------
c     Bsplvd is called to compute the components of fbasis(k,i,j)
c     associated the first collocation point. Now ileft = kcol + nconti.
      call bsplvdri(xbs,kcol+nconti,xcol(1),kcol+nconti,fbasis(1,1,1),3)

c     Uinit is called to evaluate the first npde components at the
c     left boundary point, and save in y.
c     BACOLRI --> BACOLRIVEC
      call uinit(xcol(1), y(1), npde)

c     Makeing use of the fact that only the first bspline has a nonzero
c     value at the left end point, set up the top block in work.
      do 40 i = 1, npde
         ii = (i-1) * npde + i
         work(iabdtp-1+ii) = fbasis(1,1,1)
   40 continue

c-----------------------------------------------------------------------
c     The nint blocks at the middle of the matrix will now be set up.
      do 80 i = 1, nint

c     Make use the fact that there are kcol collocation points in each
c     subinterval to find the value of ileft.
         ileft = kcol + nconti + (i - 1) * kcol

         do 70 j = 1, kcol

c     ii is the position in xcol of the j-th collocation point of the
c     i-th subinterval.
            ii = (i-1) * kcol + 1 + j

c     jj is the position in the y vector where the values for the
c     right hand side of the initial conditions, evaluated at the ii-th
c     collocation point are stored.
            jj = (ii - 1) * npde + 1

c     compute information for ii-th collocation point.
            call bsplvdri(xbs,kcol+nconti,xcol(ii),ileft,
     &          fbasis(1,1,ii),3)
c     BACOLRI --> BACOLRIVEC
            call uinit(xcol(ii), y(jj), npde)

            do 60 l = 1, kcol + nconti

c     generate the subblock in abdblk corresponding to the ii-th
c     collocation point.
c
               ll = (l-1)*npde*npde*kcol + (i-1)*nels + (j-1)*npde
               do 50 m = 1, npde
                  mm = ll + (m -1)*npde*kcol + m
                  abdblk(mm) = fbasis(l,1,ii)
   50          continue
   60       continue
   70    continue
   80 continue

c-----------------------------------------------------------------------
c     Now, set up the bottom block, using the fact that only the
c     last bspline basis function is non-zero at the right end point.
c     Simultaneously, set up the corresponding part of the right hand
c     side.
c
      call bsplvdri(xbs,kcol+nconti,xcol(ncpts),ncpts,
     &            fbasis(1,1,ncpts),3)
C     BACOLRI --> BACOLRIVEC
      ii = neq - npde + 1
      call uinit(xcol(ncpts), y(ii), npde)
C       call uinitvec(xcol(1), y(1), ncpts, npde)

C     Copy all of the points back into place.
C       call dcopy(neq, y(1), 1, work(vecoff), 1)

C     Put all of the points back into place
C       do 201 i = 1, ncpts
C         do 211 j = 1, npde
C           y(npde*(i-1)+j) = work(vecoff+(j-1)*neq/npde+i-1)
C   211   continue
C   201 continue

      do 90 i = 1, npde
         ii = ((i-1)+npde)*npde + i
         work(iabdbt-1+ii) = fbasis(kcol+nconti,1,ncpts)
   90 continue

c-----------------------------------------------------------------------
c     Copy the middle of the collocation matrix into temporary storage.
      call dcopyri(nint*nels,abdblk,1,work(iabdbk),1)

c     Check whether both boundary conditions are derichlet boundary
c     conditions. If no, copy the values at the internal collocation
c     points to work(ivcol), which will be used for newton iterations.
      if (ifglin .eq. 0) then
         call dcopyri(neq-2*npde,y(npde+1),1,work(ivcol),1)
         call dscalri(neq-2*npde,negone,work(ivcol),1)
      endif

c-----------------------------------------------------------------------
c     Generate the initial vector y(t0).
c-----------------------------------------------------------------------

c     LU decompose the matrix.
c     BACOLRI --> BACOLRILAM
       call lamdecri(neq,work(iabdtp),npde,2*npde,work(iabdbk),
     &             kcol*npde,
     &             (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &             icflag)

      if (icflag .ne. 0) goto 999

c     Solve the linear system. If derichlet boundary conditions are
c     given, this gives the basis function coefficients for the initial
c     conditions, i.e. y(t0). If not, this gives the predictor of y(t0).
c     BACOLRI --> BACOLRILAM
      call dcopyri(neq, y, 1, x, 1)

      call lamsolri(work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            x,y)

      if (icflag .ne. 0) goto 999

c     Check whether both boundary conditions are derichlet boundary
c     conditions.
      if (ifglin .eq. 1) goto 999

c-----------------------------------------------------------------------
c     Newton iteration loop.

c     Calculate (work(idelta-1+i), i = npde+1, neq-npde), which depends
c     on the nint blocks in the middle of the collocation matrix A.
      call dcopyri(neq-2*npde,work(ivcol),1,work(idelta+npde),1)
      do 130 i = 1, nint
         do 120 j = 1, kcol + nconti
            do 110 l = 1, kcol
               ll = 1+(i-1)*npde*npde*kcol*(kcol+nconti)
     &              +(j-1)*npde*npde*kcol+(l-1)*npde
               do 100 m = 1, npde
                  ii = idelta-1+npde+(i-1)*npde*kcol+(l-1)*npde+m
                  mm = (i-1)*kcol*npde+(j-1)*npde+m
                  work(ii) = work(ii) + abdblk(ll) * y(mm)
  100          continue
  110       continue
  120    continue
  130 continue

c     Copy the middle of the collocation matrix into temporary storage.
      call dcopyri(nint*nels,abdblk,1,work(iabdbk),1)

c     Update the values at the left boundary.
      call evalri(npde,kcol,kcol+2,1,ncpts,work(iu),work(iux),
     &          work(iuxx),fbasis(1,1,1),y)
      call bndxa(t0, work(iu), work(iux), work(idelta), npde)
      if (ifgfdj .lt. 2) then
          call fdbndxri(t0, work(iu), work(iux), work(idbdu),
     &                work(idbdux), work(idbdt), npde,
     &                bndxa, work(ifdwrk))
      else
          call difbxa(t0, work(iu), work(iux), work(idbdu),
     &                work(idbdux), work(idbdt), npde)
      end if

c     Set up the top block and save in work(iabdtp).
      do 150 j = 1, npde
         do 140 i = 1, npde
            ii = iabdtp - 1 + (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            work(jj) = fbasis(2,2,1) * work(idbdux-1+mm)
            work(ii) = work(idbdu-1+mm) - work(jj)
  140    continue
  150 continue

c     Update the values at the right boundary.
      call evalri(npde,kcol,ncpts,ncpts,ncpts,work(iu),work(iux),
     &          work(iuxx),fbasis(1,1,ncpts),y)
      call bndxb(t0, work(iu), work(iux), work(idelta+neq-npde), npde)
      if (ifgfdj .lt. 2) then
          call fdbndxri(t0, work(iu), work(iux), work(idbdu),
     &                work(idbdux), work(idbdt), npde,
     &                bndxb, work(ifdwrk))
      else
          call difbxb(t0, work(iu), work(iux), work(idbdu),
     &                work(idbdux), work(idbdt), npde)
      end if

c     Set up the bottom block and save in work(iabdbt).
      do 170 j = 1, npde
         do 160 i = 1, npde
            ii = iabdbt - 1 + (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            work(ii) = fbasis(kcol+1,2,ncpts) * work(idbdux-1+mm)
            work(jj) = work(idbdu-1+mm) - work(ii)
  160    continue
  170 continue

c     LU decompose the matrix.

c     BACOLRI --> BACOLRILAM
      call lamdecri(neq,work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            icflag)

      if (icflag .ne. 0) goto 999

c     Solve the corrector equation.
c     BACOLRI --> BACOLRILAM
      x = work(idelta:idelta+neq-1)
      call dcopyri(neq, work(idelta), 1, x, 1)
      call lamsolri(work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            x,work(idelta))

      if (icflag .ne. 0) goto 999

c     Now generate the corrector of y(t0).
      do 180 i = 1, neq
         y(i) = y(i) - work(idelta-1+i)
  180 continue

c-----------------------------------------------------------------------

  999 return
      end
      subroutine meshsqri(kcol, nint, x, work, h, excol, ewts)

c-----------------------------------------------------------------------
c Purpose:
c       This routine calculates the mesh size sequence, then generates
c       the collocation points and Gaussian weights for error
c       estimate.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, April 5, 2001.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 mxkcol
        parameter              (mxkcol = 10)
c                               mxkcol is the maximum number of
c                               collocation points per subinterval.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        double precision        x(nint+1)
c                               x is the spatial mesh which divides the
c                               interval [x_a, x_b] as: x_a = x(1) <
c                               x(2) < x(3) < ... < x(nint+1) = x_b.
c
c       Work Storage:
        double precision        work((kcol+3)*(kcol+3))
c                               work is a floating point work storage
c                               array.
c
c       Output:
        double precision        h(nint)
c                               h is the mesh step size sequence.
c
        double precision        excol(nint*(kcol+3))
c                               excol is the collocation point sequence
c                               which is used for error estimate.
c
        double precision        ewts(nint*(kcol+3))
c                               ewts is the gaussian weight sequence
c                               which is used for error estimate.
c
c-----------------------------------------------------------------------
c Local Variables:
        double precision        rho(mxkcol+3)
c                               rho stores the Gaussian points.
c
        double precision        wts(mxkcol+3)
c                               wts stores the Gaussian weights.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 ii
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               gauleg
c
c-----------------------------------------------------------------------
c     Calculate the mesh step size sequence.
      do 10 i = 1, nint
         h(i) = x(i+1)-x(i)
   10 continue

c     Compute the Gaussian points and Gaussian weights.
      call gaulegri(kcol+3, (kcol+3)*(kcol+3), rho, wts,
     &            work, 4)

c     Define the collocation point sequence.
      do 30 i = 1, nint
         ii = (i - 1) * (kcol + 3)
         do 20 j = 1, kcol+3
            excol(ii + j) = x(i) + h(i) * rho(j)
            ewts(ii + j) = h(i) * wts(j)
   20    continue
   30 continue

      return
      end
      subroutine radfcn(neq, t, y, fr, rpar, ipar, f, fvec, derivf,
     &                  bndxa, bndxb)

c-----------------------------------------------------------------------
c Purpose:
c       This routine calculates the vector at the right side of the
c       DAEs, which is required by RADAU5.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, May 21, 2003.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 neq
c                               neq is the number of DAEs.

        external                f
        external                fvec
        external                derivf
        external                bndxa
        external                bndxb

c
        double precision        t
c                               t is the current time.
c
        double precision        y(neq)
c                               y is the vector of bspline coefficients,
c                               including the one for radau_kcol and the
c                               one for radau_kcol+1
c
        double precision        rpar(*)
c                               rpar is the BACOLRI floating point work
c                               array.
c
        integer                 ipar(*)
c                               ipar is the BACOLRI integer work array.
c
c       Output:
        double precision        fr(neq)
c                               fr is the vector at the right side of
c                               the DAEs.
c
c-----------------------------------------------------------------------
c Direct pointers into the IPAR integer work array:
        integer                 inpde
        parameter              (inpde  =  1)
c                               ipar(inpde) = npde
c
        integer                 ikcol
        parameter              (ikcol  =  2)
c                               ipar(ikcol) = kcol.
c
        integer                 inint
        parameter              (inint  =  3)
c                               ipar(inint) = nint.
c
        integer                 incpt
        parameter              (incpt =  4)
c                               ipar(incpt) = ncpts.
c
c-----------------------------------------------------------------------
c Indirect pointers into the RPAR floating point work array:
        integer                 ixcol
        parameter              (ixcol = 22)
c                               rpar(ipar(ixcol)) stores the
c                               collocation points when using
c                               radau_kcol.
c
        integer                 iwkrj
        parameter              (iwkrj  = 30)
c                               rpar(ipar(iwkrj)) stores an additional
c                               work array required by res and jac.
c
        integer                 ibasi
        parameter              (ibasi = 31)
c                               rpar(ipar(ibasi)) stores the basis
c                               function values at the collocation
c                               points when using radau_kcol.
c                               rpar(ipar(ibasi)) contains
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts). A(k,j,i) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c-----------------------------------------------------------------------
c Local variables:
        integer                 npde
        integer                 kcol
        integer                 nint
        integer                 ncpts
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               calfcn
c
c-----------------------------------------------------------------------

      npde   = ipar(inpde)
      kcol   = ipar(ikcol)
      nint   = ipar(inint)
      ncpts = ipar(incpt)

c     Calculate the right side of the DAEs for radau_kcol.
      call calfcn(npde, kcol, nint, ncpts, neq, rpar(ipar(ixcol)),
     &            rpar(ipar(ibasi)), t, y, rpar(ipar(iwkrj)), fr,
     &            f, fvec, bndxa, bndxb)

      return
      end
      subroutine radjac(neq, t, y, dfdy, rpar, ipar, f, derivf, bndxa,
     &                  difbxa, bndxb, difbxb)
c-----------------------------------------------------------------------
c Purpose:
c       This routine calculates the Jacobian matrix for the right side
c       of the DAEs, which is required by RADAU5.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, May 22, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 neq
c                               neq is the number of DAEs.

        external                f
        external                derivf
        external                bndxa
        external                difbxa
        external                bndxb
        external                difbxb

c
        double precision        t
c                               t is the current time.
c
        double precision        y(neq)
c                               y is the vector of bspline coefficients,
c                               including the one for radau_kcol and the
c                               one for radau_kcol+1
c
        double precision        rpar(*)
c                               rpar is the BACOLRI floating point work
c                               array.
c
        integer                 ipar(*)
c                               ipar is the BACOLRI integer work array.
c
c       Output:
        double precision        dfdy(*)
c                               dfdy is the ABD Jacobian matrix for the
c                               right side of the DAE.
c
c-----------------------------------------------------------------------
c Direct pointers into the IPAR integer work array:
        integer                 inpde
        parameter              (inpde  =  1)
c                               ipar(inpde) = npde
c
        integer                 ikcol
        parameter              (ikcol  =  2)
c                               ipar(ikcol) = kcol.
c
        integer                 inint
        parameter              (inint  =  3)
c                               ipar(inint) = nint.
c
        integer                 incpt
        parameter              (incpt =  4)
c                               ipar(incpt) = ncpts.
c
        integer                 imflg6
        parameter              (imflg6 = 17)
c                               ipar(imflg6) stores the value of
c                               mflag(6) in order for it to be passed
c                               down to caljac.
c
c-----------------------------------------------------------------------
c Indirect pointers into the RPAR floating point work array:
        integer                 ixcol
        parameter              (ixcol = 22)
c                               rpar(ipar(ixcol)) stores the
c                               collocation points when using
c                               radau_kcol.
c
        integer                 iabtp
        parameter              (iabtp = 26)
c                               rpar(ipar(iabtp)) stores the top block
c                               of the ABD collocation matrices when
c                               using radau_kcol.
c
        integer                 iabbt
        parameter              (iabbt = 28)
c                               rpar(ipar(iabbt)) stores the bottom
c                               block of the ABD collocation matrices
c                               when using radau_kcol.
c
        integer                 iwkrj
        parameter              (iwkrj  = 30)
c                               rpar(ipar(iwkrj)) stores an additional
c                               work array required by res and jac.
c
        integer                 ibasi
        parameter              (ibasi = 31)
c                               rpar(ipar(ibasi)) stores the basis
c                               function values at the collocation
c                               points when using radau_kcol.
c                               rpar(ipar(ibasi)) contains
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts). A(k,j,i) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
c-----------------------------------------------------------------------
c Local variables:
        integer                 npde
        integer                 kcol
        integer                 nint
        integer                 ncpts
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               caljac
c
c-----------------------------------------------------------------------

      npde   = ipar(inpde)
      kcol   = ipar(ikcol)
      nint   = ipar(inint)
      ncpts  = ipar(incpt)

c     Calculate the Jacobian matrix at the right side of the DAEs for
c     radau_kcol.
      call caljacri(npde, kcol, nint, ncpts, neq, rpar(ipar(ixcol)),
     &            rpar(ipar(ibasi)), rpar(ipar(iabtp)),
     &            rpar(ipar(iabbt)), t, y, rpar(ipar(iwkrj)), dfdy,
     &            ipar(imflg6), f, derivf, bndxa, difbxa, bndxb, difbxb)

      return
      end
      subroutine radmas(am, rpar, ipar)

c-----------------------------------------------------------------------
c Purpose:
c       This routine calculates the mass-matrix which is required by
c       the implicit Runge-Kutta DAE solver RADAU5.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, May 17, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        double precision        rpar(*)
c                               rpar is the BACOLRI floating point work
c                               array.
c
        integer                 ipar(*)
c                               ipar is the BACOLRI integer work array.
c
c       Output:
        double precision        am(*)
c                               am is the ABD mass-matrix.
c
c-----------------------------------------------------------------------
c Direct pointers into the IPAR integer work array:
        integer                 inpde
        parameter              (inpde  =  1)
c                               ipar(inpde) = npde
c
        integer                 ikcol
        parameter              (ikcol  =  2)
c                               ipar(ikcol) = kcol.
c
        integer                 inint
        parameter              (inint  =  3)
c                               ipar(inint) = nint.
c
c-----------------------------------------------------------------------
c Indirect pointers into the RPAR floating point work array:
        integer                 iabbk
        parameter              (iabbk = 27)
c                               rpar(ipar(iabbk)) stores the nint
c                               blocks in the middle of the ABD
c                               collocation matrices when using
c                               radau_kcol.
c
c-----------------------------------------------------------------------
c Local variables:
        integer                 npde
        integer                 kcol
        integer                 nint
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               calmas
c
c-----------------------------------------------------------------------

      npde   = ipar(inpde)
      kcol   = ipar(ikcol)
      nint   = ipar(inint)

c     Calculate the mass-matrix for radau_kcol.
      call calmas(npde, kcol, nint, rpar(ipar(iabbk)), am)

      return
      end
      subroutine reinitri(npde, kcol, kold, nint, ninold, ncpts, neq,
     &                  neqold, x, xold, yold, work, lw, ipivot, h,
     &                  xbs, xcol, fbasis, y, abdblk, icflag)

c-----------------------------------------------------------------------
c Purpose:
c       This routine performs the initialization tasks after remeshing:
c
c               calculating the mesh step size sequence,
c               generating the piecewise polynomial space breakpoint
c               sequence,
c               calculating the collocation point sequence,
c               calculating the B-spline basis functions,
c               constructing abdblk of the collocation matrices and
c               calculating the bspline coefficients at the last step.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, July 15, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval after
c                               remeshing.
c
        integer                 kold
c                               kold is the number of collocation points
c                               to be used in each subinterval before
c                               remeshing.
c
        integer                 nint
c                               nint is the number of subintervals after
c                               remeshing.
c
        integer                 ninold
c                               ninold is the number of subintervals
c                               before remeshing.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        integer                 neq
c                               neq=npde*(kcol*nint+nconti) is
c                               the number of bspline
c                               coefficients after remeshing.
c
        integer                 neqold
c                               neqold=npde*(kold*ninold+nconti) is
c                               the number of bspline
c                               coefficients before remeshing.
c
        double precision        x(nint+1)
c                               x is the spatial mesh after remeshing.
c
        double precision        xold(ninold+1)
c                               xold is the spatial mesh before
c                               remeshing.
c
        double precision        yold(neqold)
c                               yold is the vector of bspline
c                               coefficients at the last time step.
c
        double precision        h(nint)
c                               h is the mesh step size sequence.
c
        integer                 lw
c                               lw is the size of the work storage
c                               array and must satisfy:
c                               lw >= 2*npde*npde*nconti+
c                                     npde*npde*kcol*(kcol+nconti)*nint
c                                     +(kold+nconti)+kold*(ninold+1)
c                                     +2*nconti
c                               Since nint >= ninold/2 and kcol >=
c                               kold+1, it implies that lw >= 3*neqold.
c
c       Work Storage:
        double precision        work(lw)
c                               work is a floating point work storage
c                               array of size lw.
c
        integer                 ipivot(neq)
c                               pivoting information from the
c                               factorization of the temporary matrix.
c
c       Output:
        double precision        xcol(ncpts)
c                               The sequence of collocation points on
c                               the interval [x_a, x_b].
c
        double precision        xbs(ncpts+kcol+nconti)
c                               The breakpoint sequence.
c
        double precision        fbasis(kcol+nconti, 3, ncpts)
c                               Basis function values at the collocation
c                               points.
c
        double precision        y(neq)
c                               y is the vector of bspline coefficients
c                               at the last time step after remeshing.
c
        double precision        abdblk(npde*npde*nint*kcol
     &                                 *(kcol+nconti))
c                               The nint blocks in the middle of
c                               the matrix A.
c
        integer                 icflag
c                               This is the status flag from COLROW
c                               which is called by crdcmp.
c                               icflag =  0, indicates non-singularity.
c                               icflag = -1, indicates singularity.
c                               icflag =  1, indicates invalid input.
c-----------------------------------------------------------------------
c Local Variables:
        integer                 ileft
c                               breakpoint information.
c
        integer                 nels
c                               the number of elements in one
c                               collocation block of work.
c
c       Pointers into the floating point work array:
        integer                 iabdtp
c                               work(iabdtp) contains a copy of the top
c                               block which is required since crdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbk
c                               work(iabdbk) contains a copy of abdblk
c                               which is required since crdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbt
c                               work(iabdbt) contains a copy of the
c                               bottom block which is required since
c                               crdcmp overwrites the input collocation
c                               matrix.
c
        integer                 ivwork
c                               work(ivwork) is the work storage
c                               required by values.
c
        double precision        x_in(neq)
c                               Array for copying RHS of linear system
c                               to be solved with LAMSOL. This is done
c                               since LAMSOL uses different arrays for
c                               input and output. This will be used for
c                               input.
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 l
        integer                 m
        integer                 ii
        integer                 ll
        integer                 mm
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bsplvd
c                               colpnt
c                               crdcmp
c                               crslve
c                               values
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c       double precision:
c                               dcopy
c
c-----------------------------------------------------------------------

c     Generate the piecewise polynomial space breakpoint sequence,
c     and calculates the collocation point sequences.
      call colpntri(kcol, nint, ncpts, x, h, work, xcol, xbs)

c-----------------------------------------------------------------------
      nels = npde*npde*kcol*(kcol+nconti)

c     Set the pointers into the floating point work array.
      iabdtp = 1
      iabdbk = iabdtp + npde*npde*nconti
      iabdbt = iabdbk + nint*nels
      ivwork = iabdbt + npde*npde*nconti

c-----------------------------------------------------------------------
c     Call values to calculate the values at xcol and at the last
c     time step. Then save in y.
      call valuesri(kold, xcol, ninold, xold, npde, ncpts, 0, y, yold,
     &            work(ivwork))

c-----------------------------------------------------------------------
c     Initialize abdblk, the top block and the bottom block to zero.
      do 10 i = 1, npde * npde * nconti
         work(iabdtp-1+i) = zero
         work(iabdbt-1+i) = zero
   10 continue
      do 20 i = 1, nint*nels
         abdblk(i) = zero
   20 continue

c-----------------------------------------------------------------------
c     Bsplvd is called to compute the components of fbasis(k,i,j)
c     associated the first collocation point. Now ileft = kcol + nconti.
      call bsplvdri(xbs,kcol+nconti,xcol(1),kcol+nconti,fbasis(1,1,1),3)

c     Makeing use of the fact that only the first bspline has a nonzero
c     value at the left end point, set up the top block in work.
      do 30 i = 1, npde
         ii = (i-1) * npde + i
         work(iabdtp-1+ii) = fbasis(1,1,1)
   30 continue

c-----------------------------------------------------------------------
c     The nint blocks at the middle of the matrix will now be set up.
      do 70 i = 1, nint

c     Make use the fact that there are kcol collocation points in each
c     subinterval to find the value of ileft.
         ileft = kcol + nconti + (i - 1) * kcol

         do 60 j = 1, kcol

c     ii is the position in xcol of the j-th collocation point of the
c     i-th subinterval.
            ii = (i-1) * kcol + 1 + j

c     compute information for ii-th collocation point.
            call bsplvdri(xbs,kcol+nconti,xcol(ii),ileft,
     &           fbasis(1,1,ii),3)

            do 50 l = 1, kcol + nconti

c     generate the subblock in abdblk corresponding to the ii-th
c     collocation point.
c
               ll = (l-1)*npde*npde*kcol + (i-1)*nels + (j-1)*npde
               do 40 m = 1, npde
                  mm = ll + (m-1)*npde*kcol + m
                  abdblk(mm) = fbasis(l,1,ii)
   40          continue
   50       continue
   60    continue
   70 continue

c-----------------------------------------------------------------------
c     Now, set up the bottom block, using the fact that only the
c     last bspline basis function is non-zero at the right end point.
c     Simultaneously, set up the corresponding part of the right hand
c     side.
c
      call bsplvdri(xbs,kcol+nconti,xcol(ncpts),ncpts,
     &            fbasis(1,1,ncpts),3)
      do 80 i = 1, npde
         ii = ((i-1)+npde)*npde + i
         work(iabdbt-1+ii) = fbasis(kcol+nconti,1,ncpts)
   80 continue

c-----------------------------------------------------------------------
c     Copy the middle of the collocation matrix into temporary storage.
      call dcopyri(nels*nint,abdblk,1,work(iabdbk),1)

c-----------------------------------------------------------------------
c     Generate the vector y.

c     LU decompose the matrix.

c     BACOLRI --> BACOLRILAM
      call lamdecri(neq,work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            icflag)

      if (icflag .ne. 0) go to 999

c     Solve the linear system. This gives the basis function
c     coefficients for the initial conditions, i.e. y(t0).
c     BACOLRI --> BACOLRILAM
      call dcopyri(neq, y, 1, x_in, 1)
      call lamsolri(work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,
     &            ipivot,x_in,y)
      if (icflag .ne. 0) go to 999

  999 return
      end
      subroutine remeshri(istart, icount, nintmx, ninold, errrat,
     &                  errint, irshfg, nint, kcol, x, work)

c-----------------------------------------------------------------------
c Purpose:
c       This routine generates a new mesh by equidistributing the error
c       in each subinterval.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, July 11, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        double precision        point5
        parameter              (point5 = 0.5d0)
c
        double precision        one
        parameter              (one    = 1.0d0)
c
        double precision        two
        parameter              (two    = 2.0d0)
c
        double precision        saffac
        parameter              (saffac = 0.2d0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 istart
c                               istart is a flag to begin the code.
c                               istart = 0, it is the initial step;
c                                      = 1, it is not the initial step.
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        integer                 nintmx
c                               the maximal number of subintervals that
c                               the user requires.
c
        integer                 ninold
c                               ninold is the number of subintervals
c                               before remeshing.
c
        double precision        errrat
c                               errrat is the value of the largest
c                               component of rpar(ipar(iercom)).
c
        double precision        errint(ninold)
c                               errint is the error estimate at
c                               each subintervals.
c
c       Output:
        integer                 irshfg
c                               irshfg is a flag for redefining all the
c                               pointers.
c                               irshfg = 0, initial call or continuation
c                                           calls;
c                                      = 1, remesh.
c
c       In-output:
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c                               As input, it is the value before
c                               remeshing; as output, it is the value
c                               after remeshing.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               ninmx >= nint >= 1.
c                               As input, it is the value before
c                               remeshing; as output, it is the value
c                               after remeshing.
c
        double precision        x(nintmx+1)
c                               x is the spatial mesh. As input, it is
c                               the value before remeshing; as output,
c                               it is the value after remeshing.
c
c       Work storage:
        double precision        work(2*ninold+1)
c
c-----------------------------------------------------------------------
c Local Variables:
        double precision        aerr
        double precision        berr
c
c       Pointers into the floating point work array:
        integer                 ierror
c                               work(ierror-1+i) is the L2-norm error
c                               estimate at the first i subintervals.
c
        integer                 ixold
c                               work(ixold) contains a copy of mesh
c                               points before remeshing.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
c
c-----------------------------------------------------------------------
c Functions used:
c                               dble
c                               int
c
c-----------------------------------------------------------------------

c     Set the pointers into the floating point work array.
      ierror = 1
      ixold = ierror + ninold

c-----------------------------------------------------------------------
c     Update icount, irshfg and nint.
      icount = icount + 1

c     If this is the first remesh at the current step which is not the
c     initial step.
      if ((icount .eq. 1) .and. (istart .eq. 1)) then
         irshfg = 1
         goto 20
      endif

c     Update errrat.
      errrat = (errrat/saffac) ** (one/dble(kcol+2))

c     Set the upper bound and lower bound of the ratio of nint over
c     ninold.
      if (errrat .gt. two) then
         errrat = two
      else
         if (errrat .lt. point5) then
            errrat = point5
         endif
      endif

      nint = int(ninold * errrat)

c     The code does not allow nint = ninold.
      if (nint .eq. ninold) then
         nint = nint + 1
      endif

c     Stop now if nint > nintmx and let bacolri pick up the error.
      if (nint .gt. nintmx) return

   20 continue

c-----------------------------------------------------------------------
c     Update work(ixold) to be the mesh before remeshing.
      do 30 i = 1, ninold + 1
         work(ixold-1+i) = x(i)
   30 continue

c-----------------------------------------------------------------------
c     Store work(i) to be the sum of the error at the first i
c     subintervals.
      work(ierror) = errint(1)
      do 40 i = ierror-1+2, ninold
         work(i) = errint(i) + work(i-1)
   40 continue

c     Let aerr to be the mean value of errint(i).
      aerr = work(ninold)/dble(nint)

c     Equidistribute the mesh points.
      berr = aerr
      j = 1

      do 60 i = 2, nint
   50    continue
         if (berr .gt. work(j)) then
            j = j + 1
            goto 50
         else
            if (j .eq. 1) then
               x(i) = work(ixold) + (work(ixold-1+2) - work(ixold))
     &                * berr/work(1)
            else
               x(i) = work(ixold-1+j) + (work(ixold-1+j+1) -
     &                work(ixold-1+j)) * (berr - work(j-1))/errint(j)
            endif
         endif
         berr = berr + aerr
   60 continue

      x(1) = work(ixold)
      x(nint+1) = work(ixold-1+ninold+1)

      return
      end
      subroutine solout(neq, t0, y, itol, rpar, ipar, irtrn)

c-----------------------------------------------------------------------
c Purpose:
c       This routine is called at each successful time step. It decides
c       whether a remeshing is necessary or not.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, Jun 5, 2003.
c
c------------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 neq
c                               neq=npde*ncpts is the number of bspline
c                               coefficients (or DAEs).
c
        double precision        t0
c                               t0 is the current time, t0 <= tout.
c
        double precision        y(neq)
c                               y is the vector of bspline
c                               coefficients at the current time.
c
        integer                 itol
c                               itol = 0, both rtol and atol are
c                                         scalars;
c                                    = 1, both rtol and atol are
c                                         vectors.
c
        double precision        rpar(*)
c                               rpar is the BACOLRI floating point work
c                               array.
c
        integer                 ipar(*)
c                               ipar is the BACOLRI integer work array.
c
c       output:
        integer                 irtrn
c                               irtrn is a status flag for remesh.
c                               irtrn = 0, indicates no need remeshing.
c                               irtrn = 1, indicates need remeshing.
c
c-----------------------------------------------------------------------
c Direct pointers into the IPAR integer work array:
        integer                 inpde
        parameter              (inpde  =  1)
c                               ipar(inpde) = npde
c
        integer                 ikcol
        parameter              (ikcol  =  2)
c                               ipar(ikcol) = kcol.
c
        integer                 inint
        parameter              (inint  =  3)
c                               ipar(inint) = nint.
c
        integer                 iest
        parameter              (iest   =  8)
c                               ipar(iest) = mflag(5), for LOI or SCI.
c
        integer                 istalr
        parameter              (istalr = 12)
c                               ipar(istalr) is the number of accepted
c                               steps after the last successful
c                               remeshing.
c
        integer                 irshfg
        parameter              (irshfg = 14)
c                               ipar(irshfg) is a flag for redefining
c                               all the pointers.
c                               ipar(irshfg) = 0, the initial step or
c                                                 any step not needing
c                                                 remesh;
c                                            = 1, a step needing remesh.
c
        integer                 icount
        parameter              (icount = 15)
c                               ipar(icount) is the number of remeshing
c                               times at the current step.
c
        integer                 istart
        parameter              (istart = 16)
c                               ipar(istart) is a flag to begin the
c                               code.
c                               ipar(istart) = 0, the initial step;
c                                            = 1, not the initial step.
c
        integer                 iatol
        parameter              (iatol  = 32)
c                               rpar(ipar(iatol)) = atol.
c
        integer                 irtol
        parameter              (irtol  = 33)
c                               rpar(ipar(irtol)) = rtol.
c
c-----------------------------------------------------------------------
c Direct pointers into the RPAR floating point work array:
        integer                 ierrat
        parameter              (ierrat =  3)
c                               rpar(ierrat) = the value of the largest
c                               component of rpar(ipar(iercom)).
c
        integer                 it0
        parameter              (it0    =  4)
c                               rpar(it0)    = t0 at the last accepted
c                               time step.
c
c-----------------------------------------------------------------------
c Indirect pointers into the RPAR floating point work array:
        integer                 ih
        parameter              (ih     = 21)
c                               rpar(ipar(ih)) stores the mesh step
c                               size sequence.
c
        integer                 ixbs
        parameter              (ixbs  = 23)
c                               rpar(ipar(ixbs)) stores the breakpoint
c                               sequence when using radau_{kcol}.
c
        integer                 iexcol
        parameter              (iexcol = 34)
c                               rpar(ipar(iexcol)) stores the
c                               collocation points which are used for
c                               error estimate.
c
        integer                 iewts
        parameter              (iewts  = 35)
c                               rpar(ipar(iewts)) stores the gaussian
c                               weights which are used for error
c                               estimate.
c
        integer                 iebas
        parameter              (iebas = 36)
c                               rpar(ipar(iebas)) stores the values
c                               of the nonzero basis functions at
c                               rpar(ipar(iexcol)) when using
c                               dassl_{kcol}.
c
        integer                 iecoef
        parameter              (iecoef = 37)
c                               rpar(ipar(iebas2)) stores the values
c                               of the nonzero basis functions at
c                               rpar(ipar(iexcol)) when using
c                               dassl_{kcol+1}.
c
        integer                 iercom
        parameter              (iercom = 38)
c                               rpar(ipar(iercom)) stores the error
c                               estimate for each component.
c
        integer                 ierint
        parameter              (ierint = 39)
c                               rpar(ipar(ierint)) stores the error
c                               estimate at each subinterval.
c
        integer                 iework
        parameter              (iework = 40)
c                               rpar(ipar(iework)) stores the floating
c                               point work array for errest.
c
        integer                 ix
        parameter              (ix     = 50)
c                               rpar(ipar(ix)) is the copy of x in rpar
c                               for use by errest.
c
        integer                 iypre
        parameter              (iypre  = 53)
c                               rpar(ipar(iypre)) stores the values of
c                               rpar(ipar(iy2)) at the previous 6 steps.
c                               It is required for a hot restart after
c                               remeshing.
c
c-----------------------------------------------------------------------
c Local variables:
        integer                 npde, kcol, nint
c
        integer                 necpts
c                               necpts=(kcol+3)*nint is the total number
c                               of collocation points used for
c                               error estimate.
c
        integer                 lenerr
c                               lenerr is the size of the floating point
c                               work array used by ERREST.
c                               lenerr>=2*npde*necpts+npde*nint.
c
        integer                 lencof
c                               Length of the work storage array for
c                               Hermite-Birkhoff coefficients and
c                               values of B-spline basis functions used
c                               inside errest's interpolation
c                               subroutines. They only need to be
c                               recalculated when the mesh changes.
c                               The SCI estimate requires more storage
c                               than the LOI estimate here.
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c                               dcopy
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               errest
c
c-----------------------------------------------------------------------

      kcol = ipar(ikcol)
      nint = ipar(inint)
      npde = ipar(inpde)

c     Calculate the number of quadrature points used for error
c     estimate and the extra storage requirements of errest.
c     The SCI scheme needs slightly more storage here.
      if (ipar(iest) .eq. 0) then
         necpts = (kcol + 2) * nint
         lenerr = (2 * necpts + nint) * npde
     &          + ((kcol - 3) * nint + (nint + 1) * 2) * npde
         lencof = (kcol + 1) * (kcol + 2)
     &          + (kcol + nconti) * (kcol - 3) * nint
     &          + (kcol + nconti) * 2 * (nint + 1)
      elseif (ipar(iest) .eq. 1) then
         necpts = (kcol + 3) * nint
         lenerr = (2 * necpts + nint) * npde
     &          + ((kcol - 2) * nint + (nint + 1) * 2) * npde
         lencof = (kcol + 4) * (kcol + 3) * nint
     &          + (kcol + nconti) * (kcol - 2) * nint
     &          + (kcol + nconti) * 2 * (nint + 1)
      endif

      call errestri(kcol, nint, npde, neq, necpts, ipar(icount),
     &            rpar(ipar(iexcol)), rpar(ipar(iewts)),
     &            rpar(ipar(ixbs)), y, ipar(istart), itol,
     &            rpar(ipar(iatol)), rpar(ipar(irtol)), lenerr,
     &            rpar(ipar(iework)), rpar(ipar(iebas)),
     &            rpar(ipar(iecoef)), lencof, rpar(ierrat),
     &            rpar(ipar(ierint)), rpar(ipar(iercom)), irtrn,
     &            rpar(ipar(ix)), rpar(ipar(ih)), ipar(iest))

      if (irtrn .ne. 0) goto 100

c     The current step is accepted.
      if (ipar(icount) .ne. 0) then
         ipar(istalr) = 1
         ipar(irshfg) = 0
      else
         ipar(istalr) = ipar(istalr) + 1
      endif

c     Update the backup information.
      call dcopyri(neq, y, 1, rpar(ipar(iypre)), 1)
      ipar(icount) = 0
      ipar(istart) = 1
      rpar(it0)    = t0

  100 continue

      return
      end

      subroutine valuesri(kcol, xsol, nint, x, npde, npts, nderiv,
     &                  usol, y, work)

c-----------------------------------------------------------------------
c Purpose:
c     This routine computes the solution u and the first nderv
c     derivatives of u at the npts points xsol. It then returns the
c     values in the array usol.

c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c-----------------------------------------------------------------------

c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 kcol
c       kcol is the number of collocation points to be used in
c       each subinterval.
c
        integer                 npts
c       npts is the number of points in the x vector.
c
        double precision        xsol(npts)
c       xsol is an arbitrary set of spatial points at which the solution
c       and the first nderv derivative values are to be calculated.
c
        integer                 nint
c       nint >= 1 is the number of subintervals defined by the spatial
c       mesh x.
c
        double precision        x(nint+1)
c       x is the spatial mesh which divides the interval [x_a,x_b] into
c       x_a = x(1) < x(2) < x(3) < ... < x(nint+1) = x_b.
c
        integer                 npde
c       npde is the number of components in the system of PDEs.
c       npde > 0.
c
        integer                 nderiv
c       nderiv is the number of derivatives of the solution which are
c       to be calculated.
c
        double precision        y(npde*(nint*kcol+nconti))
c       y is the vector of bspline coefficients at the final time step.
c
c       output:
        double precision        usol(npde, npts, nderiv+1)
c       usol is the solution and the spatial derivatives up to the
c       nderiv-th derivative at the given points and at the final time
c       step.
c
c       Work Storage:
        double precision        work((kcol+nconti)*(nderiv+1)
     *                                 +kcol*(nint+1)+2*nconti)
c       work is a floating point work storage array.
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        integer                 ileft
c                               breakpoint information.
c
        integer                 mflag
c                               mflag is required by subroutine
c                               interv.
c
        integer                 ilo
c                               ilo is required by subroutine
c                               interv.
c       Pointers into the floating point work array:
        integer                 ixbs
c                               work(ixbs) contains the breakpoint
c                               sequence.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 k
        integer                 m
        integer                 ii
        integer                 mj
        integer                 mm
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bsplvd
c                               interv
c
c-----------------------------------------------------------------------

c     set up the value for ileft, mflag and ncpts.
      ileft = 0
      mflag = -2
      ncpts = nint * kcol + nconti

c     set the pointer into the floating point work array
      ixbs  = (kcol+nconti)*(nderiv+1) + 1

c     Store the piecewise polynomial space breakpoint sequence in
c     work(ixbs).
c
      do 10 i = 1, kcol + nconti
         work(ixbs-1+i) = x(1)
         work(ixbs-1+i+ncpts) = x(nint+1)
   10 continue
      do 30 i = 2, nint
         ii = (i-2) * kcol + kcol + nconti
         do 20 k = 1, kcol
            work(ixbs-1+ii+k) = x(i)
   20 continue
   30 continue

      do 70 i = 1, npts
c
c     interv is called to compute ileft. bsplvd is called to compute
c     the values of the basis function at the required point.
         call intervri(work(ixbs), ncpts, xsol(i), ileft, mflag, ilo)
         call bsplvdri(work(ixbs),kcol+nconti,xsol(i),ileft,work,
     &               nderiv+1)
         ii = ileft - kcol - nconti
         do 60 j = 1, nderiv + 1
            do 50 k = 1, npde
               usol(k,i,j) = zero
               do 40 m = 1, kcol + nconti
                  mm = (m + ii - 1) * npde
                  mj = (j - 1) * (kcol + nconti) + m
                  usol(k,i,j) = usol(k,i,j) + y(mm+k) * work(mj)
   40          continue
   50       continue
   60    continue
   70 continue
      return
      end
      subroutine lowintri(kcol, nint, npde, istart, icount, xbs,
     &                  x, y, h, hd, g, bassc, basm, u, um, usol)
c-----------------------------------------------------------------------
c     This subroutine evaluates a Hermite-Birkhoff (H-B) interpolant of
c     order kcol (degree kcol-1) at kcol+2 Gaussian quadrature points on
c     the interval (0,1). (These values will be used elsewhere to obtain
c     an L2-norm error estimate.)
c
c     Last modified by Jack Pew, September 1, 2011.
c
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c Constants
      integer                 mxkcol
      parameter              (mxkcol = 10)
c     bacoli can have up to 10 collocation points per subinterval.
c
      integer                 nconti
      parameter              (nconti = 2)
c     nconti continuity conditions are imposed at the internal mesh
c     points of the collocation solution.
c
      double precision        zero
      parameter              (zero = 0.0d0)
c
c-----------------------------------------------------------------------
c Input
      integer                 kcol
c     Number of collocation points per subinterval for the collocation
c     solution. The collocation solution is a piecewise polynomial of
c     order kcol+1 (degree kcol.) The values obtained from this
c     interpolant are of order kcol. The interpolant is a piecewise
c     polynomial of degree kcol-1.
c
      integer                 nint
c     Number of subintervals on the current mesh.
c
      integer                 npde
c     Number of PDEs in the system being solved.
c
      integer                 istart
c     istart indicates whether or not the integration has just begun.
c         istart = 0, it is the initial step;
c                = 1, it is not the initial step.
c
      integer                 icount
c     icount is the number of remeshings that have been performed at
c     the current time step.
c     istart and icount are used to determine whether certain
c     coefficients need to be calculated, depending upon whether the
c     mesh has changed.
c
      double precision        xbs((nint+1)*kcol+nconti*2)
c     The breakpoint sequence. This is required by BSPLVD to compute the
c     values of the B-spline basis functions.
c
      double precision        x(nint+1)
c     The current mesh point sequence.
c
      double precision        y(npde*(nint*kcol+nconti))
c     B-spline coefficients computed by DASSL.
c
c-----------------------------------------------------------------------
c Inoutput
      double precision        h(2, kcol+2)
      double precision        hd(2, kcol+2)
      double precision        g(kcol-3, kcol+2)
c     Coefficients of the Hermite-Birkhoff interpolant.
c     h and hd coefficients multiply the values and derivatives of the
c     collocation solution at mesh points, g coefficients multiply the
c     values of the collocation solution at the internal points.
c     For the LOI scheme, no points external to the subinterval are
c     used, and these coefficient do not depend upon the mesh.
c     This means they really only need to be computed once,
c     but then the values would need to be copied from place to place
c     whenever nint changed and rpar grew or shrank. This copying would
c     potentially by of some confusion, as it would have to be done in
c     the main BACOLI loop. So while there is unnecessary work done here
c     (to recompute these coefficients,) we do it anyway since it takes
c     an insignificant fraction of the total CPU time.
c
      double precision        bassc(kcol+nconti, (kcol-3)*nint)
c     Values of all nonzero B-spline basis functions at the
c     interpolated points inside the subintervals.
c     These are stored between remeshings.
c
      double precision        basm(kcol+nconti, 2, nint+1)
c     Values of all nonzero B-spline basis functions and their first
c     derivatives at the mesh points. kcol of every kcol+nconti of these
c     are actually zero, but the trouble to reduce the storage used by
c     this vector seemed too great.
c     These are stored between remeshings.
c
c-----------------------------------------------------------------------
c Output
      double precision        usol(npde, nint*(kcol+2))
c     Array of evaluations of the Hermite-Birkhoff interpolant at the
c     quadrature points used for the computation of an L2-norm error
c     estimate elsewhere.
c
c-----------------------------------------------------------------------
c Work storage
      double precision        u(npde, (kcol-3)*nint)
c     Work array for storage of values of the collocation solution at
c     the internal points used in the Hermite-Birkhoff interpolant.
c
      double precision        um(2, npde, nint+1)
c     Work array for storage of values of the collocation solution, and
c     its first derivative, at the mesh points used in the Hermite-
c     Birkhoff interpolant.
c
c-----------------------------------------------------------------------
c Local variables
      double precision        scpt(mxkcol-2)
c     The kcol-2 internal interpolation points, over (0,1), used in the
c     construction of the Hermite-Birkhoff interpolant.
c     These points are not actually superconvergent for the
c     computed solution, as they would be in the SCI case. They are what
c     would be superconvergent points, if kcol were one lower.
c     This choice of points is important in order to have the leading
c     terms of the interpolation and collocation errors match.
c
      double precision        theta(mxkcol+2)
c     Quadrature points for a subinterval, which are the Gauss points
c     on (0,1).
c
      double precision        gwts(mxkcol+2)
c     Gaussian weights. Unused here.
c
      double precision        gwork((mxkcol+2)*(mxkcol+2))
c     Work storage for gauleg.
c
      integer                 ileft
c     Breakpoint information. ileft is used in conjunction with the
c     breakpoint sequence for the evaluation of B-spline basis functions
c     by BSPLVD.
c
      double precision        hx
c     Width of the ith subinterval, x(i+1) - x(i).
c
c-----------------------------------------------------------------------
c Loop indices
      integer                 i, j, k, m, ii, jj, mm
c
c-----------------------------------------------------------------------
c Subroutines Called
c                             bsplvd
c                             gauleg
c                             lowcoeff
c                             setpts
c
c-----------------------------------------------------------------------

c     Check whether the H-B coefficients and B-spline basis function
c     evaluations need to be computed. Previous values may be reused if
c     the mesh has not changed since the last step.
      if ((icount .ne. 0) .or. (istart .eq. 0)) then

c        Set thetas (weights will not be computed since last parameter
c        in the call to gauleg is set to 2.)
         call gaulegri(kcol+2, (kcol+2)*(kcol+2), theta, gwts, gwork, 2)

c        H-B coefficients are computed for the Gaussian quadrature
c        points. They are the same on every subinterval, so only one set
c        is needed.
c        Note that LOWCOEFF is called with kcol-1 instead of kcol. (This
c        is since for a given kcol and a collocation solution
c        corresponding to this given kcol, we actually want to construct
c        a LOI scheme that is one order lower.)
         do j = 1, kcol+2
            call lowcoeffri(theta(j), kcol-1, h(1,j), hd(1,j), g(1,j))
         end do

c        Call BSPLVD to obtain values of B-spline basis functions at the
c        internal interpolated points.
c        The subroutine INTERV does not need to be called to find ileft,
c        since we know which subinterval we are in already.
c
c        Note that setpts is called with kcol-1, so the points returned
c        are not in fact superconvergent for the computed collocation
c        solution. It is important that these points be used however so
c        that the leading terms of the interpolation and collocation
c        errors match.
         ileft = nconti
         do i = 1, nint
            call setptsri(scpt, kcol-1, i, x, nint)
            ileft = ileft + kcol
            do j = 1, kcol-3
               jj = (i-1) * (kcol-3) + j
               call bsplvdri(xbs, kcol+nconti, scpt(j), ileft,
     &                     bassc(1, jj), 1)
            end do
         end do

c        Call BSPLVD for the values and first derivatives of the
c        B-spline basis functions at the nint+1 mesh points.
c       (ileft does not follow the pattern at the rightmost mesh point,
c        so we need a separate call for the last mesh point.)
         ileft = nconti
         do i = 1, nint
            ileft = ileft + kcol
            call bsplvdri(xbs, kcol+nconti, x(i), ileft,
     &                  basm(1, 1, i), 2)
         end do
c        ileft = kcol+nconti+kcol*(nint-1)
         call bsplvdri(xbs, kcol+nconti, x(nint+1), ileft,
     &               basm(1, 1, nint+1), 2)

      endif

c-----------------------------------------------------------------------
c     Compute the values of the collocation solution at the internal
c     interpolated points and mesh points by multiplying the values of
c     the basis functions at those points by the B-spline coefficients
c     computed by DASSL.
c
c     The end result of these two steps (evaluating the basis functions
c     and later multiplying by the coefficients) is the same as calling
c     the subroutine VALUES, but this approach avoids a large number
c     of unnecessary calls to BSPLVD by saving some of the required
c     values.
      ii = 0
      do i = 1, nint
         do j = 1, kcol-3
            jj = (i-1) * (kcol-3) + j
            do k = 1, npde
               u(k, jj) = zero
               mm = ii
               do m = 1, kcol+nconti
c                 mm = (m + (i-1) * kcol - 1) * npde
                  u(k, jj) = u(k, jj) + y(mm+k) * bassc(m, jj)
                  mm = mm + npde
               end do
            end do
         end do
         ii = ii + kcol * npde
      end do
c     All but two of the values of the basis functions at the mesh
c     points are zero, so um gets special treatment.
      mm = 0
      do i = 1, nint
         do k = 1, npde
            um(1, k, i) = y(mm+k) * basm(1, 1, i)
     &             + y(npde+mm+k) * basm(2, 1, i)
            um(2, k, i) = y(mm+k) * basm(1, 2, i)
     &             + y(npde+mm+k) * basm(2, 2, i)
         end do
         mm = mm + kcol * npde
      end do
c     It is the final two, rather than the first two, values in the
c     first dimension of basm that are nonzero at the right endpoint.
c     Thus, this case needs to be handled separately.
      do k = 1, npde
c        mm = nint*kcol*npde from the previous loop
         um(1, k, nint+1) = y(mm+k) * basm(kcol+nconti-1, 1, nint+1)
     &               + y(npde+mm+k) * basm(kcol+nconti,   1, nint+1)
         um(2, k, nint+1) = y(mm+k) * basm(kcol+nconti-1, 2, nint+1)
     &               + y(npde+mm+k) * basm(kcol+nconti,   2, nint+1)
      end do

c-----------------------------------------------------------------------
c     With the solution values at our interpolation points extracted
c     from the B-spline interpolant (i.e. the collocation solution,) and
c     the Hermite-Birkhoff coefficients precomputed, we can now evaluate
c     Hermite-Birkhoff interpolant at the quadrature points for an error
c     estimate. (Recall that the Hermite-Birkhoff interpolant uses
c     collocation solution values as its interpolation data.)
c
c     The values are returned in the usol array.
      do i = 1, nint
         hx = (x(i+1) - x(i))
         ii = (i-1) * (kcol-3)
         do j = 1, kcol+2
            jj = (i-1) * (kcol+2) + j
            do m = 1, npde

               usol(m, jj) =  h(1, j) * um(1, m, i)
     &                     +  h(2, j) * um(1, m, i+1)
     &                     + hd(1, j) * um(2, m, i)   * hx
     &                     + hd(2, j) * um(2, m, i+1) * hx

               do k = 1, kcol-3
                  usol(m, jj) = usol(m, jj) + g(k, j) * u(m, ii+k)
               end do

            end do
         end do
      end do

      return
      end
      subroutine lowcoeffri(x, kcol, h, hd, g)
c-----------------------------------------------------------------------
c     This subroutine computes the Hermite-Birkhoff coefficients for a
c     given point at which the interpolant is to be evaluated.
c     They do not depend on the current subinterval.
c     This code uses a barycentric technique for phi, and exploits the
c     symmetry of the interpolation points.
c     Horner's rule is not used for phi, since that formulation lead to
c     some loss of precision.
c
c     Last modified by Jack Pew, September 1, 2011.
c
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c Input
      double precision        x
c     A point at which to evaluate the interpolant.
c     Typically a Gaussian quadrature point.
c
      integer                 kcol
c     The order of the polynomial is kcol+1.
c     Note that a value of kcol-1 is passed to this subroutine, however,
c     and therefore the interpolation points used here are not in fact
c     superconvergent (aside from the mesh points.) They would be
c     superconvergent for kcol-1, but are in effect arbitrary points
c     for kcol.
c
c-----------------------------------------------------------------------
c Output
      double precision        h(2)
      double precision        hd(2)
      double precision        g(kcol-2)
c     Hermite-Birkhoff coefficients for the given x.
c
c-----------------------------------------------------------------------
c Local variables
      double precision        phi
c     phi(x) = product(x-w_i), i=1,kcol-2
c     where w_i are the interpolation points.
c     A barycentric technique is applied to replace
c     phi(x)_j with (phi(x))/(x-w_j).
c
      double precision        eta1, eta2, eta
c     These factors are common between the coefficients.
c
c-----------------------------------------------------------------------
c Loop indices
      integer                 i
c
c-----------------------------------------------------------------------

c     precompute etas, which correspond to eta_1^2, eta_2^2 and eta^2
         eta1 = (x-1.d0)*(x-1.d0)
         eta2 = x*x
         eta = eta1*eta2

c     select which polynomial to evaluate for the H-B coefficients
      if (kcol .eq. 3) then
         phi = x-0.5d0
         h(1) = -2.d0*(1.d0+4.d0*x)*eta1*phi
         h(2) = 2.d0*(5.d0-4.d0*x)*eta2*phi
         hd(1) = -2.d0*x*eta1*phi
         hd(2) = 2.d0*(x-1.d0)*eta2*phi
         g(1) = 16.d0*eta

      elseif (kcol .eq. 4) then
         phi = (x-0.3110177634953864d0)*(x-0.6889822365046136d0)
         h(1) = (0.4666666666666667d1+0.3111111111111111d2*x)*eta1*phi
         h(2) = (0.3577777777777778d2-0.3111111111111111d2*x)*eta2*phi
         hd(1) = 0.4666666666666667d1*x*eta1*phi
         hd(2) = (0.4666666666666667d1*x-0.4666666666666667d1)*eta2*phi
         g(1) = -0.5761858410762886d2*eta
         g(2) = -g(1)*(x-0.3110177634953864d0)
         g(1) = g(1)*(x-0.6889822365046136d0)

      elseif (kcol .eq. 5) then
         phi = (x-0.2113248654051871d0)*(x-0.5d0)
     #         *(x-0.7886751345948129d0)
         h(1) = (-0.12d2-0.12d3*x)*eta1*phi
         h(2) = (0.132d3-0.12d3*x)*eta2*phi
         hd(1) = -0.12d2*x*eta1*phi
         hd(2) = (0.12d2*x-0.12d2)*eta2*phi
         g(1) = 0.216d3*eta*phi
         g(2) = -192.d0*eta*phi/(x-0.5d0)
         g(3) = g(1)/(x-0.7886751345948129d0)
         g(1) = g(1)/(x-0.2113248654051871d0)

      elseif (kcol .eq. 6) then
         phi = (x-0.1526267046965671d0)*(x-0.3747185964571342d0)
     #         *(x-0.6252814035428658d0)*(x-0.8473732953034329d0)
         h(1) = (0.33d2+0.462d3*x)*eta1*phi
         h(2) = (0.495d3-0.462d3*x)*eta2*phi
         hd(1) = 0.33d2*x*eta1*phi
         hd(2) = (-0.33d2+0.33d2*x)*eta2*phi
         g(1) = -0.8197591840300222d3*eta*phi
         g(2) = 0.6925405260332655d3*eta*phi
         g(3) = -g(2)/(x-0.6252814035428658d0)
         g(4) = -g(1)/(x-0.8473732953034329d0)
         g(1) = g(1)/(x-0.1526267046965671d0)
         g(2) = g(2)/(x-0.3747185964571342d0)

      elseif (kcol .eq. 7) then
         phi = (x-0.1152723378341063d0)*(x-0.2895425974880943d0)
     #         *(x-0.5d0)*(x-0.7104574025119057d0)
     #         *(x-0.8847276621658937d0)
         h(1) = (-0.9533333333333333d2-0.1779555555555556d4*x)*eta1*phi
         h(2) = (0.1874888888888889d4-0.1779555555555556d4*x)*eta2*phi
         hd(1) = -0.9533333333333333d2*x*eta1*phi
         hd(2) = (-0.9533333333333333d2+0.9533333333333333d2*x)*eta2*phi
         g(1) = 0.3131255164938139d4*eta*phi
         g(2) = -0.257196627604925d4*eta*phi
         g(3) = 0.2440533333333333d4*eta*phi/(x-0.5d0)
         g(4) = g(2)/(x-0.7104574025119057d0)
         g(5) = g(1)/(x-0.8847276621658937d0)
         g(1) = g(1)/(x-0.1152723378341063d0)
         g(2) = g(2)/(x-0.2895425974880943d0)

      elseif (kcol .eq. 8) then
         phi = (x-0.9007700226825652d-1)*(x-0.2296976813063206d0)
     #         *(x-0.405661288754607d0)*(x-0.594338711245393d0)
     #         *(x-0.7703023186936794d0)*(x-0.9099229977317435d0)
         h(1) = (0.286d3+0.6864d4*x)*eta1*phi
         h(2) = (0.715d4-0.6864d4*x)*eta2*phi
         hd(1) = 0.286d3*x*eta1*phi
         hd(2) = (-0.286d3+0.286d3*x)*eta2*phi
         g(1) = -0.1201314415336355d5*eta*phi
         g(2) = 0.9696023778080566d4*eta*phi
         g(3) = -0.8929458911121293d4*eta*phi
         g(4) = -g(3)/(x-0.594338711245393d0)
         g(5) = -g(2)/(x-0.7703023186936794d0)
         g(6) = -g(1)/(x-0.9099229977317435d0)
         g(1) = g(1)/(x-0.9007700226825652d-1)
         g(2) = g(2)/(x-0.2296976813063206d0)
         g(3) = g(3)/(x-0.405661288754607d0)

      elseif (kcol .eq. 9) then
         phi = (x-0.7229898685756272d-1)*(x-0.1863109301186906d0)
     #         *(x-0.3341852231986051d0)*(x-0.5d0)*(x-0.66581477680139
     #         49d0)*(x-0.8136890698813094d0)*(x-0.9277010131424373d0)
         h(1) = (-0.884d3-0.2652d5*x)*eta1*phi
         h(2) = (0.27404d5-0.2652d5*x)*eta2*phi
         hd(1) = -0.884d3*x*eta1*phi
         hd(2) = (-0.884d3+0.884d3*x)*eta2*phi
         g(1) = 0.4624522420172525d5*eta*phi
         g(2) = -0.3688889968687672d5*eta*phi
         g(3) = 0.3332824691372289d5*eta*phi
         g(4) = -0.3232914285714286d5*eta*phi/(x-0.5d0)
         g(5) = g(3)/(x-0.6658147768013949d0)
         g(6) = g(2)/(x-0.8136890698813094d0)
         g(7) = g(1)/(x-0.9277010131424373d0)
         g(1) = g(1)/(x-0.7229898685756272d-1)
         g(2) = g(2)/(x-0.1863109301186906d0)
         g(3) = g(3)/(x-0.3341852231986051d0)

      else
         h(1) = 0.d0
         h(2) = 0.d0
         hd(1) = 0.d0
         hd(2) = 0.d0
         do i =1,7
           g(i) = 0.d0
         end do
      endif

      return
      end
      subroutine setptsri(scpt, kcol, i, x, nint)
c-----------------------------------------------------------------------
c     This subroutine fills the scpt array with the superconvergent
c     points for the given mesh subinterval and kcol.
c
c     Last modified by Jack Pew, September 1, 2011.
c
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c Constants
      integer                 mxkcol
      parameter              (mxkcol = 10)
c     bacoli can have up to 10 collocation points per subinterval.
c
c-----------------------------------------------------------------------
c Input
      integer                 kcol
c     The value of kcol for which the superconvergent points are located.
c
      integer                 nint
c     Number of subintervals in the current mesh.
c
      double precision        x(nint+1)
c     The current mesh.
c
      integer                 i
c     Subinterval to which the superconvergent points for the given kcol
c     over [0,1] are mapped to.
c
c-----------------------------------------------------------------------
c Output
      double precision        scpt(mxkcol-2)
c     Superconvergent points for the given kcol mapped to the given
c     mesh subinterval.
c
c-----------------------------------------------------------------------
c Local variables
      double precision        h
c     The width of the given mesh subinterval.
c
c-----------------------------------------------------------------------

      h = (x(i+1) - x(i))

      if     (kcol .eq. 3) then
         scpt(1) = x(i) + 0.5000000000000000d0 * h

      elseif (kcol .eq. 4) then
         scpt(1) = x(i) + 0.3110177634953864d0 * h
         scpt(2) = x(i) + 0.6889822365046136d0 * h

      elseif (kcol .eq. 5) then
         scpt(1) = x(i) + 0.2113248654051871d0 * h
         scpt(2) = x(i) + 0.5000000000000000d0 * h
         scpt(3) = x(i) + 0.7886751345948129d0 * h

      elseif (kcol .eq. 6) then
         scpt(1) = x(i) + 0.1526267046965671d0 * h
         scpt(2) = x(i) + 0.3747185964571342d0 * h
         scpt(3) = x(i) + 0.6252814035428658d0 * h
         scpt(4) = x(i) + 0.8473732953034329d0 * h

      elseif (kcol .eq. 7) then
         scpt(1) = x(i) + 0.1152723378341063d0 * h
         scpt(2) = x(i) + 0.2895425974880943d0 * h
         scpt(3) = x(i) + 0.5000000000000000d0 * h
         scpt(4) = x(i) + 0.7104574025119057d0 * h
         scpt(5) = x(i) + 0.8847276621658937d0 * h

      elseif (kcol .eq. 8) then
         scpt(1) = x(i) + 0.9007700226825652d-1* h
         scpt(2) = x(i) + 0.2296976813063206d0 * h
         scpt(3) = x(i) + 0.4056612887546070d0 * h
         scpt(4) = x(i) + 0.5943387112453930d0 * h
         scpt(5) = x(i) + 0.7703023186936794d0 * h
         scpt(6) = x(i) + 0.9099229977317435d0 * h

      elseif (kcol .eq. 9) then
         scpt(1) = x(i) + 0.7229898685756272d-1* h
         scpt(2) = x(i) + 0.1863109301186906d0 * h
         scpt(3) = x(i) + 0.3341852231986051d0 * h
         scpt(4) = x(i) + 0.5000000000000000d0 * h
         scpt(5) = x(i) + 0.6658147768013949d0 * h
         scpt(6) = x(i) + 0.8136890698813094d0 * h
         scpt(7) = x(i) + 0.9277010131424373d0 * h

      elseif (kcol .eq. 10) then
         scpt(1) = x(i) + 0.5929571219129399d-1* h
         scpt(2) = x(i) + 0.1539696908715823d0 * h
         scpt(3) = x(i) + 0.2792835119457421d0 * h
         scpt(4) = x(i) + 0.4241841678533667d0 * h
         scpt(5) = x(i) + 0.5758158321466333d0 * h
         scpt(6) = x(i) + 0.7207164880542579d0 * h
         scpt(7) = x(i) + 0.8460303091284177d0 * h
         scpt(8) = x(i) + 0.9407042878087060d0 * h

      endif

      return
      end
      subroutine scintri(kcol, nint, npde, istart, icount, xbs,
     &                 x, y, h, hd, g, bassc, basm, u, um, usol)
c-----------------------------------------------------------------------
c     This subroutine evaluates a Hermite-Birkhoff (H-B) interpolant of
c     order kcol+2 (degree kcol+1) at kcol+3 Gaussian quadrature points
c     on the interval (0,1). (These values will be used elsewhere to
c     obtain an L2-norm error estimate.)
c
c     Last modified by Jack Pew, September 1, 2011.
c
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c Constants
      integer                 mxkcol
      parameter              (mxkcol = 10)
c     bacoli can have up to 10 collocation points per subinterval.
c
      integer                 nconti
      parameter              (nconti = 2)
c     nconti continuity conditions are imposed at the internal mesh
c     points of the collocation solution.
c
      double precision        zero
      parameter              (zero = 0.0d0)
c
c-----------------------------------------------------------------------
c Input
      integer                 kcol
c     Number of collocation points per subinterval for the collocation
c     solution. The collocation solution is a piecewise polynomial of
c     order kcol+1 (degree kcol.) The values obtained from this
c     interpolant are of order kcol+2. The interpolant is a piecewise
c     polynomial of degree kcol+1.
c
      integer                 nint
c     Number of subintervals in the current mesh.
c
      integer                 npde
c     Number of PDEs in the system being solved.
c
      integer                 istart
c     istart indicates whether or not the integration has just begun.
c         istart = 0, it is the initial step;
c                = 1, it is not the initial step.
c
      integer                 icount
c     icount is the number of remeshings that have been performed at
c     the current time step.
c     istart and icount are used to determine whether certain
c     coefficients need to be calculated, depending upon whether the
c     mesh has changed.
c
      double precision        xbs((nint+1)*kcol+nconti*2)
c     The breakpoint sequence. This is required by BSPLVD to compute the
c     values of the B-spline basis functions.
c
      double precision        x(nint+1)
c     The current mesh point sequence.
c
      double precision        y(npde*(nint*kcol+nconti))
c     B-spline coefficients computed by DASSL.
c
c-----------------------------------------------------------------------
c Inoutput
      double precision        h(2, kcol+3, nint)
      double precision        hd(2, kcol+3, nint)
      double precision        g(kcol, kcol+3, nint)
c     Coefficients of the Hermite-Birkhoff interpolant.
c     h and hd coefficients multiply the values and derivatives of the
c     collocation solution at mesh points, g coefficients multiply the
c     values of the collocation solution at the internal points.
c     These are stored between remeshings and depend upon mesh ratios
c     in the SCI case, since points external to the subinterval are
c     used. Thus, nint sets of these coefficients must be kept.
c
      double precision        bassc(kcol+nconti, (kcol-2)*nint)
c     Values of all nonzero B-spline basis functions at the
c     superconvergent points inside the subintervals.
c     These are stored between remeshings.
c
      double precision        basm(kcol+nconti, 2, nint+1)
c     Values of all nonzero B-spline basis functions and their first
c     derivatives at the mesh points. kcol of every kcol+nconti of these
c     are actually zero, but the trouble to reduce the storage used by
c     this vector seemed too great. These are stored between remeshings.
c
c-----------------------------------------------------------------------
c Output
      double precision        usol(npde, nint*(kcol+3))
c     Array of evaluations of the Hermite-Birkhoff interpolant at the
c     quadrature points used for the computation of an L2-norm error
c     estimate in elsewhere.
c
c-----------------------------------------------------------------------
c Work storage
      double precision        u(npde, (kcol-2)*nint)
c     Work array for storage of values of the collocation solution at
c     the internal points used in the Hermite-Birkhoff interpolant.
c
      double precision        um(2, npde, nint+1)
c     Work array for storage of values of the collocation solution, and
c     its first derivative, at the mesh points used in the Hermite-
c     Birkhoff interpolant.
c
c-----------------------------------------------------------------------
c Local variables
      double precision        scpt(mxkcol-2)
c     The kcol interpolated superconvergent points internal to each
c     subinterval. Some are used in more than one subinterval of the
c     Hermite-Birkhoff interpolant.
c
      double precision        theta(mxkcol+3)
c     Quadrature points for a subinterval, which are the Gauss points
c     on (0,1).
c
      double precision        gwts(mxkcol+3)
c     Gaussian weights. Unused here.
c
      double precision        gwork((mxkcol+3)*(mxkcol+3))
c     Work storage for gauleg.
c
      integer                 ileft
c     Breakpoint information. ileft is used in conjunction with the
c     breakpoint sequence for the evaluation of B-spline basis functions
c     by BSPLVD.
c
c-----------------------------------------------------------------------
c Loop indices
      integer                 i, j, k, m, jj, mm, ii
c
c-----------------------------------------------------------------------
c Subroutines Called
c                             bsplvd
c                             gauleg
c                             scicoeff
c                             setpts
c
c-----------------------------------------------------------------------

c     Check whether the H-B coefficients and B-spline basis function
c     evaluations need to be computed.
c     Previous values may be reused if the mesh has not changed
c     since the last step.
      if ((icount .ne. 0) .or. (istart .eq. 0)) then

c        Set thetas (weights will not be computed since last parameter
c        in the call to gauleg is set to 2.)
         call gaulegri(kcol+3, (kcol+3)*(kcol+3), theta, gwts, gwork, 2)

c        H-B coefficients are computed for the Gaussian quadrature
c        points at which the H-B interpolant is evaluated.
c        Since two points external to the subinterval are used in the
c        interpolant piece for each subinterval, all the coefficients
c        end up depending upon mesh subinteval ratios.
         call scicoeffri(theta, kcol, nint, x, h, hd, g)

c        Call BSPLVD to obtain values of the B-spline basis functions at
c        the internal superconvergent points used in the interpolant.
c        The subroutine INTERV does not need to be called to find ileft,
c        since we know which subinterval we are in already.
         ileft = nconti
         do i = 1, nint
            call setptsri(scpt, kcol, i, x, nint)
            ileft = ileft + kcol
            do j = 1, kcol-2
               jj = (i-1) * (kcol-2) + j
               call bsplvdri(xbs, kcol+nconti, scpt(j), ileft,
     &                     bassc(1, jj), 1)
            end do
         end do

c        Call BSPLVD for the values and first derivatives of the
c        B-spline basis functions at the nint+1 mesh points.
c       (ileft does not follow the pattern at the rightmost mesh point,
c        so we need a separate call for the last mesh point.)
         ileft = nconti
         do i = 1, nint
            ileft = ileft + kcol
            call bsplvdri(xbs, kcol+nconti, x(i), ileft,
     &                  basm(1, 1, i), 2)
         end do
c        ileft = kcol+nconti+kcol*(nint-1)
         call bsplvdri(xbs, kcol+nconti, x(nint+1), ileft,
     &               basm(1, 1, nint+1), 2)

      endif

c-----------------------------------------------------------------------
c     Compute the values of the collocation solution at the
c     superconvergent and mesh points by multiplying the values of the
c     basis functions at those points by the B-spline coefficients
c     computed by DASSL.
c
c     The end result of these two steps (evaluating the basis functions
c     and later multiplying by the coefficients) is the same as calling
c     the subroutine VALUES, but this approach avoids a large number
c     of unnecessary calls to BSPLVD by saving some of the required
c     values between remeshings.
      ii = 0
      do i = 1, nint
         do j = 1, kcol-2
            jj = (i-1) * (kcol-2) + j
            do k = 1, npde
               u(k, jj) = zero
               mm = ii
               do m = 1, kcol+nconti
c                 mm = (m + (i-1) * kcol - 1) * npde
                  u(k, jj) = u(k, jj) + y(mm+k) * bassc(m, jj)
                  mm = mm + npde
               end do
            end do
         end do
         ii = ii + kcol * npde
      end do
c     All but two of the values of the basis functions at the mesh
c     points are zero, so um gets special treatment.
      mm = 0
      do i = 1, nint
         do k = 1, npde
            um(1, k, i) = y(mm+k) * basm(1, 1, i)
     &             + y(npde+mm+k) * basm(2, 1, i)
            um(2, k, i) = y(mm+k) * basm(1, 2, i)
     &             + y(npde+mm+k) * basm(2, 2, i)
         end do
         mm = mm + kcol * npde
      end do
c     It is the final two, rather than the first two, values in the
c     first dimension of basm that are nonzero at the right endpoint.
c     This case needs to be handled separately.
      do k = 1, npde
c        mm = nint*kcol*npde from the previous loop
         um(1, k, nint+1) = y(mm+k) * basm(kcol+nconti-1, 1, nint+1)
     &               + y(npde+mm+k) * basm(kcol+nconti,   1, nint+1)
         um(2, k, nint+1) = y(mm+k) * basm(kcol+nconti-1, 2, nint+1)
     &               + y(npde+mm+k) * basm(kcol+nconti,   2, nint+1)
      end do

c-----------------------------------------------------------------------
c     With the solution values at our interpolation points extracted
c     from the B-spline interpolant (i.e. the collocation solution,) and
c     the Hermite-Birkhoff coefficients precomputed, we can now evaluate
c     Hermite-Birkhoff interpolant at the quadrature points for an error
c     estimate. (Recall that the Hermite-Birkhoff interpolant uses
c     collocation solution values as its interpolation data.)
c
c     The values are returned in the usol array.
      if (kcol .ne. 3) then
         do i = 1, nint
            if (i .ne. 1 .and. i .ne. nint) then
               ii = (i-1) * (kcol-2) - 1
            elseif (i .eq. 1) then
               ii = (i-1) * (kcol-2)
            else
               ii = (i-1) * (kcol-2) - 2
            endif
            do j = 1, kcol+3
               jj = (i-1) * (kcol+3) + j
               do m = 1, npde
                  usol(m, jj) =  h(1, j, i) * um(1, m, i)
     &                        +  h(2, j, i) * um(1, m, i+1)
     &                        + hd(1, j, i) * um(2, m, i)
     &                        + hd(2, j, i) * um(2, m, i+1)
                  do k = 1, kcol
                     usol(m, jj) = usol(m, jj)
     &                           + g(k, j, i) * u(m, ii + k)
                  end do
               end do
            end do
         end do
c     This has a different form with kcol=3, which uses a third mesh
c     point when interpolating at boundary subintervals.
      else
         do i = 1, nint
            do j = 1, kcol+3
               jj = (i-1) * (kcol+3) + j
               do m = 1, npde
                  if (i .ne. 1 .and. i .ne. nint) then
                     usol(m, jj) = g(1, j, i) * u(m, i-1)
     &                           + g(2, j, i) * u(m, i)
     &                           + g(3, j, i) * u(m, i+1)
                  elseif (i .eq. 1) then
                     usol(m, jj) = g(1, j, 1) * u(m, 1)
     &                           + g(2, j, 1) * u(m, 2)
     &                           + g(3, j, 1) * um(1, m, 3)
                  else
                     usol(m, jj) = g(1, j, nint) * um(1, m, nint-1)
     &                           + g(2, j, nint) * u(m, nint-1)
     &                           + g(3, j, nint) * u(m, nint)
                  endif
                  usol(m, jj) = usol(m, jj)
     &                        +  h(1, j, i) * um(1, m, i)
     &                        +  h(2, j, i) * um(1, m, i+1)
     &                        + hd(1, j, i) * um(2, m, i)
     &                        + hd(2, j, i) * um(2, m, i+1)
               end do
            end do
         end do
      endif

      return
      end
      subroutine scicoeffri(x, kcol, nint, mesh, h, hd, g)
c-----------------------------------------------------------------------
c     This subroutine computes the Hermite-Birkhoff coefficients for the
c     kcol+3 quadrature points per subinterval at which the interpolant
c     is to be evaluated. They differ between subintervals since points
c     from adjacent subintervals are used to construct the interpolant.
c
c     Note that the hd coefficients are also scaled by the subinterval
c     width here.
c
c     Last modified by Jack Pew, September 1, 2011.
c
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c Constants
      integer                 mxkcol
      parameter              (mxkcol = 10)
c     Maximum value kcol may take.
c
c-----------------------------------------------------------------------
c Input
      double precision        x(kcol+3)
c     Points relative to each subinterval at which to evaluate the
c     interpolant.
c
      integer                 kcol
c     The value of kcol for which the superconvergent points are located.
c
      integer                 nint
c     Number of mesh subintervals.
c
      double precision        mesh(nint+1)
c     The mesh array.
c
c-----------------------------------------------------------------------
c Output
      double precision        h(2, kcol+3, nint)
      double precision        hd(2, kcol+3, nint)
      double precision        g(kcol, kcol+3, nint)
c     H-B coefficients evaluated for the current mesh at the points x
c     on each subinterval.
c
c-----------------------------------------------------------------------
c Local variables
      double precision        phi
c     phi(x) = product(x-w_i), i=1,kcol
c     where w_i are the superconvergent points.
c     A barycentric technique is applied to replace
c     phi(x)_j with (phi(x))/(x-w_j).
c
      double precision        a, aa, b, bb, hx
c     Superconvergent points from adjacent intervals,
c     and the subinterval size.
c     a and aa are the nearest two such points from the left,
c     b and bb are the nearest two such points from the right,
c     and hx is the current subinterval size.
c     While hx is not strictly a part of any coefficient, it is
c     multiplied here to save a little work.
c
      double precision        eta(3, mxkcol+3), gamma1, gamma2
c     Other precomputed factors.
c
c-----------------------------------------------------------------------
c Loop indices
      integer                 i, j, k
c
c-----------------------------------------------------------------------

c     Precompute etas, which correspond to eta_1^2, eta_2^2 and eta^2
      do j = 1, kcol+3
         eta(1, j) = (x(j)-1.d0) * (x(j)-1.d0)
         eta(2, j) = x(j) * x(j)
         eta(3, j) = eta(1, j) * eta(2, j)
      end do

c     Evaluate the coefficients for a given kcol and mesh at the
c     kcol+3 points in x.
c     This is greatly complicated over the LOI case by the use of points
c     external to a subinterval, bringing mesh ratios into the mix, and
c     necessitating a distinction between internal and boundary
c     subintervals.
      if (kcol .eq. 3) then
         do i = 1, nint
            hx = (mesh(i+1)-mesh(i))
            if (i .ne. 1 .and. i .ne. nint) then
               a = -0.5d0*(mesh(i)-mesh(i-1))/hx
               b = 1.d0+0.5d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -1.d0/a-4.d0-1.d0/b
               gamma2 = 1.d0/(1.d0-a)+4.d0+1.d0/(1.d0-b)
               h(1,kcol+3,i) = -2.d0/a/b
               h(2,kcol+3,i) = 2.d0/(1.d0-a)/(1.d0-b)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((a-0.5d0)*(a-b)*a*a
     #            *(a-1.d0)*(a-1.d0))
               g(2,kcol+3,i) = 16.d0/(0.5d0-a)/(0.5d0-b)
               g(3,kcol+3,i) = 1.d0/((b-a)*(b-0.5d0)*b*b*(b-1.d0)
     #            *(b-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-a)*(x(j)-0.5d0)*(x(j)-b)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
               end do

            elseif (i .eq. 1) then
               b = 1.d0+0.5d0*(mesh(i+2)-mesh(i+1))/hx
               bb = 1.d0+(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -4.d0-1.d0/b-1.d0/bb
               gamma2 = 4.d0+1.d0/(1.d0-b)+1.d0/(1.d0-bb)
               h(1,kcol+3,i) = -2.d0/b/bb
               h(2,kcol+3,i) = 2.d0/(1.d0-b)/(1.d0-bb)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 16.d0/(0.5d0-bb)/(0.5d0-b)
               g(2,kcol+3,i) = 1.d0/((b-bb)*(b-0.5d0)*b*b*(b-1.d0)
     #            *(b-1.d0))
               g(3,kcol+3,i) = 1.d0/((bb-b)*(bb-0.5d0)*bb*bb*(bb-1.d0)
     #            *(bb-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-0.5d0)*(x(j)-b)*(x(j)-bb)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-bb)
               end do

            else
               aa = -(mesh(i)-mesh(i-1))/hx
               a = -0.5d0*(mesh(i)-mesh(i-1))/hx
               gamma1 = -1.d0/aa-1.d0/a-4.d0
               gamma2 = 1.d0/(1.d0-aa)+1.d0/(1.d0-a)+4.d0
               h(1,kcol+3,i) = -2.d0/aa/a
               h(2,kcol+3,i) = 2.d0/(1.d0-aa)/(1.d0-a)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((aa-0.5d0)*(aa-a)*aa*aa
     #            *(aa-1.d0)*(aa-1.d0))
               g(2,kcol+3,i) = 1.d0/((a-0.5d0)*(a-aa)*a*a
     #            *(a-1.d0)*(a-1.d0))
               g(3,kcol+3,i) = 16.d0/(0.5d0-a)/(0.5d0-aa)

               do j = 1, kcol+3
                  phi = (x(j)-aa)*(x(j)-a)*(x(j)-0.5d0)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-aa)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
               end do
            endif
         end do

c-----------------------------------------------------------------------
      elseif (kcol .eq. 4) then
         do i = 1, nint
            hx = (mesh(i+1)-mesh(i))
            if (i .ne. 1 .and. i .ne. nint) then
               a = -0.3110177634953864d0*(mesh(i)-mesh(i-1))/hx
               b = 1.d0+0.3110177634953864d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -1.d0/a-0.6666666666666667d1-1.d0/b
               gamma2 = 1.d0/(1.d0-a)+0.6666666666666667d1+1.d0
     #            /(1.d0-b)
               h(1,kcol+3,i) = 0.4666666666666667d1/a/b
               h(2,kcol+3,i) = 0.4666666666666667d1/(1.d0-a)/(1.d0-b)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/(a-0.3110177634953864d0)/(a-0.6889
     #            822365046136d0)/(a-b)/a/a/(a-1.d0)/(a-1.d0)
               g(2,kcol+3,i) = -0.5761858410762886d2/(0.3110177634953
     #            864d0-a)/(0.3110177634953864d0-b)
               g(3,kcol+3,i) = 0.5761858410762886d2/(0.68898223650461
     #            36d0-a)/(0.6889822365046136d0-b)
               g(4,kcol+3,i) = 1.d0/(b-a)/(b-0.3110177634953864d0)
     #            /(b-0.6889822365046136d0)/b/b/(b-1.d0)/(b-1.d0)

               do j = 1, kcol+3
                  phi = (x(j)-a)*(x(j)-0.3110177634953864d0)*(x(j)
     #                  -0.6889822365046136d0)*(x(j)-b)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3110177634953864d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6889822365046136d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
               end do

            elseif (i .eq. 1) then
               b = 1.d0+0.3110177634953864d0*(mesh(i+2)-mesh(i+1))/hx
               bb = 1.d0+0.6889822365046136d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -0.6666666666666667d1-1.d0/b-1.d0/bb
               gamma2 = 0.6666666666666667d1+1.d0/(1.d0-b)+1.d0
     #            /(1.d0-bb)
               h(1,kcol+3,i) = 0.4666666666666667d1/b/bb
               h(2,kcol+3,i) = 0.4666666666666667d1/(1.d0-b)/(1.d0-bb)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = -0.5761858410762886d2/(0.3110177634953
     #            864d0-b)/(0.3110177634953864d0-bb)
               g(2,kcol+3,i) = 0.5761858410762886d2/(0.68898223650461
     #            36d0-b)/(0.6889822365046136d0-bb)
               g(3,kcol+3,i) = 1.d0/(b-bb)/(b-0.3110177634953864d0)
     #            /(b-0.6889822365046136d0)/b/b/(b-1.d0)/(b-1.d0)
               g(4,kcol+3,i) = 1.d0/(bb-b)/(bb-0.3110177634953864d0)
     #            /(bb-0.6889822365046136d0)/bb/bb/(bb-1.d0)/(bb-1.d0)

               do j = 1, kcol+3
                  phi = (x(j)-0.3110177634953864d0)*(x(j)
     #                  -0.6889822365046136d0)*(x(j)-b)*(x(j)-bb)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3110177634953864d0)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6889822365046136d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-bb)
               end do

            else
               aa = -0.6889822365046136d0*(mesh(i)-mesh(i-1))/hx
               a = -0.3110177634953864d0*(mesh(i)-mesh(i-1))/hx
               gamma1 = -1.d0/aa-1.d0/a-0.6666666666666667d1
               gamma2 = 1.d0/(1.d0-aa)+1.d0/(1.d0-a)
     #            +0.6666666666666667d1
               h(1,kcol+3,i) = 0.4666666666666667d1/aa/a
               h(2,kcol+3,i) = 0.4666666666666667d1/(1.d0-aa)/(1.d0-a)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/(aa-a)/(aa-0.3110177634953864d0)
     #            /(aa-0.6889822365046136d0)/aa/aa/(aa-1.d0)/(aa-1.d0)
               g(2,kcol+3,i) = 1.d0/(a-aa)/(a-0.3110177634953864d0)
     #            /(a-0.6889822365046136d0)/a/a/(a-1.d0)/(a-1.d0)
               g(3,kcol+3,i) = -0.5761858410762886d2/(0.3110177634953
     #            864d0-a)/(0.3110177634953864d0-aa)
               g(4,kcol+3,i) = 0.5761858410762886d2/(0.68898223650461
     #            36d0-a)/(0.6889822365046136d0-aa)

               do j = 1, kcol+3
                  phi = (x(j)-aa)*(x(j)-a)*(x(j)-0.3110177634953864d0)
     #                  *(x(j)-0.6889822365046136d0)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-aa)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3110177634953864d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6889822365046136d0)
               end do
            endif
         end do

c-----------------------------------------------------------------------
      elseif (kcol .eq. 5) then
         do i = 1, nint
            hx = (mesh(i+1)-mesh(i))
            if (i .ne. 1 .and. i .ne. nint) then
               a = -0.2113248654051871d0*(mesh(i)-mesh(i-1))/hx
               b = 1.d0+0.2113248654051871d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -1.d0/a-0.1d2-1.d0/b
               gamma2 = 1.d0/(1.d0-a)+0.1d2+1.d0/(1.d0-b)
               h(1,kcol+3,i) = -0.12d2/a/b
               h(2,kcol+3,i) = 0.12d2/(1.d0-a)/(1.d0-b)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((a-0.2113248654051871d0)*(a-0.5d0)
     #            *(a-0.7886751345948129d0)*(a-b)*a*a*(a-1.d0)*(a-1.d0))
               g(2,kcol+3,i) = 0.216d3/(0.2113248654051871d0-a)
     #            /(0.2113248654051871d0-b)
               g(3,kcol+3,i) = -192.d0/(0.5d0-a)/(0.5d0-b)
               g(4,kcol+3,i) = 0.216d3/(0.7886751345948129d0-a)
     #            /(0.7886751345948129d0-b)
               g(5,kcol+3,i) = 1.d0/((b-a)*(b-0.2113248654051871d0)*(b-
     #            0.5d0)*(b-0.7886751345948129d0)*b*b*(b-1.d0)*(b-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-a)*(x(j)-0.2113248654051871d0)*(x(j)
     #                  -0.5d0)*(x(j)-0.7886751345948129d0)*(x(j)-b)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2113248654051871d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7886751345948129d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
               end do

            elseif (i .eq. 1) then
               b = 1.d0+0.2113248654051871d0*(mesh(i+2)-mesh(i+1))/hx
               bb = 1.d0+0.5d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -0.1d2-1.d0/b-1.d0/bb
               gamma2 = 0.1d2+1.d0/(1.d0-b)+1.d0/(1.d0-bb)
               h(1,kcol+3,i) = -0.12d2/b/bb
               h(2,kcol+3,i) = 0.12d2/(1.d0-b)/(1.d0-bb)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 0.216d3/(0.2113248654051871d0-b)
     #            /(0.2113248654051871d0-bb)
               g(2,kcol+3,i) = -192.d0/(0.5d0-b)/(0.5d0-bb)
               g(3,kcol+3,i) = 0.216d3/(0.7886751345948129d0-b)
     #            /(0.7886751345948129d0-bb)
               g(4,kcol+3,i) = 1.d0/((b-0.2113248654051871d0)*(b-0.5d0)*
     #            (b-0.7886751345948129d0)*(b-bb)*b*b*(b-1.d0)*(b-1.d0))
               g(5,kcol+3,i) = 1.d0/((bb-0.2113248654051871d0)*(bb
     #            -0.5d0)*(bb-0.7886751345948129d0)*(bb-b)*bb*bb
     #            *(bb-1.d0)*(bb-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-0.2113248654051871d0)*(x(j)-0.5d0)*(x(j)
     #                  -0.7886751345948129d0)*(x(j)-b)*(x(j)-bb)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2113248654051871d0)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7886751345948129d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-bb)
               end do

            else
               aa = -0.5d0*(mesh(i)-mesh(i-1))/hx
               a = -0.2113248654051871d0*(mesh(i)-mesh(i-1))/hx
               gamma1 = -1.d0/aa-1.d0/a-0.1d2
               gamma2 = 1.d0/(1.d0-aa)+1.d0/(1.d0-a)+0.1d2
               h(1,kcol+3,i) = -0.12d2/aa/a
               h(2,kcol+3,i) = 0.12d2/(1.d0-aa)/(1.d0-a)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((aa-0.2113248654051871d0)*(aa
     #            -0.5d0)*(aa-0.7886751345948129d0)*(aa-a)*aa*aa
     #            *(aa-1.d0)*(aa-1.d0))
               g(2,kcol+3,i) = 1.d0/((a-0.2113248654051871d0)*(a-0.5d0)*
     #            (a-0.7886751345948129d0)*(a-aa)*a*a*(a-1.d0)*(a-1.d0))
               g(3,kcol+3,i) = 0.216d3/(0.2113248654051871d0-a)
     #            /(0.2113248654051871d0-aa)
               g(4,kcol+3,i) = -192.d0/(0.5d0-a)/(0.5d0-aa)
               g(5,kcol+3,i) = 0.216d3/(0.7886751345948129d0-a)
     #            /(0.7886751345948129d0-aa)

               do j = 1, kcol+3
                  phi = (x(j)-aa)*(x(j)-a)*(x(j)-0.2113248654051871d0)
     #                  *(x(j)-0.5d0)*(x(j)-0.7886751345948129d0)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-aa)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2113248654051871d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7886751345948129d0)
               end do
            endif
         end do

c-----------------------------------------------------------------------
      elseif (kcol .eq. 6) then
         do i = 1, nint
            hx = (mesh(i+1)-mesh(i))
            if (i .ne. 1 .and. i .ne. nint) then
               a = -0.1526267046965671d0*(mesh(i)-mesh(i-1))/hx
               b = 1.d0+0.1526267046965671d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -1.d0/a-0.14d2-1.d0/b
               gamma2 = 1.d0/(1.d0-a)+0.14d2+1.d0/(1.d0-b)
               h(1,kcol+3,i) = 0.33d2/a/b
               h(2,kcol+3,i) = 0.33d2/(1.d0-a)/(1.d0-b)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((a-0.1526267046965671d0)*(a-0.37471
     #            85964571342d0)*(a-0.6252814035428658d0)*(a-0.84737329
     #            53034329d0)*(a-b)*a*a*(a-1.d0)*(a-1.d0))
               g(2,kcol+3,i) = -0.8197591840300222d3/(0.152626704696567
     #            1d0-a)/(0.1526267046965671d0-b)
               g(3,kcol+3,i) = 0.6925405260332655d3/(0.37471859645713
     #            42d0-a)/(0.3747185964571342d0-b)
               g(4,kcol+3,i) = -0.6925405260332655d3/(0.625281403542865
     #            8d0-a)/(0.6252814035428658d0-b)
               g(5,kcol+3,i) = 0.8197591840300222d3/(0.847373295303432
     #            9d0-a)/(0.8473732953034329d0-b)
               g(6,kcol+3,i) = 1.d0/((b-a)*(b-0.1526267046965671d0)*(b
     #            -0.3747185964571342d0)*(b-0.6252814035428658d0)*(b
     #            -0.8473732953034329d0)*b*b*(b-1.d0)*(b-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-a)*(x(j)-0.1526267046965671d0)*(x(j)-0.3
     #                  747185964571342d0)*(x(j)-0.6252814035428658d0)
     #                  *(x(j)-0.8473732953034329d0)*(x(j)-b)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1526267046965671d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3747185964571342d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6252814035428658d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8473732953034329d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
               end do

            elseif (i .eq. 1) then
               b = 1.d0+0.1526267046965671d0*(mesh(i+2)-mesh(i+1))/hx
               bb = 1.d0+0.3747185964571342d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -0.14d2-1.d0/b-1.d0/bb
               gamma2 = 0.14d2+1.d0/(1.d0-b)+1.d0/(1.d0-bb)
               h(1,kcol+3,i) = 0.33d2/b/bb
               h(2,kcol+3,i) = 0.33d2/(1.d0-b)/(1.d0-bb)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = -0.8197591840300222d3/(0.152626704696567
     #            1d0-b)/(0.1526267046965671d0-bb)
               g(2,kcol+3,i) = 0.6925405260332655d3/(0.37471859645713
     #            42d0-b)/(0.3747185964571342d0-bb)
               g(3,kcol+3,i) = -0.6925405260332655d3/(0.625281403542865
     #            8d0-b)/(0.6252814035428658d0-bb)
               g(4,kcol+3,i) = 0.8197591840300222d3/(0.847373295303432
     #            9d0-b)/(0.8473732953034329d0-bb)
               g(5,kcol+3,i) = 1.d0/((b-bb)*(b-0.1526267046965671d0)*(b
     #            -0.3747185964571342d0)*(b-0.6252814035428658d0)*(b
     #            -0.8473732953034329d0)*b*b*(b-1.d0)*(b-1.d0))
               g(6,kcol+3,i) = 1.d0/((bb-b)*(bb-0.1526267046965671d0)*
     #            (bb-0.3747185964571342d0)*(bb-0.6252814035428658d0)*
     #            (bb-0.8473732953034329d0)*bb*bb*(bb-1.d0)*(bb-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-0.1526267046965671d0)*(x(j)-0.374718
     #                  5964571342d0)*(x(j)-0.6252814035428658d0)
     #                  *(x(j)-0.8473732953034329d0)*(x(j)-b)*(x(j)-bb)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1526267046965671d0)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3747185964571342d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6252814035428658d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8473732953034329d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-bb)
               end do

            else
               aa = -0.3747185964571342d0*(mesh(i)-mesh(i-1))/hx
               a = -0.1526267046965671d0*(mesh(i)-mesh(i-1))/hx
               gamma1 = -1.d0/aa-1.d0/a-0.14d2
               gamma2 = 1.d0/(1.d0-aa)+1.d0/(1.d0-a)+0.14d2
               h(1,kcol+3,i) = 0.33d2/aa/a
               h(2,kcol+3,i) = 0.33d2/(1.d0-aa)/(1.d0-a)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((aa-0.1526267046965671d0)*(aa-0.374
     #            7185964571342d0)*(aa-0.6252814035428658d0)*(aa-0.8473
     #            732953034329d0)*(aa-a)*aa*aa*(aa-1.d0)*(aa-1.d0))
               g(2,kcol+3,i) = 1.d0/((a-0.1526267046965671d0)*(a-0.37471
     #            85964571342d0)*(a-0.6252814035428658d0)*(a-0.84737329
     #            53034329d0)*(a-aa)*a*a*(a-1.d0)*(a-1.d0))
               g(3,kcol+3,i) = -0.8197591840300222d3/(0.152626704696567
     #            1d0-a)/(0.1526267046965671d0-aa)
               g(4,kcol+3,i) = 0.6925405260332655d3/(0.37471859645713
     #            42d0-a)/(0.3747185964571342d0-aa)
               g(5,kcol+3,i) = -0.6925405260332655d3/(0.625281403542865
     #            8d0-a)/(0.6252814035428658d0-aa)
               g(6,kcol+3,i) = 0.8197591840300222d3/(0.847373295303432
     #            9d0-a)/(0.8473732953034329d0-aa)

               do j = 1, kcol+3
                  phi = (x(j)-aa)*(x(j)-a)*(x(j)-0.1526267046965671d0)
     #                  *(x(j)-0.3747185964571342d0)*(x(j)-0.625281403
     #                  5428658d0)*(x(j)-0.8473732953034329d0)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-aa)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1526267046965671d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3747185964571342d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6252814035428658d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8473732953034329d0)
               end do
            endif
         end do

c-----------------------------------------------------------------------
      elseif (kcol .eq. 7) then
         do i = 1, nint
            hx = (mesh(i+1)-mesh(i))
            if (i .ne. 1 .and. i .ne. nint) then
               a = -0.1152723378341063d0*(mesh(i)-mesh(i-1))/hx
               b = 1.d0+0.1152723378341063d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -1.d0/a-0.1866666666666667d2-1.d0/b
               gamma2 = 1.d0/(1.d0-a)+0.1866666666666667d2+1.d0
     #            /(1.d0-b)
               h(1,kcol+3,i) = -0.9533333333333333d2/a/b
               h(2,kcol+3,i) = 0.9533333333333333d2/(1.d0-a)/(1.d0-b)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((a-0.1152723378341063d0)*(a-0.28954
     #            25974880943d0)*(a-0.5d0)*(a-0.7104574025119057d0)
     #            *(a-0.8847276621658937d0)*(a-b)*a*a*(a-1.d0)*(a-1.d0))
               g(2,kcol+3,i) = 0.3131255164938139d4/(0.115272337834106
     #            3d0-a)/(0.1152723378341063d0-b)
               g(3,kcol+3,i) = -0.257196627604925d4/(0.289542597488094
     #            3d0-a)/(0.2895425974880943d0-b)
               g(4,kcol+3,i) = 0.2440533333333333d4/(0.5d0-a)/(0.5d0-b)
               g(5,kcol+3,i) = -0.257196627604925d4/(0.710457402511905
     #            7d0-a)/(0.7104574025119057d0-b)
               g(6,kcol+3,i) = 0.3131255164938139d4/(0.884727662165893
     #            7d0-a)/(0.8847276621658937d0-b)
               g(7,kcol+3,i) = 1.d0/((b-a)*(b-0.1152723378341063d0)*(b-
     #            0.2895425974880943d0)*(b-0.5d0)*(b-0.710457402511905
     #            7d0)*(b-0.8847276621658937d0)*b*b*(b-1.d0)*(b-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-a)*(x(j)-0.1152723378341063d0)*(x(j)-0.28
     #                  95425974880943d0)*(x(j)-0.5d0)*(x(j)-0.71045740
     #                  25119057d0)*(x(j)-0.8847276621658937d0)*(x(j)-b)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1152723378341063d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2895425974880943d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7104574025119057d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8847276621658937d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
               end do

            elseif (i .eq. 1) then
               b = 1.d0+0.1152723378341063d0*(mesh(i+2)-mesh(i+1))/hx
               bb = 1.d0+0.2895425974880943d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -0.1866666666666667d2-1.d0/b-1.d0/bb
               gamma2 = 0.1866666666666667d2+1.d0/(1.d0-b)
     #            +1.d0/(1.d0-bb)
               h(1,kcol+3,i) = -0.9533333333333333d2/b/bb
               h(2,kcol+3,i) = 0.9533333333333333d2/(1.d0-b)/(1.d0-bb)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 0.3131255164938139d4/(0.115272337834106
     #            3d0-bb)/(0.1152723378341063d0-b)
               g(2,kcol+3,i) = -0.257196627604925d4/(0.289542597488094
     #            3d0-bb)/(0.2895425974880943d0-b)
               g(3,kcol+3,i) = 0.2440533333333333d4/(0.5d0-bb)/(0.5d0-b)
               g(4,kcol+3,i) = -0.257196627604925d4/(0.710457402511905
     #            7d0-bb)/(0.7104574025119057d0-b)
               g(5,kcol+3,i) = 0.3131255164938139d4/(0.884727662165893
     #            7d0-bb)/(0.8847276621658937d0-b)
               g(6,kcol+3,i) = 1.d0/((b-bb)*(b-0.1152723378341063d0)*(b-
     #            0.2895425974880943d0)*(b-0.5d0)*(b-0.710457402511905
     #            7d0)*(b-0.8847276621658937d0)*b*b*(b-1.d0)*(b-1.d0))
               g(7,kcol+3,i) = 1.d0/((bb-b)*(bb-0.1152723378341063d0)
     #            *(bb-0.2895425974880943d0)*(bb-0.5d0)*(bb-0.71045740
     #            25119057d0)*(bb-0.8847276621658937d0)*bb*bb*(bb-1.d0)
     #            *(bb-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-0.1152723378341063d0)*(x(j)-0.289542597488
     #                  0943d0)*(x(j)-0.5d0)*(x(j)-0.7104574025119057d0)
     #                  *(x(j)-0.8847276621658937d0)*(x(j)-b)*(x(j)-bb)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1152723378341063d0)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2895425974880943d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7104574025119057d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8847276621658937d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-bb)
               end do

            else
               aa = -0.2895425974880943d0*(mesh(i)-mesh(i-1))/hx
               a = -0.1152723378341063d0*(mesh(i)-mesh(i-1))/hx
               gamma1 = -1.d0/aa-1.d0/a-0.1866666666666667d2
               gamma2 = 1.d0/(1.d0-aa)+1.d0/(1.d0-a)
     #            +0.1866666666666667d2
               h(1,kcol+3,i) = -0.9533333333333333d2/aa/a
               h(2,kcol+3,i) = 0.9533333333333333d2/(1.d0-aa)/(1.d0-a)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((aa-0.1152723378341063d0)*(aa-0.28
     #            95425974880943d0)*(aa-0.5d0)*(aa-0.7104574025119057d0)
     #            *(aa-0.8847276621658937d0)*(aa-a)*aa*aa*(aa-1.d0)
     #            *(aa-1.d0))
               g(2,kcol+3,i) = 1.d0/((a-0.1152723378341063d0)*(a-0.28954
     #            25974880943d0)*(a-0.5d0)*(a-0.7104574025119057d0)*(a
     #            -0.8847276621658937d0)*(a-aa)*a*a*(a-1.d0)*(a-1.d0))
               g(3,kcol+3,i) = 0.3131255164938139d4/(0.115272337834106
     #            3d0-a)/(0.1152723378341063d0-aa)
               g(4,kcol+3,i) = -0.257196627604925d4/(0.289542597488094
     #            3d0-a)/(0.2895425974880943d0-aa)
               g(5,kcol+3,i) = 0.2440533333333333d4/(0.5d0-a)/(0.5d0-aa)
               g(6,kcol+3,i) = -0.257196627604925d4/(0.710457402511905
     #            7d0-a)/(0.7104574025119057d0-aa)
               g(7,kcol+3,i) = 0.3131255164938139d4/(0.884727662165893
     #            7d0-a)/(0.8847276621658937d0-aa)

               do j = 1, kcol+3
                  phi = (x(j)-aa)*(x(j)-a)*(x(j)-0.1152723378341063d0)*(
     #                  x(j)-0.2895425974880943d0)*(x(j)-0.5d0)*(x(j)-0.
     #                  7104574025119057d0)*(x(j)-0.8847276621658937d0)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-aa)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1152723378341063d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2895425974880943d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7104574025119057d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8847276621658937d0)
               end do
            endif
         end do

c-----------------------------------------------------------------------
      elseif (kcol .eq. 8) then
         do i = 1, nint
            hx = (mesh(i+1)-mesh(i))
            if (i .ne. 1 .and. i .ne. nint) then
               a = -0.9007700226825652d-1*(mesh(i)-mesh(i-1))/hx
               b = 1.d0+0.9007700226825652d-1*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -1.d0/a-0.24d2-1.d0/b
               gamma2 = 1.d0/(1.d0-a)+0.24d2+1.d0/(1.d0-b)
               h(1,kcol+3,i) = 0.286d3/a/b
               h(2,kcol+3,i) = 0.286d3/(1.d0-a)/(1.d0-b)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((a-0.9007700226825652d-1)*(a-0.229
     #            6976813063206d0)*(a-0.405661288754607d0)*(a-0.5943387
     #            11245393d0)*(a-0.7703023186936794d0)*(a-0.90992299773
     #            17435d0)*(a-b)*a*a*(a-1.d0)*(a-1.d0))
               g(2,kcol+3,i) = -0.1201314415336355d5/(0.900770022682565
     #            2d-1-a)/(0.9007700226825652d-1-b)
               g(3,kcol+3,i) = 0.9696023778080566d4/(0.229697681306320
     #            6d0-a)/(0.2296976813063206d0-b)
               g(4,kcol+3,i) = -0.8929458911121293d4/(0.40566128875460
     #            7d0-a)/(0.405661288754607d0-b)
               g(5,kcol+3,i) = 0.8929458911121293d4/(0.594338711245393d0
     #            -a)/(0.594338711245393d0-b)
               g(6,kcol+3,i) = -0.9696023778080565d4/(0.770302318693679
     #            4d0-a)/(0.7703023186936794d0-b)
               g(7,kcol+3,i) = 0.1201314415336355d5/(0.909922997731743
     #            5d0-a)/(0.9099229977317435d0-b)
               g(8,kcol+3,i) = 1.d0/((b-a)*(b-0.9007700226825652d-1)*(b
     #            -0.2296976813063206d0)*(b-0.405661288754607d0)*(b-0.5
     #            94338711245393d0)*(b-0.7703023186936794d0)*(b-0.90992
     #            29977317435d0)*b*b*(b-1.d0)*(b-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-a)*(x(j)-0.9007700226825652d-1)*(x(j)-0.2
     #                  296976813063206d0)*(x(j)-0.405661288754607d0)
     #                  *(x(j)-0.594338711245393d0)*(x(j)-0.7703023186
     #                  936794d0)*(x(j)-0.9099229977317435d0)*(x(j)-b)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9007700226825652d-1)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2296976813063206d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.405661288754607d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.594338711245393d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7703023186936794d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9099229977317435d0)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
               end do

            elseif (i .eq. 1) then
               b = 1.d0+0.9007700226825652d-1*(mesh(i+2)-mesh(i+1))/hx
               bb = 1.d0+0.2296976813063206d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -0.24d2-1.d0/b-1.d0/bb
               gamma2 = 0.24d2+1.d0/(1.d0-b)+1.d0/(1.d0-bb)
               h(1,kcol+3,i) = 0.286d3/b/bb
               h(2,kcol+3,i) = 0.286d3/(1.d0-b)/(1.d0-bb)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = -0.1201314415336355d5/(0.900770022682565
     #            2d-1-bb)/(0.9007700226825652d-1-b)
               g(2,kcol+3,i) = 0.9696023778080566d4/(0.229697681306320
     #            6d0-bb)/(0.2296976813063206d0-b)
               g(3,kcol+3,i) = -0.8929458911121293d4/(0.40566128875460
     #            7d0-bb)/(0.405661288754607d0-b)
               g(4,kcol+3,i) = 0.8929458911121293d4/(0.594338711245393d0
     #            -bb)/(0.594338711245393d0-b)
               g(5,kcol+3,i) = -0.9696023778080565d4/(0.770302318693679
     #            4d0-bb)/(0.7703023186936794d0-b)
               g(6,kcol+3,i) = 0.1201314415336355d5/(0.909922997731743
     #            5d0-bb)/(0.9099229977317435d0-b)
               g(7,kcol+3,i) = 1.d0/((b-bb)*(b-0.9007700226825652d-1)*(b
     #            -0.2296976813063206d0)*(b-0.405661288754607d0)*(b-0.5
     #            94338711245393d0)*(b-0.7703023186936794d0)*(b-0.90992
     #            29977317435d0)*b*b*(b-1.d0)*(b-1.d0))
               g(8,kcol+3,i) = 1.d0/((bb-b)*(bb-0.9007700226825652d-1)*
     #            (bb-0.2296976813063206d0)*(bb-0.405661288754607d0)*(bb
     #            -0.594338711245393d0)*(bb-0.7703023186936794d0)*(bb-0.
     #            9099229977317435d0)*bb*bb*(bb-1.d0)*(bb-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-0.9007700226825652d-1)*(x(j)-0.2296976813
     #                  063206d0)*(x(j)-0.405661288754607d0)*(x(j)-0.59
     #                  4338711245393d0)*(x(j)-0.7703023186936794d0)
     #                  *(x(j)-0.9099229977317435d0)*(x(j)-b)*(x(j)-bb)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9007700226825652d-1)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2296976813063206d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.405661288754607d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.594338711245393d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7703023186936794d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9099229977317435d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-bb)
               end do

            else
               aa = -0.2296976813063206d0*(mesh(i)-mesh(i-1))/hx
               a = -0.9007700226825652d-1*(mesh(i)-mesh(i-1))/hx
               gamma1 = -1.d0/aa-1.d0/a-0.24d2
               gamma2 = 1.d0/(1.d0-aa)+1.d0/(1.d0-a)+0.24d2
               h(1,kcol+3,i) = 0.286d3/aa/a
               h(2,kcol+3,i) = 0.286d3/(1.d0-aa)/(1.d0-a)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((aa-0.9007700226825652d-1)*(aa-0.2
     #            296976813063206d0)*(aa-0.405661288754607d0)*(aa-0.594
     #            338711245393d0)*(aa-0.7703023186936794d0)*(aa-0.90992
     #            29977317435d0)*(aa-a)*aa*aa*(aa-1.d0)*(aa-1.d0))
               g(2,kcol+3,i) = 1.d0/((a-0.9007700226825652d-1)*(a-0.229
     #            6976813063206d0)*(a-0.405661288754607d0)*(a-0.5943387
     #            11245393d0)*(a-0.7703023186936794d0)*(a-0.90992299773
     #            17435d0)*(a-aa)*a*a*(a-1.d0)*(a-1.d0))
               g(3,kcol+3,i) = -0.1201314415336355d5/(0.900770022682565
     #            2d-1-a)/(0.9007700226825652d-1-aa)
               g(4,kcol+3,i) = 0.9696023778080566d4/(0.229697681306320
     #            6d0-a)/(0.2296976813063206d0-aa)
               g(5,kcol+3,i) = -0.8929458911121293d4/(0.40566128875460
     #            7d0-a)/(0.405661288754607d0-aa)
               g(6,kcol+3,i) = 0.8929458911121293d4/(0.594338711245393d0
     #            -a)/(0.594338711245393d0-aa)
               g(7,kcol+3,i) = -0.9696023778080565d4/(0.770302318693679
     #            4d0-a)/(0.7703023186936794d0-aa)
               g(8,kcol+3,i) = 0.1201314415336355d5/(0.909922997731743
     #            5d0-a)/(0.9099229977317435d0-aa)

               do j = 1, kcol+3
                  phi = (x(j)-aa)*(x(j)-a)*(x(j)-0.9007700226825652d-1)
     #                  *(x(j)-0.2296976813063206d0)*(x(j)-0.4056612887
     #                  54607d0)*(x(j)-0.594338711245393d0)*(x(j)-0.770
     #                  3023186936794d0)*(x(j)-0.9099229977317435d0)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-aa)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9007700226825652d-1)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2296976813063206d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.405661288754607d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.594338711245393d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7703023186936794d0)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9099229977317435d0)
               end do
            endif
         end do

c-----------------------------------------------------------------------
      elseif (kcol .eq. 9) then
         do i = 1, nint
            hx = (mesh(i+1)-mesh(i))
            if (i .ne. 1 .and. i .ne. nint) then
               a = -0.7229898685756272d-1*(mesh(i)-mesh(i-1))/hx
               b = 1.d0+0.7229898685756272d-1*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -1.d0/a-0.3d2-1.d0/b
               gamma2 = 1.d0/(1.d0-a)+0.3d2+1.d0/(1.d0-b)
               h(1,kcol+3,i) = -0.884d3/a/b
               h(2,kcol+3,i) = 0.884d3/(1.d0-a)/(1.d0-b)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((a-0.7229898685756272d-1)*(a-0.186
     #            3109301186906d0)*(a-0.3341852231986051d0)*(a-0.5d0)*
     #            (a-0.6658147768013949d0)*(a-0.8136890698813094d0)*(a
     #            -0.9277010131424373d0)*(a-b)*a*a*(a-1.d0)*(a-1.d0))
               g(2,kcol+3,i) = 0.4624522420172525d5/(0.722989868575627
     #            2d-1-a)/(0.7229898685756272d-1-b)
               g(3,kcol+3,i) = -0.3688889968687672d5/(0.186310930118690
     #            6d0-a)/(0.1863109301186906d0-b)
               g(4,kcol+3,i) = 0.3332824691372289d5/(0.3341852231986051
     #            d0-a)/(0.3341852231986051d0-b)
               g(5,kcol+3,i) = -0.3232914285714286d5/(0.5d0-a)/(0.5d0-b)
               g(6,kcol+3,i) = 0.3332824691372289d5/(0.665814776801394
     #            9d0-a)/(0.6658147768013949d0-b)
               g(7,kcol+3,i) = -0.3688889968687672d5/(0.813689069881309
     #            4d0-a)/(0.8136890698813094d0-b)
               g(8,kcol+3,i) = 0.4624522420172525d5/(0.927701013142437
     #            3d0-a)/(0.9277010131424373d0-b)
               g(9,kcol+3,i) = 1.d0/((b-a)*(b-0.7229898685756272d-1)*(b
     #            -0.1863109301186906d0)*(b-0.3341852231986051d0)*(b-0.
     #            5d0)*(b-0.6658147768013949d0)*(b-0.8136890698813094d0)
     #            *(b-0.9277010131424373d0)*b*b*(b-1.d0)*(b-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-a)*(x(j)-0.7229898685756272d-1)*(x(j)-0.1
     #                  863109301186906d0)*(x(j)-0.3341852231986051d0)
     #                  *(x(j)-0.5d0)*(x(j)-0.6658147768013949d0)*(x(j)
     #                  -0.8136890698813094d0)*(x(j)-0.927701013142437
     #                  3d0)*(x(j)-b)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7229898685756272d-1)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1863109301186906d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3341852231986051d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6658147768013949d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8136890698813094d0)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9277010131424373d0)
                  g(9,j,i) = g(9,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
               end do

            elseif (i .eq. 1) then
               b = 1.d0+0.7229898685756272d-1*(mesh(i+2)-mesh(i+1))/hx
               bb = 1.d0+0.1863109301186906d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -0.3d2-1.d0/b-1.d0/bb
               gamma2 = 0.3d2+1.d0/(1.d0-b)+1.d0/(1.d0-bb)
               h(1,kcol+3,i) = -0.884d3/b/bb
               h(2,kcol+3,i) = 0.884d3/(1.d0-b)/(1.d0-bb)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 0.4624522420172525d5/(0.722989868575627
     #            2d-1-bb)/(0.7229898685756272d-1-b)
               g(2,kcol+3,i) = -0.3688889968687672d5/(0.186310930118690
     #            6d0-bb)/(0.1863109301186906d0-b)
               g(3,kcol+3,i) = 0.3332824691372289d5/(0.3341852231986051
     #            d0-bb)/(0.3341852231986051d0-b)
               g(4,kcol+3,i) = -0.3232914285714286d5/(0.5d0-bb)
     #            /(0.5d0-b)
               g(5,kcol+3,i) = 0.3332824691372289d5/(0.665814776801394
     #            9d0-bb)/(0.6658147768013949d0-b)
               g(6,kcol+3,i) = -0.3688889968687672d5/(0.813689069881309
     #            4d0-bb)/(0.8136890698813094d0-b)
               g(7,kcol+3,i) = 0.4624522420172525d5/(0.927701013142437
     #            3d0-bb)/(0.9277010131424373d0-b)
               g(8,kcol+3,i) = 1.d0/((b-bb)*(b-0.7229898685756272d-1)*(b
     #            -0.1863109301186906d0)*(b-0.3341852231986051d0)*(b-0.
     #            5d0)*(b-0.6658147768013949d0)*(b-0.8136890698813094d0)
     #            *(b-0.9277010131424373d0)*b*b*(b-1.d0)*(b-1.d0))
               g(9,kcol+3,i) = 1.d0/((bb-b)*(bb-0.7229898685756272d-1)*
     #            (bb-0.1863109301186906d0)*(bb-0.3341852231986051d0)*
     #            (bb-0.5d0)*(bb-0.6658147768013949d0)*(bb-0.8136890698
     #            813094d0)*(bb-0.9277010131424373d0)*bb*bb*(bb-1.d0)
     #            *(bb-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-0.7229898685756272d-1)*(x(j)-0.1
     #                  863109301186906d0)*(x(j)-0.3341852231986051d0)
     #                  *(x(j)-0.5d0)*(x(j)-0.6658147768013949d0)*(x(j)
     #                  -0.8136890698813094d0)*(x(j)-0.927701013142437
     #                  3d0)*(x(j)-b)*(x(j)-bb)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7229898685756272d-1)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1863109301186906d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3341852231986051d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6658147768013949d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8136890698813094d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9277010131424373d0)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
                  g(9,j,i) = g(9,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-bb)
               end do

            else
               aa = -0.1863109301186906d0*(mesh(i)-mesh(i-1))/hx
               a = -0.7229898685756272d-1*(mesh(i)-mesh(i-1))/hx
               gamma1 = -1.d0/aa-1.d0/a-0.3d2
               gamma2 = 1.d0/(1.d0-aa)+1.d0/(1.d0-a)+0.3d2
               h(1,kcol+3,i) = -0.884d3/aa/a
               h(2,kcol+3,i) = 0.884d3/(1.d0-aa)/(1.d0-a)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((aa-0.7229898685756272d-1)*(aa-0.1
     #            863109301186906d0)*(aa-0.3341852231986051d0)*(aa-0.
     #            5d0)*(aa-0.6658147768013949d0)*(aa-0.813689069881309
     #            4d0)*(aa-0.9277010131424373d0)*(aa-a)*aa*aa*(aa-1.d0)
     #            *(aa-1.d0))
               g(2,kcol+3,i) = 1.d0/((a-0.7229898685756272d-1)*(a-0.186
     #            3109301186906d0)*(a-0.3341852231986051d0)*(a-0.5d0)*
     #            (a-0.6658147768013949d0)*(a-0.8136890698813094d0)*(a
     #            -0.9277010131424373d0)*(a-aa)*a*a*(a-1.d0)*(a-1.d0))
               g(3,kcol+3,i) = 0.4624522420172525d5/(0.722989868575627
     #            2d-1-a)/(0.7229898685756272d-1-aa)
               g(4,kcol+3,i) = -0.3688889968687672d5/(0.186310930118690
     #            6d0-a)/(0.1863109301186906d0-aa)
               g(5,kcol+3,i) = 0.3332824691372289d5/(0.3341852231986051
     #            d0-a)/(0.3341852231986051d0-aa)
               g(6,kcol+3,i) = -0.3232914285714286d5/(0.5d0-a)
     #            /(0.5d0-aa)
               g(7,kcol+3,i) = 0.3332824691372289d5/(0.665814776801394
     #            9d0-a)/(0.6658147768013949d0-aa)
               g(8,kcol+3,i) = -0.3688889968687672d5/(0.813689069881309
     #            4d0-a)/(0.8136890698813094d0-aa)
               g(9,kcol+3,i) = 0.4624522420172525d5/(0.927701013142437
     #            3d0-a)/(0.9277010131424373d0-aa)

               do j = 1, kcol+3
                  phi = (x(j)-aa)*(x(j)-a)*(x(j)-0.7229898685756272d-1)
     #                  *(x(j)-0.1863109301186906d0)*(x(j)-0.3341852231
     #                  986051d0)*(x(j)-0.5d0)*(x(j)-0.665814776801394
     #                  9d0)*(x(j)-0.8136890698813094d0)
     #                  *(x(j)-0.9277010131424373d0)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-aa)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7229898685756272d-1)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1863109301186906d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3341852231986051d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6658147768013949d0)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8136890698813094d0)
                  g(9,j,i) = g(9,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9277010131424373d0)
               end do
            endif
         end do

c-----------------------------------------------------------------------
      elseif (kcol .eq. 10) then
         do i = 1, nint
            hx = (mesh(i+1)-mesh(i))
            if (i .ne. 1 .and. i .ne. nint) then
               a = -0.5929571219129399d-1*(mesh(i)-mesh(i-1))/hx
               b = 1.d0+0.5929571219129399d-1*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -1.d0/a-0.3666666666666667d2-1.d0/b
               gamma2 = 1.d0/(1.d0-a)+0.3666666666666667d2+1.d0
     #            /(1.d0-b)
               h(1,kcol+3,i) = 0.2799333333333333d4/a/b
               h(2,kcol+3,i) = 0.2799333333333333d4/(1.d0-a)/(1.d0-b)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((a-0.5929571219129399d-1)*(a-0.1539
     #            696908715823d0)*(a-0.2792835119457421d0)*(a-0.4241841
     #            678533667d0)*(a-0.5758158321466333d0)*(a-0.7207164880
     #            542579d0)*(a-0.8460303091284177d0)*(a-0.9407042878087
     #            06d0)*(a-b)*a*a*(a-1.d0)*(a-1.d0))
               g(2,kcol+3,i) = -0.1785201442945021d6/(0.592957121912939
     #            9d-1-a)/(0.5929571219129399d-1-b)
               g(3,kcol+3,i) = 0.1412225097382635d6/(0.153969690871582
     #            3d0-a)/(0.1539696908715823d0-b)
               g(4,kcol+3,i) = -0.1259240998337037d6/(0.279283511945742
     #            1d0-a)/(0.2792835119457421d0-b)
               g(5,kcol+3,i) = 0.1197516586185479d6/(0.424184167853366
     #            7d0-a)/(0.4241841678533667d0-b)
               g(6,kcol+3,i) = -0.1197516586185479d6/(0.575815832146633
     #            3d0-a)/(0.5758158321466333d0-b)
               g(7,kcol+3,i) = 0.1259240998337037d6/(0.720716488054257
     #            9d0-a)/(0.7207164880542579d0-b)
               g(8,kcol+3,i) = -0.1412225097382635d6/(0.846030309128417
     #            7d0-a)/(0.8460303091284177d0-b)
               g(9,kcol+3,i) = 0.1785201442945021d6/(0.94070428780870
     #            6d0-a)/(0.940704287808706d0-b)
               g(10,kcol+3,i) = 1.d0/((b-a)*(b-0.5929571219129399d-1)*(b
     #            -0.1539696908715823d0)*(b-0.2792835119457421d0)*(b-0.
     #            4241841678533667d0)*(b-0.5758158321466333d0)*(b-0.720
     #            7164880542579d0)*(b-0.8460303091284177d0)*(b-0.940704
     #            287808706d0)*b*b*(b-1.d0)*(b-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-a)*(x(j)-0.5929571219129399d-1)*(x(j)-0.1
     #               539696908715823d0)*(x(j)-0.2792835119457421d0)
     #               *(x(j)-0.4241841678533667d0)*(x(j)-0.5758158321466
     #               333d0)*(x(j)-0.7207164880542579d0)*(x(j)-0.8460303
     #               091284177d0)*(x(j)-0.940704287808706d0)*(x(j)-b)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5929571219129399d-1)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1539696908715823d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2792835119457421d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.4241841678533667d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5758158321466333d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7207164880542579d0)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8460303091284177d0)
                  g(9,j,i) = g(9,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.940704287808706d0)
                  g(10,j,i) = g(10,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
               end do

            elseif (i .eq. 1) then
               b = 1.d0+0.5929571219129399d-1*(mesh(i+2)-mesh(i+1))/hx
               bb = 1.d0+0.1539696908715823d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -0.3666666666666667d2-1.d0/b-1.d0/bb
               gamma2 = 0.3666666666666667d2+1.d0/(1.d0-b)
     #            +1.d0/(1.d0-bb)
               h(1,kcol+3,i) = 0.2799333333333333d4/b/bb
               h(2,kcol+3,i) = 0.2799333333333333d4/(1.d0-b)/(1.d0-bb)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = -0.1785201442945021d6/(0.592957121912939
     #            9d-1-bb)/(0.5929571219129399d-1-b)
               g(2,kcol+3,i) = 0.1412225097382635d6/(0.153969690871582
     #            3d0-bb)/(0.1539696908715823d0-b)
               g(3,kcol+3,i) = -0.1259240998337037d6/(0.279283511945742
     #            1d0-bb)/(0.2792835119457421d0-b)
               g(4,kcol+3,i) = 0.1197516586185479d6/(0.424184167853366
     #            7d0-bb)/(0.4241841678533667d0-b)
               g(5,kcol+3,i) = -0.1197516586185479d6/(0.575815832146633
     #            3d0-bb)/(0.5758158321466333d0-b)
               g(6,kcol+3,i) = 0.1259240998337037d6/(0.720716488054257
     #            9d0-bb)/(0.7207164880542579d0-b)
               g(7,kcol+3,i) = -0.1412225097382635d6/(0.846030309128417
     #            7d0-bb)/(0.8460303091284177d0-b)
               g(8,kcol+3,i) = 0.1785201442945021d6/(0.94070428780870
     #            6d0-bb)/(0.940704287808706d0-b)
               g(9,kcol+3,i) = 1.d0/((b-bb)*(b-0.5929571219129399d-1)*(b
     #            -0.1539696908715823d0)*(b-0.2792835119457421d0)*(b-0.
     #            4241841678533667d0)*(b-0.5758158321466333d0)*(b-0.720
     #            7164880542579d0)*(b-0.8460303091284177d0)*(b-0.940704
     #            287808706d0)*b*b*(b-1.d0)*(b-1.d0))
               g(10,kcol+3,i) = 1.d0/((bb-b)*(bb-0.5929571219129399d-1)
     #            *(bb-0.1539696908715823d0)*(bb-0.2792835119457421d0)
     #            *(bb-0.4241841678533667d0)*(bb-0.5758158321466333d0)
     #            *(bb-0.7207164880542579d0)*(bb-0.8460303091284177d0)
     #            *(bb-0.940704287808706d0)*bb*bb*(bb-1.d0)*(bb-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-0.5929571219129399d-1)*(x(j)-0.1539696908
     #               715823d0)*(x(j)-0.2792835119457421d0)*(x(j)-0.4241
     #               841678533667d0)*(x(j)-0.5758158321466333d0)*(x(j)
     #               -0.7207164880542579d0)*(x(j)-0.8460303091284177d0)
     #               *(x(j)-0.940704287808706d0)*(x(j)-b)*(x(j)-bb)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5929571219129399d-1)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1539696908715823d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2792835119457421d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.4241841678533667d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5758158321466333d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7207164880542579d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8460303091284177d0)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.940704287808706d0)
                  g(9,j,i) = g(9,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
                  g(10,j,i) = g(10,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-bb)
               end do

            else
               aa = -0.1539696908715823d0*(mesh(i)-mesh(i-1))/hx
               a = -0.5929571219129399d-1*(mesh(i)-mesh(i-1))/hx
               gamma1 = -1.d0/aa-1.d0/a-0.3666666666666667d2
               gamma2 = 1.d0/(1.d0-aa)+1.d0/(1.d0-a)
     #            +0.3666666666666667d2
               h(1,kcol+3,i) = 0.2799333333333333d4/aa/a
               h(2,kcol+3,i) = 0.2799333333333333d4/(1.d0-aa)/(1.d0-a)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((aa-0.5929571219129399d-1)*(aa-0.1
     #            539696908715823d0)*(aa-0.2792835119457421d0)*(aa-0.42
     #            41841678533667d0)*(aa-0.5758158321466333d0)*(aa-0.720
     #            7164880542579d0)*(aa-0.8460303091284177d0)*(aa-0.9407
     #            04287808706d0)*(aa-a)*aa*aa*(aa-1.d0)*(aa-1.d0))
               g(2,kcol+3,i) = 1.d0/((a-0.5929571219129399d-1)*(a-0.1539
     #            696908715823d0)*(a-0.2792835119457421d0)*(a-0.4241841
     #            678533667d0)*(a-0.5758158321466333d0)*(a-0.7207164880
     #            542579d0)*(a-0.8460303091284177d0)*(a-0.9407042878087
     #            06d0)*(a-aa)*a*a*(a-1.d0)*(a-1.d0))
               g(3,kcol+3,i) = -0.1785201442945021d6/(0.592957121912939
     #            9d-1-a)/(0.5929571219129399d-1-aa)
               g(4,kcol+3,i) = 0.1412225097382635d6/(0.153969690871582
     #            3d0-a)/(0.1539696908715823d0-aa)
               g(5,kcol+3,i) = -0.1259240998337037d6/(0.279283511945742
     #            1d0-a)/(0.2792835119457421d0-aa)
               g(6,kcol+3,i) = 0.1197516586185479d6/(0.424184167853366
     #            7d0-a)/(0.4241841678533667d0-aa)
               g(7,kcol+3,i) = -0.1197516586185479d6/(0.575815832146633
     #            3d0-a)/(0.5758158321466333d0-aa)
               g(8,kcol+3,i) = 0.1259240998337037d6/(0.720716488054257
     #            9d0-a)/(0.7207164880542579d0-aa)
               g(9,kcol+3,i) = -0.1412225097382635d6/(0.846030309128417
     #            7d0-a)/(0.8460303091284177d0-aa)
               g(10,kcol+3,i) = 0.1785201442945021d6/(0.94070428780870
     #            6d0-a)/(0.940704287808706d0-aa)

               do j = 1, kcol+3
                  phi = (x(j)-aa)*(x(j)-a)*(x(j)-0.5929571219129399d-1)
     #               *(x(j)-0.1539696908715823d0)*(x(j)-0.2792835119457
     #               421d0)*(x(j)-0.4241841678533667d0)*(x(j)-0.5758158
     #               321466333d0)*(x(j)-0.7207164880542579d0)*(x(j)
     #               -0.8460303091284177d0)*(x(j)-0.940704287808706d0)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-aa)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5929571219129399d-1)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1539696908715823d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2792835119457421d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.4241841678533667d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5758158321466333d0)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7207164880542579d0)
                  g(9,j,i) = g(9,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8460303091284177d0)
                  g(10,j,i) = g(10,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.940704287808706d0)
               end do
            endif
         end do
      endif

      return
      end

      subroutine fdderivf(t, x, u, ux, uxx, dfdu, dfdux, dfduxx, npde,
     &                    f, work)
c-----------------------------------------------------------------------
c         This subroutine produces central differences approximations 
c         for this first and second partial derivitives of the the PDE 
c         system f using in BACOLI models, ut = f(t, x, u, ux, uxx). 
c         This subroutine is called by CALJAC in place of the user
c         provided routine DERIVF.
c-----------------------------------------------------------------------

c     Input:
      integer          npde
      double precision t
      double precision x
      double precision u(npde)
      double precision ux(npde)
      double precision uxx(npde)
      external         f
c     Output:
      double precision dfdu(npde,npde)
      double precision dfdux(npde,npde)
      double precision dfduxx(npde,npde)
c     Locals:
      external         d1mach
      double precision del, d1, d2, delinv, oldval
c     Work storage:
      double precision work(npde,2)
c     Machine dependent variables:
      double precision squr, uround, d1mach
c     Loop indices
      integer i, j
c----------------------------------------------------------------------

c     We apply increments on the order of size sqrt uround

c     uround is the largest relative spacing on this machine, this is
c     the smallest double precision floating point number possible on
c     machine.
      uround = d1mach(4)

      squr = sqrt(uround)

c     first, dfdu
      do 10 j = 1, npde
          oldval = u(j)
          del = squr * max(abs(u(j)), uround)

c         d1 and d2 are what is *actually* added or subtracted
          d1 = (u(j) + del) - u(j)
          d2 = (u(j) - del) - u(j)
          u(j) = u(j) + d1
          call f(t, x, u, ux, uxx, work(1,1), npde)
          u(j) = oldval + d2
          call f(t, x, u, ux, uxx, work(1,2), npde)
          u(j) = oldval
          delinv = 1d0/(d1-d2)
          do 15 i = 1, npde
              dfdu(i,j) = (work(i,1) - work(i,2)) * delinv
   15     continue
   10 continue

c     second, dfdux
      do 20 j = 1, npde
          oldval = ux(j)
          del = squr * max(abs(ux(j)), uround)

c         d1 and d2 are what is *actually* added or subtracted
          d1 = (ux(j) + del) - ux(j)
          d2 = (ux(j) - del) - ux(j)
          ux(j) = ux(j) + d1
          call f(t, x, u, ux, uxx, work(1,1), npde)
          ux(j) = oldval + d2
          call f(t, x, u, ux, uxx, work(1,2), npde)
          ux(j) = oldval
          delinv = 1d0/(d1-d2)
          do 25 i = 1, npde
              dfdux(i,j) = (work(i,1) - work(i,2)) * delinv
   25     continue
   20 continue

c     last, dfduxx 
      do 30 j = 1, npde
          oldval = uxx(j)
          del = squr * max(abs(uxx(j)), uround)

c         d1 and d2 are what is *actually* added or subtracted
          d1 = (uxx(j) + del) - uxx(j)
          d2 = (uxx(j) - del) - uxx(j)
          uxx(j) = uxx(j) + d1
          call f(t, x, u, ux, uxx, work(1,1), npde)
          uxx(j) = oldval + d2
          call f(t, x, u, ux, uxx, work(1,2), npde)
          uxx(j) = oldval
          delinv = 1d0/(d1-d2)
          do 35 i = 1, npde
              dfduxx(i,j) = (work(i,1) - work(i,2)) * delinv
   35     continue
   30 continue

      return
      end

      subroutine fdbndxri(t, u, ux, dbdu, dbdux, dbdt, npde, bndx, work)
c-----------------------------------------------------------------------
c         This subroutine produces central differences approximations
c         for the first partial derivitives of the left and right
c         boundary conditions using in BACOLI models, b(t,u,ux) = 0.
c         This subroutine is called by INIY in place of
c         user subroutines DIFBXA and DIFBXB in order to simplify
c         the user interface when finite difference approximated Jacobi
c         are used (ie, when mflag(6) = 0).
c-----------------------------------------------------------------------

c     Input:
      integer          npde
      double precision t
      double precision u(npde)
      double precision ux(npde)
      external         bndx
c     Output:
      double precision dbdu(npde,npde)
      double precision dbdux(npde,npde)
      double precision dbdt(npde)
c     Locals:
      external         d1mach
      double precision del, d1, d2, delinv, oldval
c     Work storage:
      double precision work(npde,2)
c     Machine dependent variables:
      double precision squr, uround, d1mach
c     Loop indices
      integer i, j
c-----------------------------------------------------------------------

c     We apply increments on the order of size sqrt uround

c     uround is the largest relative spacing on this machine, this is 
c     the smallest double precision floating point number possible on
c     this machine.
      uround = d1mach(4)

      squr = sqrt(uround)

c     first, dbdu
      do 10 j = 1, npde
          oldval = u(j)
          del = squr * max(abs(u(j)), uround)

c         d1 and d2 are what is *actually* added or subtracted
          d1 = (u(j) + del) - u(j)
          d2 = (u(j) - del) - u(j)
          u(j) = u(j) + d1
          call bndx(t, u, ux, work(1,1), npde)
          u(j) = oldval + d2
          call bndx(t, u, ux, work(1,2), npde)
          u(j) = oldval
          delinv = 1d0/(d1-d2)
          do 15 i = 1, npde
              dbdu(i,j) = (work(i,1) - work(i,2)) * delinv
   15     continue
   10 continue

c     second, dbdux
      do 20 j = 1, npde
          oldval = ux(j)
          del = squr * max(abs(ux(j)), uround)

c         d1 and d2 are what is *actually* added or subtracted
          d1 = (ux(j) + del) - ux(j)
          d2 = (ux(j) - del) - ux(j)
          ux(j) = ux(j) + d1
          call bndx(t, u, ux, work(1,1), npde)
          ux(j) = oldval + d2
          call bndx(t, u, ux, work(1,2), npde)
          ux(j) = oldval
          delinv = 1d0/(d1-d2)
          do 25 i = 1, npde
              dbdux(i,j) = (work(i,1) - work(i,2)) * delinv
   25     continue
   20 continue

c     last, dbdt
      del = squr * max(abs(t), uround)
c     d1 and d2 are what is *actually* added or subtracted
      d1 = (t + del) - t
      d2 = (t - del) - t
      call bndx(t+d1, u, ux, work(1,1), npde)
      call bndx(t+d2, u, ux, work(1,2), npde)
      delinv = 1d0/(d1-d2)
      do 30 i = 1, npde
          dbdt(i) = (work(i,1) - work(i,2)) * delinv
   30 continue

      return
      end