! Program grid_implicit_sparse.
!
! Program written by María Martínez-Barbeito and Pere Colet.
!
! Integrates power grid dynamics using a first order semi-implicit Euler method:
! x(t+dt)=x(t)+dt*F(x(t+dt)) approximated as x(t+dt)=x(t)+dx with dx solution of
! the linear set of equations (I-dt*J)*dx=F(x)*dt where J is the Jacobian.
!
! Takes advantage of the sparse structure of the grid (<2% non-zero terms in I-dtJ).
! The program uses PARDISO routines included in Intel Math Kernel Libraries (MKL) to
! solve the sparse linear set of equations. There is also an open source version as
! part of the PARDISO project: https://www.pardiso-project.org/
!
!-------------------------------------------------------------------------------
!
! Power grid with N nodes and L links
!
! A node can be a generator or a consumer.
! A link is a transmission line connecting two nodes (undirected network).
!
! POWER PLANTS: 4 equations and 4 variables (phase,freq.,Pm,Ps)
! CONSUMERS: 2 equations and 2 variables (phase,freq.)
!
! Variables stored in one vector x:
! x = ( phase_1, ..., phase_N, freq_1, ..., freq_N,
!       Pm_1, ..., Pm_Nplants, Ps_1, ..., Ps_Nplants )
! Pm_k is the mechanical power for plant on node plantNode(k).
!
! Parameters read from parameters.dat
! Inital condition read from file or given by the dispatch.
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
program grid_implicit_sparse
implicit none
!
integer :: N,Nplants,L,Nt,Nassets ! number of: nodes, power plants, transmission lines, variables, assets
integer :: Ndspch,nstep,slowOps ! number of: dispatches, time steps within a dispatch, slow switchOffs
integer :: i,j,k,idspch,istep   ! counters
!
double precision :: twoPi ! 2*pi number (=2*3.14)
double precision :: t0,t,dt,twrite,writeStep !initial time, time, integration step, time of last writing, writing step
double precision :: wR,D,DwR    ! reference frequecy, freq. dependend load, D/wR.
double precision :: start,finish     !cpu time counters
double precision :: epsilon,tau_ou,p_ou,coeff_ou  ! noise variables
!
double precision :: P_wind,P_wind0 ! wind power
!
integer, dimension(:,:), allocatable :: Ml   ! initial and final nodes of a line
double precision, dimension(:), allocatable :: B     ! line capacity
double precision, dimension(:), allocatable :: itau  !time scale in dPm/dt
double precision, dimension(:), allocatable :: kwR,kwRon,kwRoff,konoff,kwRoff0,kwRstp  !secondary control
double precision, dimension(:), allocatable :: PGRwR,PGRwRon,PGRwRoff,PGRwRoff0,PGRwRstp  !primary control
double precision, dimension(:), allocatable :: wR22HPG,wR22HPGon,wR22HPGoff,wR22HPGoff0,wR22HPGstp !inertia
double precision, dimension(:), allocatable :: on_off ! which plants are working in the reference case
!
double precision, dimension(:), allocatable :: lambda,Pref    ! Power plants forcing parameters
character (len=10), dimension(:), allocatable  :: plantStatus  !status of power plants
integer, dimension(:), allocatable :: plantNode,plantId ! Node where a plant is located and Id of plant located at a given node
integer, dimension(:), allocatable :: slowOpId ! plant id for which a slow switch off is taking place
!
double precision, dimension(:), allocatable :: Pl,dPl,Pleff,xi_ou,Gaussian  ! loads
double precision, dimension(:), allocatable :: P_assets,P_assets0,dP_assets0 ! asset power and change on a time step.
double precision, dimension(:), allocatable :: P_basket ! power from assets (including wind), line balance and primary control
integer, dimension (:), allocatable :: assetNode  ! node where asset is located
character (len=10), dimension(:), allocatable  :: assetKind  ! type of asset
!
double precision, dimension(:), allocatable :: x  ! vector of variables
double precision, dimension(:), allocatable :: deltax,deltaxs  ! F*dt and dx given by pardiso
double precision, dimension(:), allocatable :: IdtJ  ! nonzero coefficients of I-dt*J by rows.
integer, dimension(:), allocatable :: idxcol  ! col index: IdtJ(k) belongs to column idxcol(k) of I-dt*J
integer, dimension(:), allocatable :: idxrow  ! row index: VIJ(idxrow(i))) is the first element in row i of I-dt*J
integer, dimension(:), allocatable :: idxFreqEqPhase,idxFreqEqFreq  ! idxFreqEqPhase(i) and idxFreqEqFreq(i) locate in IdtJ partial freq-eq_i respect to phase_i and freq_i
integer, dimension(:), allocatable :: idxLineA,idxLineB !for line l connecting i and j, idxLineA(l) locates in IdtJ partial freq-eq_i respect to phase_j and idxLineB(l) partial freq-eq_j respect to phase_i
!
character (len=120) :: initFile
character (len=120) :: dispFile
character (len=120) :: resFile
!
include 'mkl_pardiso.fi'   ! definitions needed by MKL pardiso.
TYPE(MKL_PARDISO_HANDLE) pt(64)   ! internal solver memory pointer (integer*8 vector with 64 elements).
integer :: maxfct, mnum, mtype, phase, error, msglvl,nrhs
integer, dimension (64) :: iparm  ! array of parameters for pardiso options
integer :: nzcoeffs    ! number of non-zero coefficients of I-dtJ
integer :: idum(1)
double precision ::  ddum(1)

!-------------------------------------------------------------------------------
!
call cpu_time(start)
!
twoPi=8.d0*datan(1.d0)
!
!-------------------------------------------------------------------------------
!
call parameters_pardiso   ! set pardiso default parameters
!
call read_and_allocate
!
call build_sparse
!
call initialization
!
!call write_data
!
!------------------------------- Dynamical evolution ---------------------------
!
do idspch=1,Ndspch  ! OUTER LOOP over dispatches with data assimilation.
  call dispatch
  phase=13   ! Pardiso has 3 phases: analysis, factorization, and solving. phase=13 means go from analysis to solving. Analysis is slow and only neded for the first call with a given matrix. After it is fine do only factorization+solving (phase=23). Here analysis is done every dispatch.
!
  do istep=1,nstep  ! INNER LOOP over dt steps during a dispatch interval.
!   use load and assets power at t+dt. Thus, Pl and P_assets are updated before dynamics
    P_assets0=P_assets0+dP_assets0
    Pl=Pl+dPl
!
    ! If wind generation exceeds demand, we only introduce the demanded amount (-solar) in the system
    P_wind0=0.0d0 !scalar
    !Calculate total wind generation
    do i=1,Nassets
      if (assetKind(i).eq.'wind') then
        P_wind0=P_wind0+P_assets0(i)
      endif
    enddo !i
    P_wind=min(P_wind0,sum(Pl)-sum(P_assets0)+P_wind0)
    P_assets=P_assets0
    do i=1,Nassets
      if (assetKind(i).eq.'wind') then
        P_assets(i)=P_wind*P_assets0(i)/P_wind0
      endif
    enddo
!
    call rand_gaussian
    xi_ou=p_ou*xi_ou+coeff_ou*Gaussian   ! Update Orstein-Uhlenbeck noise
    Pleff=Pl*(1.d0+epsilon*xi_ou) ! Fluctuating load
!
    do k=1,slowOps  ! Slow switch-off of power plants
      call plantOp('SlowSwitchOff',slowOpId(k))
    enddo !k
!
    call dynamics  !calculate deltax and terms of I-dt*J depending on variables.
!
!   solve (I-dt*J).deltaxs = deltax
    call pardiso (pt,maxfct,mnum,mtype,phase,Nt,IdtJ,idxrow,idxcol,idum,nrhs,iparm,msglvl,deltax,deltaxs,error)
    phase=23   ! after first call to pardiso, we skip analysis.
!
    x=x+deltaxs  ! Update variables
    t=t+dt       ! Update time
!
    if (t.ge.twrite+writeStep*0.999999) call write_data
  enddo !istep
!
enddo !idspch
!
!-------------------------------------------------------------------------------
!
call cpu_time(finish)
write(*,*) 'Finished.', Ndspch*nstep, 'time-steps performed in',real(finish-start), 'seconds.'
!
!-------------------------------------------------------------------------------

contains

!------------------ pardiso parameters setting ---------------------------------
!
subroutine parameters_pardiso
  implicit none
! Values can be overwritten in parameters1.dat
! Parameters 7,9,12,15-17 and 21-64 not in use.
!
  iparm = 0  ! Set to 0 all but those set below or in parameters1.dat
  iparm(1) = 1 ! no solver default
  iparm(2) = 2 ! fill-in reordering from METIS
  iparm(3) = 1 ! number of cores          <------------------------
  iparm(4) = 81 ! Use CGS iteration with precision 10^-8.   <------
  iparm(5) = 0 ! no user fill-in reducing permutation
  iparm(6) = 0 ! =0 solution on the first n components of x
  iparm(8) = 2 ! numbers of iterative refinement steps      <------
  iparm(10) = 13 ! perturb the pivot elements with 1E-13
  iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
  iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric)
  iparm(14) = 0 ! Output: number of perturbed pivots
  iparm(18) = 0 ! Output: number of nonzeros in factor LU (reported if negative)
  iparm(19) = 0 ! Output: Mflops for LU factorization (reported if negative)
  iparm(20) = 0 ! Output: Numbers of CG Iterations
  error = 0 ! initialize error flag
  msglvl = 0 ! If 1 prints statistical information, If 0 does not.
  mtype = 11 ! Type of matrix. 11 corresponds to real unsymmetric
  maxfct=1   ! Number factorizations of identical sparse structures used simultaneously the same handle pt.
  mnum=1     ! structure used among maxfct structures with the same handle.
  nrhs=1     ! number of rhs to be solved. We only solve for one rhs, deltax
!  Initialize internal memory pointer. Only needed for FIRST call to PARDISO.
!  do i = 1, 64
!    pt(i)%DUMMY = 0
!  enddo
  return
end subroutine parameters_pardiso
!
!-------------------------------------------------------------------------------

!------------------ Construction of sparse structure for I-dt*J ----------------
!
subroutine build_sparse
!  Builts IdtJ structure, a vector with the non-zero coefficients of I-dt*J in format CSR3 (3-array variation of the Compressed Sparse Row format)
!  Evaluates vectors idxcol and idxrow and fills the constant terms of IdtJ.
!
  implicit none
  integer, dimension (:), allocatable :: link !vector to temporarily store lines that connect to a node.
  integer :: ni,nj,il,idx
  allocate (link(N))
!
! 1. Phase eq.: rows 1 to N of I-dt*J. Two non-zero constant terms per row: diagonal and partial phase-eq_i respect to freq_i.
  do i=1,N
    idxrow(i)=2*i-1  ! location in IdtJ of the first term of row i of I-dt*J
    IdtJ(2*i-1)=1.d0 ! diagonal term of I-dt*J
    idxcol(2*i-1)=i  ! I-dt*J column containing the diagonal term.
    IdtJ(2*i)=-dt    ! partial phase-eq_i respect to freq_i
    idxcol(2*i)=N+i  ! I-dt*J column containing partial phase-eq_i respect to freq_i
  enddo
  idx=2*N+1          ! location in IdtJ to start filling the next coefficients.
!
! 2. Frequency eq.: rows N+1 to 2*N of I-dt*J. Non-zero terms are the diagonal, those related to lines and, for power plants, partial respect to Pm. Evaluate only indexes. Values filled in subroutine dynamics.
  do i=1,N
    idxrow(i+N)=idx  ! location in IdtJ of the first term of row i+N of I-dt*J
    link=0
    do il=1,L        ! Go through list of transmission lines
      ni=Ml(il,1)    ! Transmission line il connects node ni...
      nj=Ml(il,2)    ! ... with node nj
      if (ni.eq.i) link(nj)=il   ! if line il starts at i, link stores il.
      if (nj.eq.i) link(ni)=-il  ! if line il ends at i, link stores -il.
    enddo
    do j=1,i-1       ! lines connecting node i with node j<i
      if (link(j).eq.0) cycle
        if (link(j).gt.0) then   ! Store location of partial freq-eq_i respect to phase_j in IdtJ...
          idxLineA(link(j))=idx  ! ...in idxLineA if line starts at i...
        else
          idxLineB(-link(j))=idx ! ...or in idxLineB if line ends at i.
        endif
        idxcol(idx)=j  ! I-dt*J column containing partial freq-eq_i respect to phase_j
        idx=idx+1
    enddo
    idxcol(idx)=i      ! I-dt*J column containing partial freq-eq_i respect to phase_i
    idxFreqEqPhase(i)=idx  ! location in IdtJ of partial freq-eq_i respect to phase_i
    idx=idx+1
    do j=i+1,N         ! lines connecting node i with node j>i
      if (link(j).eq.0) cycle
        if (link(j).gt.0) then   ! Store location of partial freq-eq_i respect to phase_j in IdtJ...
          idxLineA(link(j))=idx  ! ...in idxLineA if line starts at i...
        else
          idxLineB(-link(j))=idx ! ...or in idxLineB if line end at i.
        endif
        idxcol(idx)=j  ! I-dt*J column containing partial freq-eq_i respect to phase_j
        idx=idx+1
    enddo
    idxcol(idx)=i+N    ! I-dt*J column containing partial freq-eq_i respect to freq_i
    idxFreqEqFreq(i)=idx  ! location in IdtJ of partial freq-eq_i respect to freq_i
    idx=idx+1
    if (plantId(i).gt.0) then
      idxcol(idx)=2*N+plantId(i) ! I-dt*J column containing partial freq-eq_i respect to Pm_i
      idx=idx+1
    endif
  enddo
!
! 3. Pm equation for conventional power plants: rows 2N+1 to 2N+Nplants of I-dt*J. The non-zero terms are partial Pm-eq_k respect to freq_plantNode(k), the diagonal and partial Pm-eq_k respect to Ps_k. The last two are constant and filled here. The first is filled by plantOp subrutine.
  do k=1,Nplants
    idxrow(2*N+k)=idx            ! location in IdtJ of the first term of row 2*N+k of I-dt*J
    idxcol(idx)=N+plantNode(k)   ! I-dt*J column containing partial Pm-eq_k respect to freq_plantNode(k)
    IdtJ(idx+1)=1.d0+dt*itau(k)  ! diagonal term of I-dt*J including partial Pm-eq_k respect to Pm_k
    idxcol(idx+1)=2*N+k          ! I-dt*J column containing partial Pm-eq_k respect to Pm_k
    IdtJ(idx+2)=-dt*itau(k)      ! partial Pm-eq_k respect to Ps_k
    idxcol(idx+2)=2*N+Nplants+k  ! I-dt*J column containing partial Pm-eq_k respect to Ps_k
    idx=idx+3
  enddo
!
! 4. Ps equation for conventional power plants: rows 2N+Nplants+1 to 2N+2Nplants of I-dt*J. The non-zero terms are partial Ps-eq_k respect to freq_plantNode(k) and the diagonal. The last is constant and filled here. The first is filled by plantOp subrutine.
  do k=1,Nplants
    idxrow(2*N+Nplants+k)=idx     ! location in IdtJ of the first term of row 2*N+Nplants+k of I-dt*J
    idxcol(idx)=N+plantNode(k)    ! I-dt*J column containing partial Ps-eq_k respect to freq_plantNode(k)
    IdtJ(idx+1)=1.d0+dt*lambda(k) ! diagonal term of I-dt*J including partial Ps-eq_k respect to Ps_k
    idxcol(idx+1)=2*N+Nplants+k   ! I-dt*J column containing partial Ps-eq_k respect to Ps_k
    idx=idx+2
  enddo
  idxrow(Nt+1)=idx  ! the last term of idxrow(i) is the number of nonzero elements + 1.
!
! Check if the number of non-zero elements is correct
  if (idx.ne.nzcoeffs+1) then
    print *, 'Error counting non-zero elements.'; STOP
  endif
  deallocate (link)
  return
end subroutine build_sparse
!
!-------------------------------------------------------------------------------

!------------------ Parameter read and coefficients evaluation -----------------
!
subroutine read_and_allocate
  implicit none
  integer :: istat,i1,i2,i,il,node
  double precision :: aux,H,PG,PGc,tau,kappa,Rinv,defaultHPG
  character (len=10) :: str1
  namelist /Network/ N,Nplants,L,Nassets          !network parameters
  namelist /IntegrationParameters/ dt,t0,Ndspch,nstep,writeStep !integration parameters
  namelist /ModelParameters/ wR,D,tau_ou,epsilon,defaultHPG  !model parameters
  namelist /fileName/ initFile,dispFile,resFile

  open(10,file='parameters.dat',action='read')
!
  read (nml=Network,unit=10)
  Nt=2*N+2*Nplants ! total number of variables
  nzcoeffs=2*L+4*N+6*Nplants !=2*(L+Nt+Nplants), number of nonzero terms in I-dt*J
!
! Allocation of arrays for all nodes
  allocate(plantId(N)) ! power plant id (=<0 for consumers)
  allocate(Pl(N),dPl(N),Pleff(N),xi_ou(N),Gaussian(N))
  allocate(P_basket(N))
! Allocation of arrays for power plants
  allocate(Pref(Nplants),itau(Nplants),lambda(Nplants))
  allocate(wR22HPG(N),wR22HPGon(Nplants),wR22HPGoff(Nplants),wR22HPGstp(Nplants))
  allocate(PGRwR(Nplants),PGRwRon(Nplants),PGRwRoff(Nplants),PGRwRstp(Nplants))
  allocate(kwR(Nplants),kwRon(Nplants),kwRoff(Nplants),kwRstp(Nplants))
  allocate(on_off(Nplants),wR22HPGoff0(Nplants))
!
  allocate(plantNode(Nplants),plantStatus(Nplants),slowOpId(Nplants))
  allocate(assetNode(Nassets),assetKind(Nassets),P_assets(Nassets),P_assets0(Nassets),dP_assets0(Nassets))
! Allocation of arrays for lines
  allocate(Ml(L,2),B(L))
! Allocation of arrays for variables and flow
  allocate(x(Nt),deltax(Nt),deltaxs(Nt))
! Allocation of arrays associated to matrix I-dtJ
  allocate(IdtJ(nzcoeffs),idxcol(nzcoeffs),idxrow(Nt+1))
  allocate(idxLineA(L),idxLineB(L),idxFreqEqPhase(N),idxFreqEqFreq(N))

! Read namelists:
  read(nml=IntegrationParameters,unit=10)

  read(nml=ModelParameters,unit=10)
    if (wR.le.0.d0) then
      print *, 'wR not positive. wR =', wR; STOP
    endif
    wR=twoPi*wR ! wR is in Hz in the parameters file
    DwR=D/wR
    wR22HPG=wR*wR/(2.d0*defaultHPG) !default inertia coeff for consumer nodes.

    read (nml=fileName,unit=10)
!
    plantId=-1  ! Initialize plantId without plants.
!
! Power plants parameters
23  read(10,*,ERR=23) i1  !ignore coments until a line starting with an integer is found
    do i=1,Nplants
      read(10,*) k,node,H,PG,PGc,tau,Rinv,kappa
      if ((node.lt.1).or.(node.gt.N)) then
        print *, 'Fail to locate plant ', k, ' on node ', node, '. Out of range.', node; STOP
      else if (plantId(node).ne.-1) then
        print *, 'Fail to locate plant ', k, ' on node ', node, '. Already occupied by plant', plantId(node); STOP
      else if (tau.le.0.d0) then
        print*, 'tau not positive for plant ', k, tau; STOP
      else if ((H.le.0).or.(PG.le.0)) then
        print *, 'H or PG not positive for plant ', k, 'on node ', node, '. H, PG=', H, PG; STOP
      endif
      plantNode(k)=node     ! place plant k on node
      plantId(node)=k   ! ocupy node with plant k

      itau(k)=1.0d0/tau
      wR22HPGon(k)=wR*wR/(2.d0*H*PG)
      wR22HPGoff0(k)=wR*wR/(2.d0*defaultHPG)
      PGRwRon(k)=PGc*Rinv/wR
      kwRon(k)=kappa/wR
    enddo !i

! Line parameters
25  read(10,*,ERR=25) i1  !ignore coments until a line starting with an integer is found
    do i=1,L
      read(10,*) il,Ml(il,1),Ml(il,2),aux
      B(il)=aux
    enddo !i
    read(10,*) str1
    do i=1,Nassets
      read(10,*) j,assetNode(j),assetKind(j)
    enddo

! Pardiso parameters
    do
      read (10,*,IOSTAT=istat) i1,i2
      if (istat.eq.0) iparm(i1)=i2
      if (istat.lt.0) exit
    enddo
  close(10)
  return
end subroutine read_and_allocate
!-------------------------------------------------------------------------------

!------------------ Initialization ---------------------------------------------
!
subroutine initialization
  implicit none
  double precision :: Pe1,Pe2,time
  character (len=10) :: str1

  t=t0          ! initial time
  write(*,*) 'Implicit integration method using sparse matrices'
  write(*,*) 'dt=',real(dt),'seconds' ! print integration step

  call read_dispatch(Pl,Pref,P_assets0) ! Read data for initial dispatch

  if (initFile.eq.'autoinit') then
     x(1:2*N)=0.d0
     x(2*N+1:2*N+Nplants)=Pref
     x(2*N+Nplants+1:2*N+2*Nplants)=Pref
  else
    open(10,file=trim(initFile),action='read') ! Read initial condition
    do i=1,Nt
      read(10,*) x(i)
    enddo
    close(10)
  endif
!
! Initialize power plants
  do k=1,Nplants
    if (Pref(k).gt.0.d0) then
      call plantOp('InitOn',k)   ! set plant k as initally on
    else
      call plantOp('InitOff',k)  ! set plant k as initially offline
    endif
  enddo
!
  Pleff=Pl
  P_assets=P_assets0
!
! Initialize Orstein-Uhlenbeck noise
  p_ou=exp(-dt/tau_ou)
  coeff_ou=sqrt((1.0d0-p_ou**2)*0.5d0/tau_ou)
  call random_init(.false.,.false.)
  call rand_gaussian
  xi_ou=Gaussian/sqrt(2.d0*tau_ou)

 return
end subroutine initialization
!-------------------------------------------------------------------------------
!
!--------------------------- Gaussian random number ----------------------------
subroutine rand_gaussian
  implicit none
  integer :: N2
  double precision, dimension(2*int((N+1)/2)) :: U

  call random_number(U) ! array of uniform random numbers of mean 0 and variance 1

  N2=int((N+1)/2)

  ! Box-Muller algorithm to generate Gaussian random numbers:
  U(:N2)=sqrt(-2.d0*log(1.d0-U(:N2)))
  Gaussian(:N2)=U(:N2)*cos(twoPi*U(N2+1:2*N2))
  Gaussian(N2+1:N)=U(:N-N2)*sin(twoPi*U(N2+1:N))

  return
end subroutine rand_gaussian
!-------------------------------------------------------------------------------
!
!
!------------------ Dynamical equations for flow and Jacobian ------------------
!
subroutine dynamics
! Evaluates deltax and variable-dependent terms of I-dt*J (On rows N+1 to 2N associated to freq. eq.: diagonal, line terms and for plants Pm term)
!
! Input/Output: IdtJ. Nonzero coefficients of I-dt*J in format CSR3.
!
! Output: deltax
!
! Dynamics for all nodes:
!      d theta_i = w_i dt
!      d w_i     = ((wR^2)/2*H_i*PG_i(w_i+wR)){Pm_k - (1+D*w_i/wR)Pl_i - sum[ Bij*sin(theta_i-theta_j)]} dt
!   Consumers: Pm=0. Cross nodes: Pm=0, Pl=0.
!   Power plants, k=plantId(i), i=plantNode(k) and
!      d Pm_k = (1/tau_k)[ Ps_k - Pm_k - (PG/R)(w_i/wR) ] dt
!      d Ps_k = [-k(w_k/wR) - lambda(Ps_k-Pref_i)] dt
!
implicit none
integer :: i,il,ni,nj,k
double precision :: aux
!
  P_basket=0.d0
  do k=1,Nassets
    P_basket(assetNode(k))=P_basket(assetNode(k))+P_assets(k)
  enddo 
  IdtJ(idxFreqEqPhase)=0.d0   ! initializes sum for partial derivative freq-eq_i respect to phase_i
  do il=1,L
    ni=Ml(il,1)
    nj=Ml(il,2)
    aux=B(il)*sin(x(ni)-x(nj))     ! power transmitted by line il from ni to nj
    P_basket(ni)=P_basket(ni)-aux
    P_basket(nj)=P_basket(nj)+aux
    aux=B(il)*cos(x(ni)-x(nj))
!
    IdtJ(idxLineA(il))=-dt*wR22HPG(ni)*aux/(x(N+ni)+wR)  ! partial derivative freq-eq_ni respect to phase_nj
    IdtJ(idxFreqEqPhase(ni))=IdtJ(idxFreqEqPhase(ni))-IdtJ(idxLineA(il)) !update sum for partial derivative freq-eq_ni respect to phase_ni
    IdtJ(idxLineB(il))=-dt*wR22HPG(nj)*aux/(x(N+nj)+wR)  ! partial derivative freq-eq_nj respect to phase_ni
    IdtJ(idxFreqEqPhase(nj))=IdtJ(idxFreqEqPhase(nj))-IdtJ(idxLineB(il)) !update sum for partial derivative freq-eq_nj respect to phase_nj
  enddo
!
! For conventional power plants evaluate d Pm and d Ps for flow and partial derivative freq-eq_ni respect to Pm_i for Jacobian.
  do k=1,Nplants
    ni=plantNode(k)
    P_basket(ni)=P_basket(ni)+x(2*N+k)        ! Add Pm to P_basket
    deltax(2*N+k)=dt*itau(k)*(x(2*N+Nplants+k)-x(2*N+k)-PGRwR(k)*x(N+ni))
    deltax(2*N+Nplants+k)=-dt*kwR(k)*x(N+ni)-dt*lambda(k)*(x(2*N+Nplants+k)-Pref(k))
    IdtJ(idxrow(ni+N+1)-1)=-dt*wR22HPG(ni)/(x(N+ni)+wR)
  enddo !i
!
! For all nodes evaluate d theta and d w for flow and diagonal term of I-dt*J including partial derivative freq-eq_i respect to freq_i
  do i=1,N
    deltax(i)=dt*x(N+i)
    deltax(N+i)=dt*wR22HPG(i)*(P_basket(i)-(1.0d0+DwR*x(N+i))*Pleff(i))/(x(N+i)+wR)
    IdtJ(idxFreqEqFreq(i))=1.d0+dt*wR22HPG(i)*(P_basket(i)-(1.d0-D)*Pleff(i))/(x(N+i)+wR)**2
  enddo !i
!
  return
end subroutine dynamics
!
!-------------------------------------------------------------------------------

!------------------ Conventional Power Plant Operation -------------------------
!
subroutine plantOp(operation,iplant)
  implicit none
  integer :: iplant,ni,k
  character (LEN=*) :: operation
  ni=plantNode(iplant)
!
  if ((operation.eq.'SwitchOn').or.(operation.eq.'InitOn')) then
    wR22HPG(ni)=wR22HPGon(iplant)      ! Nominal inverse inertia coefficient
    PGRwR(iplant)=PGRwRon(iplant)      ! Nominal primary control
    kwR(iplant)=kwRon(iplant)          ! Nominal secondary control
    plantStatus(iplant)='On'
!
  else if (operation.eq.'SlowSwitchOff') then
    wR22HPG(ni)=wR22HPG(ni)+wR22HPGstp(iplant)   ! decrease inertia -> increase the coefficient in dw/dt.
    PGRwR(iplant)=PGRwR(iplant)+PGRwRstp(iplant) ! decrease primary control
    kwR(iplant)=kwR(iplant)+kwRstp(iplant)       ! decrease secondary control
    if ((wR22HPG(ni).ge.wR22HPGoff(iplant)).and.(PGRwR(iplant).le.PGRwRoff(iplant)).and.(kwR(iplant).le.kwRoff(iplant))) then
      wR22HPG(ni)=wR22HPGoff(iplant)
      PGRwR(iplant)=PGRwRoff(iplant)
      kwR(iplant)=kwRoff(iplant)
      plantStatus(iplant)='Off'
      do k=1,slowOps-1
        slowOpId(k)=slowOpId(k+1)
      enddo
      slowOps=slowOps-1
    endif
!
  else if (operation.eq.'InitOff') then
    x(2*N+iplant)=0.d0                ! Pm=0
    x(2*N+Nplants+iplant)=0.d0        ! Ps=0
    wR22HPG(ni)=wR22HPGoff(iplant)    ! inverse intertia set to off value.
    PGRwR(iplant)=PGRwRoff(iplant)    ! Primary control set to off value.
    kwR(iplant)=kwRoff(iplant)        ! Secondary contol set to off value.
    plantStatus(iplant)='Off'
!    `
  else
    print *, 'Power plant operation: ', operation, ' not defined'
    STOP
  endif
!
! update partial Pm_eq_iplant and Ps-eq_iplant respect to freq_plantNode(iplant)
  IdtJ(idxrow(2*N+iplant))=dt*itau(iplant)*PGRwR(iplant)
  IdtJ(idxrow(2*N+Nplants+iplant))=dt*kwR(iplant)
  return
end subroutine plantOp
!
!-------------------------------------------------------------------------------

!------------------ Dispatch ---------------------------------------------------
!
subroutine read_dispatch(Pl,Pref,P_assets0)
  implicit none
  double precision :: Pe1,lmbd,aux
  double precision, dimension (:) :: Pl,Pref,P_assets0
  integer :: i,k
  character (len=10) :: str1
!
  open(11,file=dispFile,action='read')
  Pl=0.d0
  Pref=0.d0
  P_assets=0.d0
  read (11,*) str1,pe1
  if (str1.ne."DISPATCH") then
    print *, 'Error in dispatch file. Found', str1, ' instead of DISPATCH' ; STOP
  endif
  read (11,*) str1
  do i=1,N+1
    read (11,*,err=31) j,Pe1
    Pl(j)=Pe1
  enddo
31 do k=1,Nplants+1
    read (11,*,err=32) j,Pe1,aux
    Pref(j)=Pe1
    on_off(j)=aux
  enddo
32 do k=1,Nplants+1
      read (11,*,err=33) j,lmbd
      if (lmbd.eq.0) then
        lambda(j)=0.0d0
      else
        lambda(j)=1.0d0/(60.0d0*lmbd) ! lmbd is in minutes in the dispatch file
      endif
    enddo
33 do i=1,Nassets+1
    read (11,*,err=34) j,Pe1
    P_assets0(j)=Pe1
   enddo
!
34 wR22HPGoff=wR22HPGon ! inertia on...
   do k=1,Nplants
     if (on_off(k).eq.0) wR22HPGoff(k)=wR22HPGoff0(k) !... or off (if the plant was offline in the reference case)
     !
     PGRwRoff(k)=on_off(k)*PGRwRon(k) ! primary control switched on or off(=0)
     kwRoff(k)=on_off(k)*kwRon(k) ! secondary control switched on or off(=0)
   enddo !k
!
return
end subroutine read_dispatch

subroutine dispatch
  double precision, dimension (:), allocatable :: Pl1,P_assets1
  integer :: k
  allocate (P_assets1(Nassets),Pl1(N))

  call read_dispatch (Pl1,Pref,P_assets1)
  dPl=(Pl1-Pl)/dfloat(nstep)        ! Load change per dt
  dP_assets0=(P_assets1-P_assets0)/dfloat(nstep)  ! Assets change per dt

  slowOps=0
  do k=1,Nplants
    if ((plantStatus(k).ne.'On').and.(Pref(k).gt.0.d0)) call plantOp('SwitchOn',k)
    if ((plantStatus(k).ne.'Off').and.(Pref(k).le.0.d0)) then
      slowOps=slowOps+1
      slowOpId(slowOps)=k
      wR22HPGstp(k)=(wR22HPGoff(k)-wR22HPGon(k))/dfloat(nstep)
      PGRwRstp(k)=(PGRwRoff(k)-PGRwRon(k))/dfloat(nstep)
      kwRstp(k)=(kwRoff(k)-kwRon(k))/dfloat(nstep)
      plantStatus(k)='Stopping'
    endif
  enddo
  return
end subroutine dispatch
!
!-------------------------------------------------------------------------------

!------------------ Write output data ------------------------------------------
subroutine write_data
  implicit none

  open(20,file=resFile,action='write')
  write(20,'(F10.5)') x(N+1)/twoPi
  twrite=t

  return
end subroutine write_data
!-------------------------------------------------------------------------------
end program grid_implicit_sparse
