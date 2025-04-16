c
c     ------------------------------------------------------------------
c     Main routine for calculating the gradients of the CNEO-TDDFT
c     (under frozen nuclear orbital assumption)excitation energies 
c     with respect to nuclear coordinates.
c
c     Workflow modified from tddft_gradient.F, with most of the functions
c     replaced by their CNEO-TDDFT counterparts.
c     ------------------------------------------------------------------
c
      logical function cneo_tddft_gradient(rtdb)
c
c     Notation:  Xia (i=occupied,a=virtual) 
c     Occupied index runs faster than the virtual index.
c
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "xc.fh"
#include "rtdb.fh"
#include "cosmo.fh"
#include "tddft_grad_util.fh"
#include "nwc_const.fh"
#include "ccneo.fh"
c
      character*255 title ! the job title
c
      integer rtdb        ! the runtime database handle
      integer ipol        ! =1 (restricted), =2 (unrestricted)
      logical tda         ! True if Tamm-Dancoff approximation
      logical oskel       ! Use of symmetry in Fock build
      logical osinglet    ! Compute singlet excited states?
      logical otriplet    ! Compute triplet excited states?
      logical lhascd
      integer nsroots     ! the number of singlet roots
      integer ihdl_sroots ! the singlet roots handle
      integer iptr_sroots ! the singlet roots list
      integer ntroots     ! the number of triplet roots
      integer ihdl_troots ! the triplet roots handle
      integer iptr_troots ! the triplet roots list
      integer nroots      ! the number of roots
      integer ihdl_roots  ! the roots handle
      integer iptr_roots  ! the roots list
      integer ntotroots   ! the total number of roots
      integer ihdl_bfao   ! the AO basis set handle
      integer ihdl_geom   ! the geometry handle
c
      integer nao         ! the number of AOs
      integer nmo(2)      ! the number of MOs
      integer nfc(2)      ! the number of frozen core orbitals
      integer nfv(2)      ! the number of frozen virtual orbitals
      integer nocc(2)     ! the number of occupied orbitals
      integer naoc(2)     ! the number of active occupied orbitals
      integer nav(2)      ! the number of active virtual orbitals
      integer nov(2)      ! nov(i) = naocc(i)*nav(i) 
      integer nat         ! the number of atoms
c
      integer ihdl_mo_e   ! mo energy memory handle
      integer iptr_mo_e   ! mo energy memory index
      integer ihdl_root_e ! excitation energy memory handle
      integer iptr_root_e ! excitation energy memory index
c
      double precision tol2e ! 2-electron integral tolerance
c
      double precision solve_thresh ! convergence threshold for solver
c
      integer solve_maxit  ! the maximum number of iterations in solver
      integer solve_maxsub ! the maximum number expansion vectors
                           ! per root
c
      character*255 fn_solve   ! filename for linear system solutions
      character*6   solve_when ! when to store solutions
c
      character*255 fn_civecs         ! filename for unrestricted vectors
      character*255 fn_civecs_singlet ! filename for singlet state vectors
      character*255 fn_civecs_triplet ! filename for triplet state vectors
      character*255 fn_movecs         ! filename for MO vectors
c
      character*40 xname(numfunc)       ! Names of X functionals
      character*40 cname(numfunc)       ! Names of C functionals
      character*40 xcname(numfunc)      ! Names of XC functionals
      double precision xfac(numfunc)    ! Prefactors for X functionals
      double precision cfac(numfunc)    ! Prefactors for C functionals
      logical lxfac(numfunc)            ! True if local X functional
      logical nlxfac(numfunc)           ! True if gc X functional
      logical lcfac(numfunc)            ! True if local C functional
      logical nlcfac(numfunc)           ! True if gc C functional
      logical xccomb(numfunc)           ! True if composite XC functional
c
      logical lhf                       ! True if HF (no DFT)
      logical ldft                      ! True if pure DFT (no HF exch)
      logical lhybrid                   ! True if hybrid DFT and HF
      logical lhasdft                   ! True if it has DFT
      logical lhashf                    ! True if it has HF
c
c     Global array for various quantities 
c     g_something(1) = alpha component (also in RHF case)
c     g_something(2) = beta  component
c
      integer g_mo(2)      ! global arrays to the MOs
      integer g_xpy(2)     ! global arrays to hold X+Y
      integer g_xmy(2)     ! global arrays to hold X-Y
      integer g_tp(2)      ! global arrays to hold T, and later P
      integer g_r(2)       ! global arrays to hold right-hand-sides
      integer g_z(2)       ! global arrays to hold solutions Ze
      integer g_w(2)       ! global arrays to hold energy weighted density matrices
      integer g_g          ! global array to hold the gradients
c
      double precision kfac ! Prefactor for HF exchange
c
      double precision tddft_grad_wall ! wall clock time
      double precision tddft_grad_cpu  ! cpu time
c
      logical cneo_tddft_energy
      external cneo_tddft_energy
      logical oroot
      logical odebug
c
      character*32 pname
c
      parameter (oskel=.false.)
      logical task_gradient
      external task_gradient
      logical task_energy
      external task_energy
c     === Nuclear variables === 
      integer ihdl_bfn(max_quantum_nuc)   ! basis handles for quantum nuc
      integer nao_n(max_quantum_nuc)      ! number of AOs for quantum nuc
      integer nmo_n(max_quantum_nuc)      ! number of MOs for quantum nuc
      integer naoc_n(max_quantum_nuc), nav_n(max_quantum_nuc)
      integer nov_n(max_quantum_nuc)
      integer g_mo_n(max_quantum_nuc)     ! global arrays to the nuclear MOs
      integer ihdl_mo_e_n, iptr_mo_e_n    ! nuclear mo energy handle/pointer
      integer g_r_n(max_quantum_nuc)      ! global arrays to hold nuclear RHSs
      integer g_z_n(max_quantum_nuc)      ! global arrays to hold nuclear solutions Zn
      double precision tol2e_cneo         ! tolerance for 2-body integrals (en, nn)
c
c MN solvation models -->
c
      integer do_cosmo_vem, istep_cosmo_vem
      integer g_vem(3)
c
      if (.not. rtdb_get
     $   (rtdb,'cosmo:do_cosmo_vem',mt_int,1,do_cosmo_vem))
     $    do_cosmo_vem = 0
c
c <-- MN solvation models
c
      pname = "cneo_tddft_gradient: "
c
      oroot = ga_nodeid().eq.0
c
c     Debug flag
      odebug = .false.
      if (rtdb_get(rtdb,'tddft_grad:debug', mt_log, 1, odebug))
     &     odebug = .true.
c
c     TDDFT gradient has not been computed successfully yet.
      cneo_tddft_gradient = .false.
c
c     Make sure we have the excited states first
      call ga_sync
      if (.not.cneo_tddft_energy(rtdb)) return
      call ga_sync
c
c     Get the start times
      tddft_grad_wall = -util_wallsec()
      tddft_grad_cpu  = -util_cpusec()
c
c     Init pstat
      call dft_pstat_init(rtdb)
c      
      call util_print_push
      call util_print_rtdb_load(rtdb,'tddft')
c
c     Print header
      if (oroot) then
         write(luout,*)
         write(luout,*)
         call util_print_centered(luout,
     +        'NWChem CNEO-TDDFT Gradient Module',40, .true.)
         write(luout,*)
         write(luout,*)
      endif
c
c     Parameter initialization
      call cneo_tdgrad_init(rtdb,ipol,tda,nat,
     +  nao,nocc,nmo,nfc,nfv,naoc,nav,nov,fn_movecs,numfunc,
     +  xname,cname,xcname,xfac,cfac,lxfac,nlxfac,lcfac,
     +  nlcfac,xccomb,lhf,ldft,lhybrid,lhasdft,lhashf,
     +  nroots, iptr_roots, ihdl_roots, fn_civecs,
     +  nsroots,iptr_sroots,ihdl_sroots,fn_civecs_singlet,
     +  ntroots,iptr_troots,ihdl_troots,fn_civecs_triplet,
     +  osinglet,otriplet, ntotroots,ihdl_bfao,ihdl_geom,tol2e,
     +  solve_thresh,solve_maxit,solve_maxsub,
     +  fn_solve,solve_when,kfac,title,lhascd,
     +  ihdl_bfn, nao_n, nmo_n, naoc_n, nav_n, nov_n,
     +  tol2e_cneo)
c
c     Cosmo initialization, if necessary
      cosmo_on = .false.
      if (rtdb_get(rtdb,'slv:cosmo',mt_log,1,cosmo_on)) then
       if(cosmo_on) then
            call errquit(pname//'cosmo not implemented',0,0)
       endif
      endif
c
c     Integral initialization
      call scf_get_fock_param(rtdb,tol2e)
      call cneo_tdgrad_intfile
c
c     Print calculation information
      call tddft_grad_print_info(luout,nroots,int_mb(iptr_roots),
     +     nsroots,int_mb(iptr_sroots),ntroots,int_mb(iptr_troots),
     +     title)
c
c     Print start times
      if (oroot) then
        if (util_print('tddft grad time',print_default)) then
          write(luout,*)
          write(luout,"(' Start at time  cpu:',f11.1,'s  wall:',f11.1,
     +                  's')")-tddft_grad_cpu,-tddft_grad_wall
          write(luout,*)
        endif
      end if
c
c     Allocate orbital energies
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_create_e"
      call tddft_grad_create_e(ipol,nao,ihdl_mo_e,iptr_mo_e)
      call tddft_grad_create_e(nuc_num,nao_n(1),ihdl_mo_e_n,
     $                                                iptr_mo_e_n)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_create_e"
c
c     Allocate MO coefficients.
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_create_mo"
      call tddft_grad_create_mo(ipol,ihdl_geom,ihdl_bfao,g_mo)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_create_mo"
c
c     Allocate (X+Y) & (X-Y)
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_create_xy"
      call tddft_grad_create_xy(tda,ipol,naoc,nav,ntotroots,g_xpy,g_xmy)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_create_xy"
c
c     Allocate T, T+(X+Y), or T+V(=P)
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_create_tp"
      call tddft_grad_create_tp(ipol,naoc,nav,ntotroots,g_tp)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_create_tp"
c
c     Allocate excitation energies
      if (oroot.and.odebug) 
     &   write(luout,*)"call tddft_grad_create_root_e"
      call tddft_grad_create_root_e(ntotroots,ihdl_root_e,iptr_root_e)
      if (oroot.and.odebug) 
     &   write(luout,*)"done tddft_grad_create_root_e"
c
c     Load the MO vectors from the file specified by fn_movecs
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_load_mo"
      call tddft_grad_load_mo(ipol,nao,fn_movecs,dbl_mb(iptr_mo_e),g_mo)
      call cneo_tdgrad_mo_n(rtdb,nuc_num,ihdl_geom,ihdl_bfn,nao_n(1),
     $                              g_mo_n,dbl_mb(iptr_mo_e_n))
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_load_mo"
c
c     Load CI vectors from TDDFT calculation: (X+Y) and (X-Y)
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_load_xy"
      call tddft_grad_load_xy(tda,ipol,naoc,nav,nov,nroots,
     +     int_mb(iptr_roots),fn_civecs,nsroots,int_mb(iptr_sroots),
     +     fn_civecs_singlet,ntroots,int_mb(iptr_troots),
     +     fn_civecs_triplet,g_xpy,g_xmy,dbl_mb(iptr_root_e))
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_load_xy"
c
c     Compute the TDDFT energy for each root using the solution vectors
      if (oroot.and.odebug)
     &    write(luout,*)"call tddft_grad_compute_energy"
      call cneo_schwarz_set(ihdl_bfao)
      call tddft_grad_compute_energy(rtdb,ihdl_geom,ihdl_bfao,tol2e,tda,
     +   ipol,nao,nfc,naoc,nocc,nav,nfv,ntotroots,
     +   dbl_mb(iptr_mo_e),g_mo,g_xpy,g_xmy,kfac,lhashf,otriplet,lhascd)
      if (oroot.and.odebug)
     &    write(luout,*)"done tddft_grad_compute_energy"
c
c     Compute the unrelaxed difference density T
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_compute_t"
      call tddft_grad_compute_t(tda,ipol,naoc,nav,ntotroots,g_xpy,g_xmy,
     +     g_tp)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_compute_t"
c
c     Allocate right hand side (R: g_r) and solution arrays (Z: g_z)
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_create_r"
      call tddft_grad_create_r(ipol,naoc,nav,ntotroots,g_r)
      call tddft_grad_create_r(ipol,naoc,nav,ntotroots,g_z)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_create_r"
c
c     Assemble right hand side: R
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_compute_r"
      call tddft_grad_compute_r(rtdb,ihdl_geom,ihdl_bfao,tol2e,tda,ipol,
     +     nao,nfc,naoc,nocc,nav,nfv,ntotroots,g_mo,g_xpy,g_xmy,g_tp,
     +     g_r,kfac,lhashf,otriplet,lhascd)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_compute_r"

      call cneo_tdgrad_r_n(nuc_num,ihdl_geom,ihdl_bfao,ihdl_bfn,
     $      nao,naoc,nav,nocc,nfc,nfv,nao_n,naoc_n,nav_n,g_mo,g_mo_n,
     $      g_tp,ntotroots,ipol,g_r_n,tol2e_cneo)
      call cneo_tdgrad_create_r_n(nuc_num,naoc_n,nav_n,ntotroots,g_z_n)
c
c     Solve the CNEO coupled (Z_e,Z_n) equation
c     g_z, g_z_n have the solutions
      if (oroot.and.odebug) write(luout,*) "call cneo_tdgrad_solve"

      call cneo_tdgrad_solve(rtdb,ihdl_geom,ihdl_bfao,ipol,
     +      ntotroots,nao,nocc,naoc,nav,nfc,nfv,nmo,g_mo,iptr_mo_e,
     +      .false.,.false.,lhashf,otriplet,lhascd,tol2e,kfac,
     +      solve_thresh,solve_maxsub,solve_maxit,
     +      ihdl_bfn,nao_n,naoc_n,nav_n,nmo_n,g_mo_n,iptr_mo_e_n,
     +      g_r,g_r_n,g_z,g_z_n,tol2e_cneo)

      if (oroot.and.odebug) write(luout,*) "done cneo_tdgrad_solve"
c
c     Deallocate right hand side: R
      if (oroot.and.odebug) write(luout,*) "call cneo_tdgrad_destroy_r"
      call tddft_grad_destroy_r(ipol,g_r)
      call cneo_tdgrad_destroy_r_n(nuc_num,g_r_n)
      if (oroot.and.odebug) write(luout,*) "done cneo_tdgrad_destroy_r"
c
c     Compute the relaxed density matrix: P=T+Z
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_compute_p"
      call tddft_grad_compute_p(ipol,naoc,nav,ntotroots,g_tp,g_z)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_compute_p"
c
c     Create Lagrange multipliers: W
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_create_tp"
      call tddft_grad_create_tp(ipol,naoc,nav,ntotroots,g_w)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_create_tp"
c
c     Compute Lagrange multipliers: W
c     W is the energy-weighted difference density matrix
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_compute_w"
      call cneo_tdgrad_compute_w(rtdb,ihdl_geom,ihdl_bfao,tol2e,tda,
     +     ipol,ntotroots,nfc,naoc,nocc,nav,nfv,nao,g_mo,g_tp,g_z,
     +     g_xpy,g_xmy,dbl_mb(iptr_mo_e),dbl_mb(iptr_root_e),g_w,
     +     kfac,lhashf,otriplet,lhascd,nuc_num,
     +     ihdl_bfn, naoc_n,nav_n,nao_n,g_mo_n,g_z_n,tol2e_cneo)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_compute_w"
c
c     Finalize integrals
      if (oroot.and.odebug)write(luout,*)"call tddft_grad_finalize_ints"
      call cneo_tdgrad_finalize_ints(rtdb,nuc_num,ihdl_bfao,ihdl_bfn)
      if (oroot.and.odebug)write(luout,*)"done tddft_grad_finalize_ints"
c
c     Create gradients array
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_create_g"
      call tddft_grad_create_g(nat,ntotroots,g_g)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_create_g"
c
c     Compute the gradients
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_compute_g"
      call cneo_tdgrad_compute_g(
     +  rtdb,ihdl_geom,ihdl_bfao,tda,nat,nao,
     +  nfc,naoc,nocc,nav,nfv,ipol,ntotroots,oskel,tol2e,g_mo,
     +  g_tp,g_w,g_xpy,g_xmy,g_g,kfac,lhashf,otriplet,
     +  nuc_num,nuc_cent,ihdl_bfn,nao_n,naoc_n,nav_n,g_mo_n,g_z_n)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_compute_g"
c
c     Deallocate arrays
      if (oroot.and.odebug) write(luout,*) "call deallocate arrays"
      call tddft_grad_destroy_tp(ipol,g_w)
      call tddft_grad_destroy_r(ipol,g_z)
      call cneo_tdgrad_destroy_r_n(nuc_num,g_z_n)
      call tddft_grad_destroy_tp(ipol,g_tp)
      call tddft_grad_destroy_xy(tda,ipol,g_xpy,g_xmy)
      call tddft_grad_destroy_mo(ipol,g_mo)
      call cneo_tdgrad_destory_mo_n(nuc_num,g_mo_n)
      call tddft_grad_destroy_root_e(ihdl_root_e)
      call tddft_grad_destroy_e(ihdl_mo_e_n)
      call tddft_grad_destroy_e(ihdl_mo_e)
      if (oroot.and.odebug) write(luout,*) "done deallocate arrays"
c MN solvation models -->
c
c We don't need to print out the excited-state gradients in the case of VEM
c
      if (do_cosmo_vem.ne.0) goto 100
c
c <-- MN solvation models
c
c     Print gradients of the excited states of interest
      if (oroot.and.odebug) 
     &     write(luout,*) "call tddft_grad_print_result"
      if (ipol.eq.2) then
        call cneo_tdgrad_print_result(ntotroots,int_mb(iptr_roots),
     +                             nat,ihdl_geom,g_g)
      else
        if (osinglet) then
          call cneo_tdgrad_print_result(ntotroots,int_mb(iptr_sroots),
     +                               nat,ihdl_geom,g_g)
        else if (otriplet) then
          call errquit(pname//'triplet not implemented',0,0)
        endif
      endif
      if (oroot.and.odebug) 
     &      write(luout,*) "done tddft_grad_print_result"
c
 100  continue
c
c     Deallocate gradients
      if (oroot.and.odebug)
     &      write(luout,*) "call tddft_grad_destroy_g"
      call tddft_grad_destroy_g(g_g)
      if (oroot.and.odebug)
     &      write(luout,*) "done tddft_grad_destroy_g"
c
c     Finalize gradients
      if (oroot.and.odebug) 
     &      write(luout,*) "call tddft_grad_finalize"
      call cneo_tdgrad_finalize(rtdb, nsroots,ihdl_sroots,ntroots,
     +            ihdl_troots,nroots,ihdl_roots,
     +            ihdl_geom,ihdl_bfao,nuc_num,ihdl_bfn)
      if (oroot.and.odebug) 
     &      write(luout,*) "done tddft_grad_finalize"
c
c     Clean up cosmo stuff, if necessary
      if (oroot.and.odebug.and.cosmo_on)
     &      write(luout,*) "call cosmo_tidy"
      if (cosmo_on) then
         call cosmo_tidy(rtdb)
         cosmo_on = .false.
      endif
      if (oroot.and.odebug.and.cosmo_on)
     &      write(luout,*) "done cosmo_tidy"
c
c     Clean up
      call ga_sync()
      call grid_cleanup(.true.)
      cneo_tddft_gradient = .true.
c
c     Print end times and wrap up
      tddft_grad_cpu  = tddft_grad_cpu  + util_cpusec()
      tddft_grad_wall = tddft_grad_wall + util_wallsec()
      if (oroot) then
        if (util_print('tddft grad time',print_default)) then
          write(LuOut,*)
          write(LuOut,"(' CNEO-TDDFT Gradient time  cpu:',f11.1,'s  
     +           wall:',f11.1,'s')")tddft_grad_cpu,tddft_grad_wall
        endif
      endif
c
c MN solvation models -->
c
c cosmo-vem: Do one more TDDFT iteration. The whole SCF--TDDFT energy--TDDFT density
c calculation will be repeated as many times as needed to achieve 
c convergence on the VEM energy as determined by hnd_coschg in src/property/cosmo.F.
c When such a convergence is achieved istep_cosmo_vem is set to 3 and the program exits the loop.  
c
      cosmo_on = .false.
      if (rtdb_get(rtdb,'slv:cosmo',mt_log,1,cosmo_on)) then
       if(cosmo_on) then
        istep_cosmo_vem = 0
        if (.not. rtdb_get
     $(rtdb,'cosmo:istep_cosmo_vem',mt_int,1,istep_cosmo_vem))
     $  call errquit('tddft_grad: cannot get istep_cosmo_vem',
     $ 0,rtdb_err)
        if (istep_cosmo_vem.eq.2.and.do_cosmo_vem.ne.0) then
         if (.not.task_gradient(rtdb)) call errquit
     $('error: task_gradient(rtdb) failed for cosmo-vem',
     & 0,geom_err)
        endif
c
c cosmo-vem: Do an additional DFT energy calculation if emission is requested overwise do nothing. 
c
        if (istep_cosmo_vem.eq.3.and.do_cosmo_vem.eq.2) then
         if (.not. rtdb_cput(rtdb, 'task:theory', 1, 'dft'))
     $ call errquit('error: failed changing theory for for cosmo-vem',
     1 0,rtdb_err)
         if (.not. rtdb_put(rtdb,'dft:converged', mt_log, 1, .false.))
     $ call errquit("task: could not write to rtdb",
     $ 0,rtdb_err)
         if (.not.task_energy(rtdb)) call errquit
     $('error: task_energy(rtdb) failed for cosmo-vem',
     & 0,geom_err)
        endif
c
c cosmo-vem: Destroy global arrays
c
        if (istep_cosmo_vem.eq.3) then
         if(rtdb_get(rtdb,'cosmo:g_vem',mt_int,3,g_vem)) then
          if (.not. ga_destroy(g_vem(1))) call errquit
     $        ('hnd_coschg: could not destroy g_vem(1)', 0, ga_err)
          if (.not. ga_destroy(g_vem(2))) call errquit
     $        ('hnd_coschg: could not destroy g_vem(2)', 0, ga_err)
          if (.not. ga_destroy(g_vem(3))) call errquit
     $        ('hnd_coschg: could not destroy g_vem(3)', 0, ga_err)
          call ga_sync()
         endif
        endif
       endif
      endif
c
c <-- MN solvation models 
c
      call dft_pstat_print
      call util_print_pop
      return
      end
c $Id$
