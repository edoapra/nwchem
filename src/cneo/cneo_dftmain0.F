c
c     == Main CNEO-DFT driver ==
      logical function cneo_dft_main0d(rtdb)
c     
c     $Id$
c     
      implicit none
c
#include "errquit.fh"
c
      integer rtdb              ! rtdb handle [input]
      Logical LResult, IOLGC, oprint, oprint_time, oprint_ioinfo
c     
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "dft_conv.fh"
#include "global.fh"
#include "util.fh"
#include "tcgmsg.fh"
#include "schwarz.fh"
#include "msgids.fh"
#include "dftps.fh"
#include "bgj.fh"
#include "ccneo.fh"
#include "ccneo_dft.fh"
c
#ifdef SCALAPACK
#include "dft_scaleig.fh"
      integer ga_cholesky,ga_llt_i
      external ga_cholesky,ga_llt_i
#endif
c     
#include "zora.fh"
#include "rel_consts.fh"
#include "x2c.fh"
#include "case.fh"
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      integer g_dens(2), g_evecs(2), g_vxc(4)
c
      integer g_dens_at(2)
c
      integer g_svecs, g_xcinv, g_oi, g_svecs_n(max_quantum_nuc)
      integer g_s, g_fock, g_tmpb, g_tmp, g_s_n(max_quantum_nuc)
      integer nvxc
      integer me, nproc, info
      double precision Enuc, Etold
      integer iter, nconstr
      integer iVxc_opt, lmiss, ivcoul_opt
      integer isvals, lsvals,  idum,n_dep
      integer isvals_n,lsvals_n
      double precision toll_s
c
      logical lsuccess, ltried_scf, ltried_cgmin, lcgmin
      logical cneo_dft_scf,atmblk, lcdft, dft_cscf_equals_scftype
      external cneo_dft_scf, dft_cscf_equals_scftype
      logical dft_scf_so
      external dft_scf_so
      logical  dft_uks_nr_solve, dft_roks_nr_solve
      external dft_uks_nr_solve, dft_roks_nr_solve
      logical  dft_cg_solve
      external dft_cg_solve
c
c     == RT_TDDFT ==
      logical lrt_tddft       ! real-time TDDFT dynamics flag
      integer g_dipx  ! ga for dipole x
      integer g_dipy  ! ga for dipole y
      integer g_dipz  ! ga for dipole z
      integer g_vel   ! velocity ga's
      integer g_ddotE ! ga for D.E
      integer alo(3),ahi(3)
      double precision origin(3)
      data origin/0d0,0d0,0d0/
c
      character*80 theory
      logical status
      logical ldirect, xc_active
      character*32 pname
      character*1 yev0
c 
      double precision eone, etwo, enrep, energy
      double precision dft_time
c
      logical largenode
c
      logical lxps
      integer i,sval_shift
c
      nexc = 1
c
c     == Preliminaries ==
      cneo_dft_main0d = .false.
      pname = "cneo_dft_main0d:"
      me = ga_nodeid()
      nproc = ga_nnodes()
      oprint = util_print('cneo dft info', print_default)
      oprint_time = util_print('cneo dft timings', print_high)
      oprint_ioinfo = util_print('io info', print_debug)
      if(ga_nodeid().eq.0.and.
     .     util_print('geometry', print_high)) then
        if(.not.geom_print(geom)) call errquit(
     .       pname//'geom print failed',0, GEOM_ERR)
      endif
      IOLGC = .TRUE.
      if (noio.eq.1) IOLGC = .FALSE.
c
      lcgmin = .false.
      lsuccess = .false.
      ltried_scf = .false.
      ltried_cgmin = .false.
c
c     initialize
      if (.not. rtdb_put(rtdb, 'cneo:energy', MT_DBL, 1, 0.d0))
     &   call errquit('dft_scf: failed to store energy in rtdb', 0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'cneo:converged', MT_LOG,1,.false.))
     &   call errquit('dft_scf: failed to store converged in rtdb', 0,
     &       RTDB_ERR)
c
c     == SCALAPACK related ==
c
#ifdef SCALAPACK
#if defined(PARALLEL_DIAG)
      call slexit
      call slexit2
      if (.not. rtdb_cget(rtdb, 'dft:scaladiag',1,yev))
#ifdef USE_PDSYEVR__
     &     yev='r'
#else
     &     yev='d'
#endif
      if (.not. rtdb_get(rtdb, 'dft:scalanb',mt_int,1,nb))
     &     nb=0
#endif
#endif
c
c     == Trash junk that may have remained in restart calcs ==
      if(rtdb_get(rtdb,'bgj:xc_active', MT_LOG, 1,xc_active)) then
        if (xc_active) then
c         Junk left from previous run, clobber it!!!
c         Leaving it causes trouble in fock_2e because it then also
c         calls fock_xc. The routine dft_fockbld does not expect to be
c         using fock_xc anywhere but computes the XC-contributions by
c         a different route. This conflict will cause the calculation
c         to fail.
          if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     &      call errquit(pname//"could not store xc_active",
     &                   0,RTDB_ERR)
        endif
      endif
c
c     == Check which solver to use ==
      if (.not.rtdb_get(rtdb, 'dft:cgmin', mt_log, 1, lcgmin))
     &   lcgmin=.false.
c
c     == Check real-time TDDFT flag ==
      lrt_tddft = .false.
c
c     == Evaluate the nuclear repulsion energy ==
      LResult = geom_nuc_rep_energy(geom, Enuc)
      if (.not.LResult) 
     & call errquit(pname//'Error in nuclear repulsion',1, BASIS_ERR)
c     
c     == Determine the characteristics of the AO Gaussian basis set ==
c     == This uses the old NWdft data structures, rather than tthe rtdb ==
      if (.not. bas_numbf(AO_bas_han, nbf_ao)) then
         call errquit(pname//'Exiting',1, BASIS_ERR)
      endif 
c
c     Do XPS calculation
      lxps = .false.
c
      iter = 0
      Etold = 0.D0
c
c     == ZORA options ==
      do_zora_get_correction=.false.
      do_zora = .false.
c
c     turn off grid nodisk with zora
      if (.not.rtdb_get(rtdb,'dft:largenode', mt_log, 1, largenode))
     &     largenode=.false.
      if (do_zora.and.largenode) then
         call errquit('cneo_dft: unsupported function', 0, 0)
      end if ! do_zora and largenode
c
c     == Get the grid cutoff or set default ==
      if(.not.rtdb_get(rtdb,'zora:cutoff',mt_dbl,1,zoracutoff))
     &  zoracutoff = 1.0d-08
      if(.not.rtdb_get(rtdb,'zora:cutoff_EFG',          ! FA
     &                 mt_dbl,1,zoracutoff_EFG))        ! FA
     &  zoracutoff_EFG = 1.0d-08                        ! FA
      if(.not.rtdb_get(rtdb,'zora:cutoff_NMR',          ! FA
     &                 mt_dbl,1,zoracutoff_NMR))        ! FA
     &  zoracutoff_NMR = 1.0d-08                        ! FA
      if(.not.rtdb_get(rtdb,'zora:do_NonRel',           ! FA
     &                 mt_log,1,do_NonRel))             ! FA
     &  do_NonRel = .false.                             ! FA
      if(.not.rtdb_get(rtdb,'zora:not_zora_scale',      ! FA
     &                 mt_log,1,not_zora_scale))        ! FA
     &  not_zora_scale = .false.                        ! FA

c ... jochen:
c     == X2C options ==
      do_x2c = .false.
! c     we use the variable zoracutoff also for the x2c option
      if(.not.rtdb_get(rtdb,'x2c:cutoff',mt_dbl,1,zoracutoff))
     &  zoracutoff = 1.0d-08
c

c
c     == Get the light speed and other constants, if specified ==
      call rel_init_con(rtdb)
      clight_au = cau
c
c     == Coulomb Attenuation Method (CAM/LC) parameters ==
      if (.not.rtdb_get(rtdb, 'dft:cam_exch', mt_log, 1,
     &   cam_exch)) cam_exch=.false.
      if (.not.rtdb_get(rtdb, 'dft:cam_srhf', mt_log, 1,
     &   cam_srhf)) cam_srhf=.false.
      if (.not.rtdb_get(rtdb, 'dft:cam_omega', mt_dbl, 1,
     &   cam_omega)) cam_omega=0.d0
      if (.not.rtdb_get(rtdb, 'dft:cam_alpha', mt_dbl, 1,
     &   cam_alpha)) cam_alpha=0.d0
      if (.not.rtdb_get(rtdb, 'dft:cam_beta', mt_dbl, 1,
     &   cam_beta)) cam_beta=0.d0
c
      ldirect = .false.
c     
c     == Select method to compute the XC potential ==
      if (XC_bas_han.eq.-99999)then
         iVxc_opt = 0
         XCFIT    = .FALSE.
      else
         iVxc_opt = 1
         XCFIT    = .TRUE.
      endif

c
c     == Done with XC setup stuff ==
c     
c     == Select method to compute Coulomb potential ==
c     
      if (CD_bas_han.eq.-99999)then
         iVcoul_opt = 0
         CDFIT      = .FALSE.
      else
         iVcoul_opt = 1
         CDFIT      = .TRUE.
      endif
c
c
c     
c     Compute the AO overlap matrix.
c     
      if (.not. rtdb_get(rtdb, 'dft:atomblock',mt_log,1,atmblk))
     &     atmblk=.true.

      if(atmblk) then
         g_s = ga_create_atom_blocked(geom, AO_bas_han, 'AO ovl')
         do i=1, nuc_num
            g_s_n(i) = ga_create_atom_blocked(geom,bas_n(i),'nuc ovl')
            call ga_zero(g_s_n(i))
         enddo
      else
         if (.not. ga_create(mt_dbl, nbf_ao, nbf_ao, 'AO ovl',
     &        0, 0, g_s))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)
         do i=1,nuc_num
            if (.not. ga_create(mt_dbl, nbf_n, nbf_n, 'nuc ovl',
     &                                 0, 0, g_s_n(i)))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)
            call ga_zero(g_s_n(i))
         enddo
      endif
      call ga_zero(g_s)
      if (.not.MA_Push_Get(MT_Dbl, nbf_ao, 'ovl eig vals', lsvals, 
     &   isvals))
     &   call errquit(pname//'Cannot allocate ovl eig vals', 0,
     &       MA_ERR)
      call dfill(nbf_ao, 0.0d0, dbl_mb(isvals), 1)
      call int_1e_ga(AO_bas_han, AO_bas_han, g_s, 'overlap', oskel)
      if (oskel)call sym_symmetrize(geom, AO_bas_han, .false., g_s)

      if (.not.MA_Push_Get(MT_Dbl, nbf_n*nuc_num, 'ovl eig vals', 
     &      lsvals_n, isvals_n))
     &   call errquit(pname//'Cannot allocate ovl eig vals', 0,
     &       MA_ERR)
         call dfill(nbf_n*nuc_num, 0.0d0, dbl_mb(isvals_n), 1)

      do i=1,nuc_num
         call int_1e_ga(bas_n(i),bas_n(i),g_s_n(i),'overlap',.false.)
      enddo
c     
c     
c     orthog BS
c     
c     Diag Overlap Matrix
c     
      if (.not. ga_duplicate(g_s, g_svecs, 'AO ovl eig vecs'))
     &     call errquit(pname//'Error creating ga',0,GA_ERR)
      call ga_zero(g_svecs)
      if (odftps) call pstat_on(ps_diag)
#if defined(PARALLEL_DIAG)
#ifdef SCALAPACK
c     we need all eigenvectors for the overlap matrix,
         yev0=yev
         if(yev.eq.'m')  yev='r'
         idum=0
         call dft_scaleig(g_s, g_svecs, Dbl_MB(isvals), idum)
         yev=yev0
#else
      call ga_diag_std(g_s, g_svecs, Dbl_MB(isvals))
#endif
#else 
      call ga_diag_std_seq(g_s, g_svecs, Dbl_MB(isvals))
#endif
      if (odftps) call pstat_off(ps_diag)

      sval_shift = 0
      do i=1,nuc_num
         if (.not. ga_duplicate(g_s_n(i), g_svecs_n(i),
     $                                 'nuc ovl eig vecs'))
     &     call errquit(pname//'Error creating ga',0,GA_ERR)
         call ga_zero(g_svecs_n(i))
         if (odftps) call pstat_on(ps_diag)
#if defined(PARALLEL_DIAG)
#ifdef SCALAPACK
c     we need all eigenvectors for the overlap matrix,
            yev0=yev
            if(yev.eq.'m')  yev='r'
            idum=0
            call dft_scaleig(g_s_n(i), g_svecs_n(i), 
     $                        Dbl_MB(isvals_n+sval_shift), idum)
            yev=yev0
#else
         call ga_diag_std(g_s_n(i),g_svecs_n(i),
     $                           Dbl_MB(isvals_n+sval_shift))
#endif
#else 
         call ga_diag_std_seq(g_s_n(i), g_svecs_n(i), 
     $                           Dbl_MB(isvals_n+sval_shift))
#endif
         if (odftps) call pstat_off(ps_diag)
         sval_shift = sval_shift + nbf_n
      enddo
c     
      if (.not. ga_duplicate(g_s, g_evecs(1), 'alpha evecs'))
     &     call errquit(pname//'Error creating ga',0,GA_ERR)
      call ga_zero(g_evecs(1))
c
      if (.not. ga_duplicate(g_s, g_dens(1), 'alpha density'))
     & call errquit(pname//'Error creating ga',0,GA_ERR)
      call ga_zero(g_dens(1))
c
      if (.not. ga_duplicate(g_s, g_dens_at(1), 'alpha atomic density'))
     & call errquit(pname//'Error creating ga',0,GA_ERR)
      call ga_zero(g_dens_at(1))
c
      if (ipol.ne.1) then
       if (.not. ga_duplicate(g_s, g_dens(2), 'beta density'))
     &  call errquit(pname//'Error creating ga',0,GA_ERR)
       call ga_zero(g_dens(2))
c
       if (.not. ga_duplicate(g_s, g_dens_at(2), 'beta atomic density'))
     &  call errquit(pname//'Error creating ga',0,GA_ERR)
       call ga_zero(g_dens_at(2))
c
       if (.not. ga_duplicate(g_s, g_evecs(2), 'beta evecs'))
     &  call errquit(pname//'Error creating ga',0,GA_ERR)
       call ga_zero(g_evecs(2))
c
       if (.not. ga_duplicate(g_s, g_vxc(2), 'Vxcs'))
     &  call errquit(pname//'Error creating ga',0,GA_ERR)
       call ga_zero(g_Vxc(2))
c
      endif 
c
      call dft_cscf_movecs(g_evecs)
c
      nVxc = ipol
      if (.not. ga_duplicate(g_s, g_fock, 'fock matrix'))
     &     call errquit(pname//'Error creating ga',0,GA_ERR)
      call ga_zero(g_fock)
      g_vxc(1) = g_fock
c
      lkeeps=.false.
c
c     == Set orbital overlap tolerance ==
      call util_tolls(rtdb,.false.,toll_s,n_dep,dbl_mb(isvals),nbf_ao)
cc     
c     == Begin iteration sequence with S**(-1/2) in g_evecs arrays ==
      if (.not. ga_duplicate(g_s, g_tmp, 'temp'))
     &     call errquit(pname//'Error creating ga',0,GA_ERR)
      call diis_bld12(toll_s, dbl_mb(isvals), g_svecs, g_evecs(1), 
     &     g_tmp,  2)
      if (.not. ga_destroy(g_tmp)) call errquit
     &   (pname//'Could not destroy g_tmp', 0, GA_ERR)
      if (ipol.gt.1) call ga_copy(g_evecs(1), g_evecs(2))
c     
      if (.not. rtdb_get(rtdb, 'dft:lcdft', mt_log, 1, lcdft))
     &     lcdft = .false.
      nconstr = 0
c
      if (oprint_time) call dft_tstamp(' Before call to CNEO_DFT_SCF.')
c 
      if (.not. rtdb_cget(rtdb, 'dft:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)
      if (odftps) call pstat_on(ps_scfdft)
c
c
c     == Call the main drivers ==


      if (theory .eq. 'dft') then
       if (.not.lrt_tddft) then
          if ((.not.lcgmin)) then
            cneo_dft_main0d = cneo_dft_scf(rtdb, Etold, Enuc,iVcoul_opt,
     &                         iVxc_opt, iter, g_dens, g_dens_at,
     &                         g_evecs, g_vxc, g_fock, g_svecs,
     &                         dbl_mb(isvals),g_xcinv,g_s,g_s_n,lcdft,
     &                         nconstr, g_svecs_n, dbl_mb(isvals_n))
          endif
       else
          cneo_dft_main0d = .false.
       end if   ! lrt_tddft
      else if (theory .eq. 'sodft') then
         cneo_dft_main0d = .false.
      endif ! theory
c
      if (odftps) call pstat_off(ps_scfdft)
c     
c     == Clean up memory usage ==
      if (.not. ga_destroy(g_s)) call errquit
     &      (pname//'Could not destroy g_s', 0, GA_ERR)
      if (me.eq.0) then
         if (CDFIT.and.IOLGC) then
            call dft_invio('CD', idum, 1, 'KILL', lmiss)
            call dft_invio('CDI', idum, 1, 'KILL', lmiss)
         endif
         if ((.not.ADFT).and.XCFIT.and.IOLGC) then
            call dft_invio('XC', idum, 1, 'KILL', lmiss)
         endif
      endif
      if (.not. ga_destroy(g_fock)) call errquit
     &   (pname//'Could not destroy g_fock', 0, GA_ERR)
      if (ipol.eq.2)then
         if (.not. ga_destroy(g_vxc(2))) call errquit
     &      (pname//'Could not destroy g_vxc(2)', 0, GA_ERR)
      endif
      if (ipol.ne.1) then
         if (.not. ga_destroy(g_evecs(2))) call errquit
     &      (pname//'Could not destroy g_evecs(2)', 0, GA_ERR)
         if (.not. ga_destroy(g_dens(2))) call errquit
     &      (pname//'Could not destroy g_dens(2)', 0, GA_ERR)
         if (.not. ga_destroy(g_dens_at(2))) call errquit
     &      (pname//'Could not destroy g_dens_at(2)', 0, GA_ERR)
      endif   ! ipol
c
      if (.not. ga_destroy(g_evecs(1))) call errquit
     &   (pname//'Could not destroy g_evecs(1)', 0, GA_ERR)
      if (.not. ga_destroy(g_dens(1))) call errquit
     &   (pname//'Could not destroy g_dens(1)', 0, GA_ERR)
      if (.not. ga_destroy(g_dens_at(1))) call errquit
     &   (pname//'Could not destroy g_dens_at(1)', 0, GA_ERR)
      if (.not. ga_destroy(g_svecs)) call errquit
     &   (pname//'Could not destroy g_svecs', 0, GA_ERR)

      do i=1,nuc_num
         if (.not.ga_destroy(g_s_n(i))) call errquit
     +      (pname//'could not destroy g_s_n',0,GA_ERR)
         if (.not.ga_destroy(g_svecs_n(i))) call errquit
     +      (pname//'could not destroy g_svecs_n',0,GA_ERR)
      enddo

c
      if (.not.ma_chop_stack(lsvals))
     &  call errquit(pname//'Cannot chop stack lsvals', 12, MA_ERR)

c
      return
      end
