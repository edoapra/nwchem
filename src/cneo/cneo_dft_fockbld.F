c     ============================================
c     Build Fock matrices for CNEO-DFT calculations
c
c     Elec-nuc interaction block is added on the conventional DFT Fock
c     Nuclear Fock matrices are built (hcore_nuc + en + nn)
c     ============================================
c
      subroutine cneo_dft_fockbld(rtdb,natoms,ntotel,g_ndens,
     ,     g_fock,g_focks,g_dens,g_vxc,g_movecs,g_xcinv,
     ,     g_tmp,g_nq,g_wght,g_xyz,
     ,     ivcoul_Opt,nexc,ivxc_opt,wght_GA,
     ,     n_batch,n3c_int,n3c_dbl,IOLGC,
     ,     evals, rdens_atom,dd3ceri,ii3cwhat,
     ,     cetobfr, iter, ecore, ecoul,exc, rho_n,ldmix,
     ,     test_sic,exact_pot,ecoul_aux,dif_lamda,totsic,
     ,     i_degen,g_rho_exact,g_tmp_exact, g_addit_exact,
     ,     n_levels, act_levels, oprint_time,
cc AJL/Begin/FDE
c     ,     g_rel_sf,g_frozemb)
     ,     g_rel_sf,g_frozemb,g_vxc_fde,g_xcinv_fde,g_dens_fde,eelec,
     ,     g_nfock, ecore_n, e_en, e_nn)
cc AJL/End
c     
c     $Id$
c     
      implicit none
#include "errquit.fh"
#include "cdft.fh"
#include "global.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "dftps.fh"
#include "cosmo.fh"
#include "zora.fh"
#include "x2c.fh"
#include "case.fh"
#include "cfock.fh"
#include "ccneo.fh"
#include "ccneo_dft.fh"
#include "ccneo_constrain.fh"
c#include "dimqm.fh"
c
      integer g_fock,g_dens(2),g_vxc(4),g_movecs(2) ! [handle in/values out]
      integer g_focks(2),g_ndens(max_quantum_nuc)
      integer g_xcinv
      integer g_tmp,g_nq,g_wght,g_xyz
      integer g_rho_exact,g_tmp_exact, g_addit_exact
      integer rtdb ! [in]
      integer natoms ! [in]
      logical wght_GA
c     double precision evals (*)
      integer evals(2)
      double precision rdens_atom(*)
      double precision dd3ceri(*)
      integer ii3cwhat(*)
      integer cetobfr(*)
      logical oprint_time ![in]
      integer iter ! [in]
      integer n_levels(2), act_levels
      double precision ecore, ecoul,exc(2) ! [out]
      integer ntotel,n_batch,n3c_int,n3c_dbl
      logical IOLGC, ldmix
      integer test_sic ! [in]
      integer exact_pot, iter_ini
      integer i_degen(2)  ! [in]
      double precision ecoul_aux(2),dif_lamda,totsic
      double precision rho_n ![out]
      double precision fant_a,fant_d,dft_fant
      external dft_fant
      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      integer ivcoul_Opt,nexc,ivxc_opt
c
      integer me,nscr,lscr,iscr,lcd_coef,icd_coef
c
      double precision Ebq
c
      integer g_rel_sf(*)       ! spin-free zora/x2c ga
      double precision Erel_sf  ! spin-free zora/x2c energy
c
      integer g_frozemb       ! frozen embedding ga
      double precision Efrozemb  ! energy

      integer g_cosmo

      logical no2e
c
c     timings
c
      double precision time1_ke,time2_ke
      double precision time1_pot,time2_pot
c  jbecca START-- making NWChem great again with efield
      double precision efield(3)
      double precision eelec
      logical lescale
c  jbecca END
c
c MN solvation models -->
c
      double precision gcds
      integer do_cosmo_vem
      logical do_cosmo_smd
c
c DIM/QM variables
c
      double precision g_vdim
      double precision Edim
      logical ldimqm, dimqm_go,use_dimqm
c
cc AJL/Begin/SPIN ECPs
      integer ecp_channels
cc AJL/Begin/FDE 
      integer i, j, iatom
      double precision Exc_fde(2)
      integer g_xcinv_fde
      integer g_vxc_fde(4)
c      double precision rho_n_fde ![out]
      integer g_dens_fde(2)
cc AJL/End
      double precision t0
      integer g_adft(2), g_nfock(max_quantum_nuc),g_nhcore
      logical ofirst_local, oreopen_local
      logical int2e_file_open, int2e_file_close
      external int2e_file_open, int2e_file_close
      double precision ecore_n, e_en, e_nn
      integer g_tmp_n

      if (.not. (owritefile .or. oreadfile)) then
         ofirst_local = .true.
         owritefile = .true.
         oreadfile = .false.
      else
         ofirst_local = .false.
         owritefile = .false.
         oreadfile = .true.
      endif
c
      do_cosmo_smd = .false.
      do_cosmo_vem = 0
      gcds = 0.D0
c
c
      ldimqm=.false.
      dimqm_go=.false.
c     =========================================================================

      me=ga_nodeid()
      if (odftps) call pstat_on(ps_fockb)

c ... jochen: to run tests on one-electron systems, we
c     simply disable all two-electron contributions if
c     requested by an input switch
      if(.not.rtdb_get(rtdb,'dft:no2e',MT_LOG,1,no2e))
     &   no2e = .false.
      if (no2e.and.me.eq.0) write (luout,*)
     &   'WARNING: dft:no2e set. Disabling Coulomb and XC Fock terms'


c
      Ecoul  = 0d0
      Exc(1) = 0d0
      Exc(2) = 0d0

      icd_coef = 1

c Jorge Garza/Begin
      if (ldmix) then
        iter_ini = 1
      else
        iter_ini = 3
      end if
      if (test_sic.eq.2.or.test_sic.eq.4) then
        if (iter.ge.iter_ini) then
           call xc_sicdeg(i_degen, evals, n_levels,act_levels)
        end if
      end if
c Jorge Garza/End
c     
c     Accumulate core hamiltonian into Fock matrix; compute core energy
c     
      if (oprint_time)
     &     call dft_tstamp(' calling 1-e int ')
      if (odftps) call pstat_on(ps_onel)
      call ga_zero(g_fock)
c
      if(util_print('dft timings', print_high))
     &      time1_ke=util_cpusec()   ! start kinetic time
      call int_1e_ga(ao_bas_han, ao_bas_han, g_fock,'kinetic',oskel)
      if(util_print('dft timings', print_high))
     &      time2_ke=util_cpusec()   ! end kinetic time
c
c     print kinetic energy matrix time
      if(util_print('dft timings', print_high)) then
        if (me.eq.0) then
           write(*,"(4x,'Kinetic Matrix Time:',F13.1,'s')") 
     &        time2_ke-time1_ke
        endif
      end if
c
      if(util_print('kinetic_energy', print_high)) then
         Ecore = ga_ddot(g_dens(1), g_fock)
         if(ga_nodeid().eq.0) then
            write(luout,222) Ecore
 222        format(/
     &     '           Kinetic energy =', f20.12/)
         endif
      endif
c
cc AJL/Begin/SPIN ECPs
cc    Read this value from rtdb vvvv
      if (.not.rtdb_get(rtdb, 'dft:spin_polarised_ecps',
     &  mt_int, 1, ecp_channels)) then
        ecp_channels = 1
      endif

      if (ecp_channels.eq.2.and.ipol.eq.1)
     & call errquit('dft_fockbld: spin_polarised_ecps with mult = 1', 
     &               0, BASIS_ERR)

c     == Copy Fock matrix to individual components for spin ECPs ==
      call ga_copy(g_fock, g_focks(1))
      if (ecp_channels.gt.1) call ga_copy(g_fock, g_focks(2))

      e_en = 0d0
      do i=1,nuc_num
         call ga_zero(g_fock)
         if (ofirst_local) then
            oreadfile = .false.
            owritefile = .true.
         end if
         call get_j_dm(geom,bas_e,bas_n(i),-1.0d0,tol2e_cneo,.false.,
     $         g_ndens(i),g_fock,'enints',i,0,.true.)

         e_en = e_en+ga_ddot(g_fock,g_dens(1))

         call ga_dadd(1d0,g_fock,1d0,g_focks(1),g_focks(1))
      enddo

      call util_file_name('aoints',.true.,.false.,int2efilename)
      if (ofirst_local) then
         oreopen_local = .false.
         owritefile = .true.
         oreadfile = .false.
      else
         oreopen_local = .true.
         owritefile = .false.
         oreadfile = .true.
      end if

      if (.not. int2e_file_open(int2efilename,0,0,tol2e_cneo,
     $      oreopen_local))
     &      call errquit('cneo_fock: open aoints failed',0,0)
      call cneo_schwarz_set(ao_bas_han)
      call cneo_reset_int2efile_print(print_debug)
c
c     if requested separate calculation of bq-electron energy (MV)
c     ----------------------------------------------------------- 
      if (rtdb_get(rtdb, 'dft:bq_energy', mt_dbl, 1,Ebq)) then
         call errquit('cneo_dft_fockbld:unsupported',0,0)
      else
        if(util_print('dft timings', print_high))
     &        time1_pot=util_cpusec()   ! start time
cc AJL/Begin/SPIN ECPs
c        call int_1e_ga(ao_bas_han, ao_bas_han, g_fock,
        call int_1e_ga(ao_bas_han, ao_bas_han, g_focks(1), 
     &                 'potential', oskel)
        if(util_print('dft timings', print_high))
     &        time2_pot=util_cpusec()   ! end time
c
c     print potential matrix time
        if(util_print('dft timings', print_high)) then
         if (me.eq.0) then
           write(*,"(4x,'Potential Matrix Time:',F13.1,'s')")
     &        time2_pot-time1_pot
         endif
        end if
      end if
c

c     == Copy Fock matrix to individual components for normal calcs ==
      if (ecp_channels.eq.1.and.ipol.gt.1) 
     &  call ga_copy(g_focks(1), g_focks(2))

      call ga_sync()
c
cc      Ecore = ga_ddot(g_dens(1), g_fock)
cc AJL: Ecore = ga_ddot(g_dens(1), g_focks(1))
cc      write(6,*)'Ecore(1:T):',ga_ddot(g_dens(1), g_focks(1))

cc    With spin-polarised ECPs Ecore will be spin dependent
cc    See Szabo and Ostlund pg. 215
cc        
      if (ecp_channels.gt.1) then
c     
c       call ga_sync()   
      else

        Ecore = ga_ddot(g_dens(1), g_focks(1))
c        write(6,*)'Ecore(1:T):',ga_ddot(g_dens(1), g_focks(1)) 
 
      endif
cc AJL/End
c
c     add Bonacic-Fantucci repulsive term
c
      fant_d=-1d0
      fant_a=-1d0
      if (odftps) call pstat_off(ps_onel)
c
c     Pre-compute reduced total density matrices over atoms
c 
      nscr = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if (.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &   call errquit('dft_scf: cannot allocate scr',0, MA_ERR)
      call util_ga_mat_reduce(nbf_ao, natoms, cetobfr, g_dens, 
     &                     ipol, rdens_atom, 'absmax', 
     &                     dbl_mb(iscr), nbf_ao_mxnbf_ce,.true.)
      if (.not.ma_pop_stack(lscr))
     &   call errquit('dft_scf: cannot pop stack',0, MA_ERR)
c
      if (ipol.gt.1)then
c     
c        Restore alpha and beta densities.
c     
         call ga_dadd(1d0, g_dens(1), -1d0, g_dens(2), g_dens(1))
c     
c        Note that g_dens(1) now contains the alpha density
c        matrix and g_dens(2) contains the beta
c
c        Pre-compute reduced alpha and beta density matrices over atoms
c 
         call dfill(ipol*natoms*natoms, 0.0d0, rdens_atom, 1)
         nscr = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
         if (.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &      call errquit('dft_scf: cannot allocate scr',0, MA_ERR)
         call util_ga_mat_reduce(nbf_ao, natoms, cetobfr, 
     &                       g_dens, ipol, rdens_atom, 
     &                       'absmax', dbl_mb(iscr), nbf_ao_mxnbf_ce,
     A        .true.)
         if (.not.ma_pop_stack(lscr))
     &      call errquit('dft_scf: cannot pop stack',0, MA_ERR)
      endif
c
c     DIM/QM JEM
c     Calculate DIM potential and energy
c     We calculate the DIM potential at each point from the total density
c     This potential is used to perturb both alpha and beta spins, so only
c     one call to dim_grid_quadv0 is needed
      if (ldimqm .and. dimqm_go) then
        if(.not.ga_duplicate(g_vxc(1), g_vdim, 'dim:1'))
     $    call errquit("dft_fockbld: can't create dim:1", 99, GA_ERR)
        call ga_zero(g_vdim)
        call dim_grid_quadv0(rtdb, g_vdim)
c        call ga_dadd(1.d0, g_vxc(1), 1.d0, g_fock, g_fock)
        Edim = ga_ddot(g_dens(1), g_vdim)
        if(ipol > 1) then
          Edim = Edim + ga_ddot(g_dens(2), g_vdim)
        end if
        call ga_sync()
      end if
c      
      if (CDFIT)then
         call errquit('cneo_dft_fockbld:unsupported',0,0)
      else
         Ecoul = 0d0
      endif
c     
c     
c     == Compute the XC potential and energy ==
      g_vxc(1) = g_tmp
      call ga_zero(g_vxc(1))
      rho_n = 0.0d0
cc AJL/Begin
cc Can't remember why I implemented this? To be tested.
cc      if (ipol.gt.1) call ga_zero(g_vxc(2))
      if (ipol.gt.1) call ga_copy(g_vxc(1), g_vxc(2))
cc AJL/End
      if (oprint_time)call dft_tstamp(' Before call to GETVXC.  ')
c
      if (odftps) call pstat_on(ps_xc)
c
c     == Calculate the HF exchange (exact Coulomb), DFT exchange-correlation ==
c     ! ADFT: pass fitting coefficients to build fitted density on the grid
      call xc_getv(rtdb, Exc, Ecoul,nExc, iVxc_opt, g_xcinv, g_dens, 
     &             g_vxc, IOLGC, g_wght, g_xyz, g_nq,wght_GA, rho_n,
     &             rdens_atom, cetobfr, natoms, dbl_mb(icd_coef))
c
      do i=1,ipol 
        call ga_dadd(1.d0, g_vxc(i), 1.d0, g_focks(i), g_focks(i))
      enddo


      if (odftps) call pstat_off(ps_xc)
c
      if (oprint_time)
     &     call dft_tstamp(' Before call to GETVCOUL.')
      if (odftps) call pstat_on(ps_vcoul)

      ! Compute contribution to Fock matrix from fitted
      ! quantities
      call dft_getvc(Dbl_MB(icd_coef), dd3cERI, Ecoul,
     &               g_vxc, iVcoul_opt, n_batch, 
     &               n3c_int, ii3cwhat, n3c_dbl,
     &               IOLGC,
     &               .false., ipol)
      if (odftps) call pstat_off(ps_vcoul)
c     
c     Add V coul to Fock Matrix
c     
c     add relativistic energy contributions
c     == if frozen embedding ==
      Efrozemb = 0.d0
c     == Add exact pot contributions ==
      if (exact_pot.eq.1) then
        call xc_exact_pot(Ecoul_aux, g_dens, g_vxc, g_rho_exact,
     &                    dif_lamda, g_tmp_exact, g_addit_exact)
      end if
c
c  Jorge Garza/Begin
      if (ipol.eq.2.and.noc(2).eq.0) call ga_zero(g_vxc(2))
      totsic=0.0d00
      if (test_sic.eq.2.or.test_sic.eq.4) then
        if (iter.ge.iter_ini) then
          if (test_sic.eq.4) then
            call localize_sic(g_movecs, evals)
          end if
          call xc_sic(rtdb, nExc,iVxc_opt, g_dens, 
     &                g_vxc,  g_wght, g_xyz, g_nq, wght_GA,
     &                rdens_atom, cetobfr, natoms,
     &                g_movecs, totsic, i_degen, n_levels)
          Exc(1) = Exc(1) + totsic
        end if
      end if
c  Jorge Garza/End
c
cc AJL/Begin/SPIN ECPs
cc     == Copy Fock matrix to individual components ==
c      call ga_copy(g_fock, g_focks(1))
c      if (ipol.gt.1) 
c     &     call ga_copy(g_fock, g_focks(2))
cc AJL/End
c
c     == Add exchange-correlation to the Fock matrix ==       
      call ga_dadd(1.d0, g_vxc(1), 1.d0, g_focks(1), g_focks(1))
      if (ipol.gt.1) 
     &     call ga_dadd(1.d0, g_vxc(2), 1.d0, g_focks(2), g_focks(2))

c jbecca START - damnit, quit taking this out
C
C     Compute -D.E contribution due to external applied electric field
C     (if dft:efield is stored in rtdb).  Note that we can potentially
C     apply different fields to the alpha and beta spins this way.
C
C

c
c     DIM/QM JEM
c     Add DIM potential to the Fock matrix
      if (ldimqm .and. dimqm_go) then
c  jbecca START: scaling potential here to test some stuff.
c              originally 1.d0
        call ga_dadd(1.d0, g_vdim, 1.d0, g_focks(1), g_focks(1))
c       DIM potential is from total density, so same for both alpha and beta
        if (ipol.gt.1)
     &    call ga_dadd(1.d0, g_vdim, 1.d0, g_focks(2), g_focks(2))
        if (.not. ga_destroy(g_vdim))
     &    call errquit("destroy vdim failed", 9, GA_ERR)
      end if


      if (oprint_time)
     &     call dft_tstamp(' End of parallel region. ')
c     
      if (odftps) call pstat_off(ps_fockb)
      if (.not. int2e_file_close(.true.))
     $       call errquit('cneo:get_j_dm: close file',0,0)
      call cneo_reset_int2efile_print(print_default)
c
cc AJL/Begin/SPIN ECPs (Debug)
c      if(util_print('Intermediate Fock Matrix', print_debug)) then
c        if (ga_nodeid() .eq. 0) then
c          write(6,*)
c          call util_print_centered(6,
c     $         'Intermediate Fock Matrix',40,.true.)
c          write(6,*)
c          call util_flush(6)
c        end if
c        call ga_print(g_focks(1))
c        if (ipol.eq.2) call ga_print(g_focks(2))
cc AJL/End
c
c


      if (ipol .gt.1)
     $   call ga_dadd(1d0,g_dens(1),1d0,g_dens(2),g_dens(1))

      e_nn = 0d0
      ecore_n = 0d0

      if (.not. ga_duplicate(g_nfock(1),g_tmp_n,'tmp n')) call errquit
     $      ('fockbld:ga_duplicate?',0,0)

      do i=1,nuc_num
         iatom = nuc_cent(i)
         call ga_zero(g_nfock(i))
         call get_j_dm(geom,bas_n(i),bas_e,-1.0d0,tol2e_cneo,.false.,
     $         g_dens(1), g_nfock(i), 'neints', i,0,.true.)

         do j=1,nuc_num
            if (j.ne.i) then
               if (ofirst_local .and. (j.gt.i)) then
                  oreadfile = .false.
                  owritefile = .true.
               end if
               call ga_zero(g_tmp_n)
               call get_j_dm(geom,bas_n(i),bas_n(j),1.0d0,
     $                  tol2e_cneo,.false.,g_ndens(j),g_tmp_n,'nnints',
     $                  i,j,.true.)
               call ga_dadd(1d0,g_tmp_n,1d0,g_nfock(i),g_nfock(i))
               e_nn = e_nn+ga_ddot(g_tmp_n, g_ndens(i))
            endif
         enddo
         call hcore_nuc(geom, bas_n(i), g_nhcore, iatom)
         call ga_dadd(1d0, g_nhcore, 1d0, g_nfock(i), g_nfock(i))
         if (.not. ga_destroy(g_nhcore)) call errquit
     $         ('fockbld:ga_destroy?',0,0)
         ecore_n = ecore_n+ga_ddot(g_nfock(i), g_ndens(i))
         call ga_copy(g_nfock(i),g_fock0(i))
         call cneo_solve_constrain(i)
         call ga_dadd(1d0,g_fock0(i),1d0,g_fock_add(i),g_nfock(i))
      enddo

      if (util_print('cneo_fock_bld',print_debug)) then
         do i=1,ipol
            call ga_print(g_focks(i))
         enddo
         do i=1,nuc_num
            call ga_print(g_nfock(i))
         enddo
      endif

      if (ipol .gt. 1)
     $   call ga_dadd(1d0,g_dens(1),-1d0,g_dens(2),g_dens(1))

      if(.not.ga_destroy(g_tmp_n)) call errquit
     +   ('cneo_dft_fockbld:destroy g_tmp_n?',0,GA_ERR)


      end
