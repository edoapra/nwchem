c
c     main scalar DFT driver
c
      logical function cneo_dft_scf(rtdb, Etold, Enuc, iVcoul_opt, 
     &           iVxc_opt, iter, g_dens, g_dens_at, g_movecs, g_vxc,
     &           g_fock, g_svecs, svals, g_xcinv, g_s, g_s_n,
     &           lcdft, nconstr, g_svecs_n, svals_n)
c     
c     $Id$
c     
      implicit none
#include "errquit.fh"
! #include "nwc_const.fh"
c
      integer rtdb              ! [input]
      double precision Etold, Enuc
      integer iVcoul_opt
      integer iVxc_opt
      integer iter
      integer g_dens(2), g_movecs(2), g_vxc(4), 
     &        g_fock, g_svecs,
     &        g_xcinv
      integer g_frozemb
c
      integer g_dens_at(2) 
c
      double precision  toll_s
      double precision svals(*), svals_n(*) ! [in] S evals
c     
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "cdft.fh"
#include "oep.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "dft_conv.fh"
#include "dftps.fh"
#include "cosmo.fh"
#include "zora.fh"
#include "x2c.fh"
#include "modelpotential.fh"
#include "inp.fh"
#include "ccneo.fh"
#include "cscf_cneo.fh"
#include "cconstrain.fh"
! #include "cscf_cneo.fh"
c     
      Logical movecs_write, movecs_converged, movecs_read
      External movecs_write, movecs_converged, movecs_read
c
      external dft_lindep
      integer dft_lindep
c
      integer g_s, g_s_n(max_quantum_nuc) ! [in]
      integer g_svecs_n(max_quantum_nuc)
c
      logical use_nwxc,use_dimqm
      logical lvan, do_ddl, lbrak, check
      double precision dl_conv, dl2, ddl2, x1, x2, fddl, dlold
      logical lcdft
      integer g_fockv(2)
      integer nconstr, ipop, counter_cdft, cdft_maxiter
      integer iatms1(max_constr), fatms1(max_constr)
      integer iatms2(max_constr), fatms2(max_constr), dtype(max_constr)
      integer g_constr(max_constr), g_shalf
      integer l_mdiis, k_mdiis, l_hess, k_hess, l_conscr, k_conscr
      double precision constr(max_constr), multipliers(max_constr),
     $            multipliers0(max_constr)
      double precision grad_constr(max_constr)
      double precision gold_constr(max_constr),p_constr(max_constr)
c
      integer ixyz, lxyz
      logical oprint
      double precision rms(max_quantum_nuc+2), derr(max_quantum_nuc+2)
      integer nmo(2), icall(max_quantum_nuc+2), nprint
      integer mclosed(2), mopen(2)  ! used for cosmo_charges: used in same way as scf does
      integer n3c_dbl, n3c_int, n_batch
      integer l_3cwhat, k_3cwhat, l_3cERI, k_3cERI
c
cc    Jorge Garza/Begin
c
c     Variables used to obtain the exact exchange-correlation potential
c
      integer g_rho_exact(2), g_tmp_exact(2), g_addit_exact(2),
     &        g_movecs_aux(2)
      integer incre, exact_pot, split, n_levels(2), act_levels
      double precision lamda, p_lamda, dif_lamda, lamda_old,
     &                 Ecoul_aux(2)
c
      integer test_sic, condfukui, l_degen, i_degen(2)
      double precision totsic
cc    Jorge Garza/End
c     AOR begin
      double precision exdm
      integer ixdm_v, ixdm_a, ixdm_ml
c     AOR end
      Integer l_eval
      integer k_eval(2)
      integer natoms, nTotEl
      integer l_occ, k_occ
      integer i
      integer me, nproc
      integer g_tmp, g_tmp2, g_fockt,  g_wght, g_xyz,g_nq
      integer g_tmp_n
      integer nheap, nstack
      integer ispin, idone
      integer npol
      integer leneval
      integer lcntoce, icntoce, lcntobfr, icntobfr,
     &        lcetobfr, icetobfr, lrdens_atom, irdens_atom
      double precision start_wall, current_wall, elapsed_wall,
     &                 save_wall, current_cpu, start_cpu,
     &                 wall_time_reqd
      integer int_wall_time_reqd
      double precision pp
      double precision anucl_charg, anel
      double precision etnew, tol2e, tol2e_sleazy,tdots, edisp
c     
c     Note, damping, levelshifting, and diising logicals
c     are used to turn on/off these procedures per
c     iteration.  The alternative logicals nodamping, 
c     nolevelshifting, and nodiis are specified and held
c     for the entire convergence sequence.
c     
      Logical  IOLGC, mulliken
      logical converged, wght_GA
      logical oprint_parm, oprint_conv, oprint_vecs, 
     &        oprint_eval, oprint_syma, oprint_time, 
     &        oprint_info, oprint_tol, oprint_final_vecs, 
     &        oprint_energy_step, oprint_intermediate_fock,
     &        oprint_3c2e, oprint_interm_overlap, oprint_interm_S2,
     &        oprint_conv_details, oprint_sic,oprintinertia
      double precision zero, onem, one
      parameter(zero = 0.d0, one = 1.d0, onem = -one)
c
      integer ilo, ihi          ! For printing movecs analysis
c
      double precision dft_dencvg, dft_time
      external dft_dencvg
      integer l_ir, k_ir
cc AJL/Begin/FDE
c      logical last_time_energy,noscf,frozemb
      logical last_time_energy,noscf
cc AJL/End
      logical fon, reorth
      character*7 name
      character*4 scftype
      character*7 vecs_or_dens
      character*255 basis_name, basis_trans,blob

      integer nopen, nclosed, g_denso(2),k_tmp,l_tmp,
     .     l_gdiis,k_gdiis,ipolmod, g_ndenso(max_quantum_nuc)
      logical dft_mem3c
      external dft_mem3c
c !!! BGJ
      logical cphf_poliz, do_poliz
      external cphf_poliz
c !!! BGJ
c
c     ----- cosmo -----
c
      logical odbug
      logical ocosmo_got_gasphase
      double precision egas
      double precision esol
      double precision ecosmo ! intra cosmo charges energy
      logical cosmo_esp
      character*(255) cosmo_file
c
c     slow switching of annealing
c
      double precision smear_in,ss,dft_rabuck,pstrace
      integer n_dep
      integer n_rabuck
      logical rabuck_act,spinset
      external dft_rabuck
c
      integer g_focks(2)
      double precision nel_fon(4)
      integer nmo_fon(4), ncore_fon(2)
      integer g_movad,nmo_adapt(2),
     ,     ncanorg,g_gmovecs(2)
      logical gotsmear,quickguess,diisreset
      logical llldb,llldb_out,util_statldb,staticguess
      external util_statldb
      integer ndisrst
c
c     vdw contrib
      double precision dum
      logical xc_chkdispauto
      external xc_chkdispauto
      logical disp
c
c     == zora related ==
      logical dft_zora_read, dft_zora_write, dft_zora_inquire_file
      external dft_zora_read, dft_zora_write, dft_zora_inquire_file
      character*255 zorafilename
      integer g_zora_scale_sf(2)
      double precision ener_scal
      integer icalczora
      logical ofinite,status,Knucl ! FA-added-09-26-11
c
c     virtual spectrum
      logical lvspec
      integer owstart(2), owend(2)
      character*255 vspecfilename
c
c     timings 
      double precision time1,time2

c     for Gaussian Nuclear Model
      integer iatom,l_zetanuc, k_zetanuc
      external get_zetanuc_arr
      integer dft_nonzvecs
      external dft_nonzvecs
c -----FA: 04-21-11: for Gaussian Nuclear Model --- END

c ... jochen: MO analysis print threshold
      double precision tanalyze
cc     AJL/Begin/FDE
      integer g_vxc_fde(4)
      integer g_xcinv_fde
      integer g_dens_fde(2)
      logical ldimqm, dimqm_go
      double precision edimqm
      double precision eelec
cc     AJL/End
c
      gotsmear(ss)=abs(ss).gt.1d-9
c
c ... exact decoupling and other relativistic Hamiltonians
      integer g_rel(2)
      integer  g_x2c_u(2), g_x2c_u2c, g_x2c_c2u
c
c     misc
c
      logical debug_fon
      logical do_zora_recalc, zora_recalc
      external zora_recalcc
      logical x2c_noatm
      external x2c_noatm

c MN solvation models -->
c
      double precision gcds, gspol, espol, gspoldyn, gspolneq
      double precision gstote, estote
      integer do_cosmo_vem, istep_cosmo_vem
      logical do_cosmo_smd

C     RTDB flag for printing out all matrices
      logical lprint_mats  !print matrices to stdout?
      logical lsave_mats   !save matrices to file?
      logical molden, forceatguess, save_evals, densmat
      integer back_to_atguess
      character*256 file_mat1,movecs_in_org
      integer ma_type, no_sflip,l_at_flip,k_at_flip
      character*26 date

      integer g_ndens(max_quantum_nuc), g_odens
      integer g_nfock(max_quantum_nuc)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      double precision ecore_n, e_en, e_nn, sval_start
      character*256 movecs_out_n
      character*5 str_i

c     == initialize DIIS call counter ==
      icall = 0
c
      gcds = 0.d0
      do_cosmo_vem = 0
      istep_cosmo_vem = 0
      do_cosmo_smd = .false.
c
c <-- MN solvation models
c
c     =========================================================================
c
c
      k_3cERI  = 0
      k_3cwhat = 0
c
      call ecce_print_module_entry('dft')
      cneo_dft_scf = .false.
      ocosmo_got_gasphase = .false.
      lmaxov_sv = lmaxov
      oprint = util_print('information', print_low)
      oprint_info = util_print('common', print_debug)
      oprint_parm = util_print('parameters', print_default)
      oprint_3c2e = util_print('3c 2e integrals', print_default)
      oprint_conv = util_print('convergence', print_low)
      oprint_conv_details = util_print('convergence details', 
     &                                 print_high)
      oprint_vecs = util_print('intermediate vectors', print_high)
      oprint_eval = util_print('intermediate evals', print_high)
      oprint_syma = util_print('interm vector symm', print_high)
      oprint_time = util_print('dft timings', print_high)
      oprint_tol = util_print('screening parameters', print_high)
      oprint_energy_step = util_print('intermediate energy info',
     &                                print_high)
      oprint_intermediate_fock = util_print('intermediate fock matrix',
     &                                      print_high)
      oprint_interm_S2 = util_print('intermediate S2',print_high)
      oprint_interm_overlap = util_print('intermediate overlap',
     &                                      print_high)
      oprint_final_vecs = util_print('final vectors', print_high)
      oprint_sic = util_print('SIC information', print_high)
      oprintinertia = util_print('inertia', print_default).and.
     . ga_nodeid().eq.0
      odbug = util_print('cosmo', print_high)

      me = ga_nodeid()
      nproc = ga_nnodes()

c     ---------------------------------------------------------
c     check which relativistic Hamiltonian (if any) we're using
c     and write some information to the output.
c     comment: Niri's zora implementation added an array g_zora_sf
c     to the Fock matrix in dft_fockbld. This array contains
c     T(rel) - T(nrel). The nonrelativistic kinetic energy is always
c     added to the Fock matrix in dft_fockbld. 
c     I changed the array name for g_zora_sf to g_rel and use
c     it now for both methods. 
c     ---------------------------------------------------------


c     make sure we don't use zora, x2c, or another
c     two-component relativistic method at the same time


      call x2c_scf_init(rtdb)

c
c !!! BGJ
c     Store SCF hamiltonian type as DFT for use in BGJ routines
      if (.not. rtdb_put(rtdb, 'bgj:scf_type', MT_INT, 1, 2))
     $     call errquit('dft_scf: put of bgj:scf_type failed',0,
     &       RTDB_ERR)
c !!! BGJ
c
c     check if smearing is on
c
      ssmear=0.d0
      delta=1d99
      ! derr(1)=1d99
      ! derr(2)=1d99
      derr = 1d99
      n_rabuck=0
      rabuck_act=.false.
      if (.not.rtdb_get(rtdb, 'dft:diisreset', mt_int, 1,
     &   ndisrst))ndisrst=9999
cold  diisreset=.false.
      if(ndisrst.gt.0) diisreset=.true.
      if (.not.rtdb_get(rtdb, 'dft:spinset', mt_log, 1,
     &   spinset))spinset=.true.
c
c     check for FON in rtdb:
      fon=.false.
      debug_fon = .false.

c     read MO analysis threshold
      if (.not.rtdb_get(rtdb,'movecs:tanalyze',MT_DBL,1,tanalyze))
     &   tanalyze = 0.025d0


c
      if (.not.rtdb_get(rtdb, 'dft:reorth', mt_log, 1,
     &   reorth))reorth=.false.
c
c     see if levelshifting monitoring is desired
c
      if (.not. rtdb_get(rtdb, 'dft:check_shift', mt_log, 1,
     &   check_shift))then
         check_shift = .false.      
      endif
c     
      if (.not. geom_ncent(geom, natoms))
     &   call errquit('dft_scf: geom_ncent failed',73, GEOM_ERR)
      if (.not. geom_nuc_charge(geom, anucl_charg))
     &   call errquit('dft_scf: geom_nuc_charge failed', 0, GEOM_ERR)
c     
      anel = int(anucl_charg) - rcharge
c
c     Pre-compute mapping vectors
c 
      if (.not.ma_push_get
     &   (mt_int,nshells_ao,'cntoce map',lcntoce,icntoce))
     &   call errquit('dft_scf:push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     &   (mt_int,nshells_ao*2,'cntoce map',lcntobfr,icntobfr))
     &   call errquit('dft_scf:push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     &   (mt_int,natoms*2,'cntoce map',lcetobfr,icetobfr))
     &   call errquit('dft_scf:push_get failed', 13, MA_ERR)
c     
      call build_maps(ao_bas_han, int_mb(icntoce), int_mb(icntobfr), 
     &                int_mb(icetobfr), natoms, nshells_ao)

      if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, noc))
     &   call errquit('dft_scf: rtdb_get of noc failed', 0, RTDB_ERR)
      if (.not. MA_Push_Get(MT_Dbl,nbf_ao*ipol,'mo occ',l_occ,k_occ))
     &   call errquit('dft_scf: failed to alloc',999, MA_ERR)
c
c     get/set orbital overlap tolerancec

      call util_tolls(rtdb,oprint_conv_details,toll_s,n_dep,
     .     svals,nbf_ao)
      nmo(1)=nbf_ao
      if (ipol.eq.1)then
c     
c        noc(1) = No. of electrons alpha
c     
         nTotEl = 2*noc(1)
         mclosed(1)=noc(1)
         mopen(1)=0
      else
c     
         nTotEl = noc(1) + noc(2)
         nmo(2) = nmo(1)
         mclosed(1)=0
         mclosed(2)=0
         mopen(1)=noc(1)
         mopen(2)=noc(2)
      endif
c
c     AOR begin
      call xc_xdm_init(rtdb,ixdm_v,ixdm_ml)
c     AOR end
c    Jorge Garza/Begin
      call xc_sicinit(rtdb, test_sic, condfukui, exact_pot,
     &                l_degen, i_degen, noc, act_levels)
c     Jorge Garza/End
      use_nwxc = util_module_avail("nwxc")
      if (use_nwxc) then
         call nwxc_getvals("nwxc_is_on",use_nwxc)
      endif
      use_dimqm = util_module_avail("dimqm")
      if (use_dimqm) then
         call dimqm_init(rtdb)
         call dimqm_used(ldimqm)
         call dimqm_getgo(dimqm_go)
      else
         ldimqm=.false.
         dimqm_go=.false.
         edimqm=0d0
      endif
c
c     Set aside some memory for reduced density matrix
c
      ipolmod=ipol
      if(test_sic.ne.0) ipolmod=2
      if (.not.MA_Push_Get(MT_Dbl,ipolmod*natoms*natoms,'rdens_atom',
     &   lrdens_atom,irdens_atom))
     &   call errquit('dft_scf: cannot allocate rdens_atom',0, MA_ERR)
c     
c     determine pattern of orbitals' occupancy
c     
      if (ipol .eq. 1)then
c
c        RHF occupations
c
         call dfill(nbf_ao, 0.0d0, dbl_mb(k_occ), 1)
         do i = 1, noc(1)
            dbl_mb(i-1+k_occ) = 2.0d0
         enddo
      else  
c
c        UHF occupations
c
         call dfill(2*nbf_ao, 0.0d0, dbl_mb(k_occ), 1)
         do i = 1, noc(1)
            dbl_mb(i-1+k_occ) = 1.0d0
         enddo
         do i = nbf_ao+1, nbf_ao+noc(2)
            dbl_mb(i-1+k_occ) = 1.0d0
         enddo
      endif 
      do i=1,nuc_num
         call dfill(nbf_n,0.0d0,dbl_mb(k_occ_n(i)),1)
         dbl_mb(k_occ_n(i)) = 1.0d0
      enddo
c     
      wght_GA = .false.
c     
c     Determine whether to fit the electronic charge density.
c     
      CDFIT = .FALSE.
      XCFIT = .FALSE.
c     
c     Define various constants.
c     
      npol = (ipol*(ipol+1))/2
c     
      itol_max = itol2e
      iaoacc_max = iaoacc
      tol_rho_max = tol_rho
      if (oprint_time)
     &     call dft_tstamp(' Before 3c-2e initialize.')
c
      mulliken = .false.
      if (imull.eq.1)mulliken = .true.
      IOLGC = .TRUE.
      if (noio.eq.1)IOLGC = .FALSE.
c     
c     Energy decomposition switch
c     
      nExc    = idecomp + 1
      Etnew = 0.d0
      tdots = 0d0
c
c     vdw  bit
c
c     activate disp if is present in rtdb
c     or include dispersion if this a functional that includes dispersion
c
      Edisp = 0.0d0
      if (.not.rtdb_get(rtdb, 'dft:disp', mt_log, 1, disp))
     &   disp=.false.
c
      if(disp.or.xc_chkdispauto())
     &      call xc_vdw(rtdb,geom,Edisp,dum,'energy')
c     
c     == scf energy convergence criterion ==
      if (.not. ga_duplicate(g_fock, g_tmp, 'tmp matrix'))
     &  call errquit('dft_scf:failed duplicate',g_fock, GA_ERR)
      if (.not. ga_duplicate(g_fock, g_focks(1), ' alpha fock'))
     &  call errquit('dft_scf:failed duplicate',g_fock, GA_ERR)
c
c MN solvation models -->
c
c     for VEM calculations if istep_cosmo_vem=2
c
      if (istep_cosmo_vem.eq.2) nodiis=.true.
      if (istep_cosmo_vem.eq.2) nolevelshifting=.true.
c
c <-- MN solvation models
c
c     == set up local convergence parameters ==
      diising = diis
      damping = damp
      levelshifting = levelshift
      keep_damp_on = .false.
      keep_levl_on = .false.
      keep_diis_on = .false.
      ndamp_input = ndamp
      rlshift_input = rlshift
      ndamp_def = 0
      rlshift_def = 0.0
      rlshift = rlshift_def
c     
      if (nodamping)damping = .false.
      if (nolevelshifting) then 
        levelshifting = .false.
        rlshift = rlshift_def
      endif
c
c     init diis handles to detect allocation
c
      k_gdiis=-1
      k_mdiis=-1
      k_hess=-1
      
      if (nodiis)then
         diising = .false.
      else
         if (.not.MA_alloc_Get(MT_int,nfock*(ipol+nuc_num)*2,
     ,        'ga handle for diis',l_gdiis,k_gdiis))
     ,        call errquit('dft_scf: cannot allocate g_diis',0, MA_ERR)
         if(nconstr.ne.0) then
            if (.not.MA_alloc_Get(MT_DBL,nfock*nconstr,
     ,           'ga handle for mdiis',l_mdiis,k_mdiis))
     ,           call errquit('dft_scf: cannot allocate mdiis',0,MA_ERR)
            if (.not.MA_alloc_Get(MT_DBL,max(nfock,nconstr)*nconstr,
     ,           'ga handle for hess',l_hess,k_hess))
     ,           call errquit('dft_scf: cannot allocate hess',0, MA_ERR)
         endif
      endif

      homo_lumo_gap = 200.0d0
      if (ncydp.ne.0)then
         damping = .true. 
         ndamp = ndamp_input
      endif
      if (ncysh.ne.0)then
         levelshifting = .true.
         rlshift = rlshift_input
      endif
      if (ncyds.ne.0)then
         diising = .true.
      endif
c     
c     
c     == Begin the SCF cycle ==
c     
c     == allocate eigenvalue array, including second pointer to beta block ==
      leneval = nbf_ao * ipol
      if (.not.MA_Push_Get(MT_Dbl,leneval,'eval',l_eval,k_eval(1)))
     &     call errquit('dft_scf: cannot allocate eval',0, MA_ERR)
      if (ipol .eq. 2)then  
         k_eval(2) = k_eval(1) + nbf_ao
      endif 
c     
c     == dump DFT parameters (if debugging) to see if they make sense ==
      if (me.eq.0.and.oprint_info)call dft_dump_info(me)

      if (ipol.eq.1)then
         scftype = 'RHF'
      elseif (ipol.eq.2)then
         scftype = 'UHF'
      endif
c
c     == needed for dftmp2 ==
      if (.not. rtdb_cput(rtdb,'scf:scftype', 1, scftype))
     &   call errquit('dft_scf: rtdb_cput scftype failed', 0, RTDB_ERR)
c
c     == allocate array for irreps ==
      if (.not.MA_Push_Get(mt_int,ipol*nbf_ao,'dft:irreps',l_ir,k_ir))
     &   call errquit('dft_scf: cannot allocate irreps',0, MA_ERR)
      call ifill(ipol*nbf_ao,1,int_mb(k_ir),1)
      nopen = mult - 1
      nclosed = (nTotEl - nopen) / 2
      if (.not. rtdb_put(rtdb,'scf:nclosed', mt_int, 1, nclosed))
     &   call errquit('dft_scf: rtdb_put nclosed failed', 0, RTDB_ERR)
c
      if (.not. bas_name(ao_bas_han, basis_name, basis_trans))
     $     call errquit('dft_scf: bas_name?', 0, BASIS_ERR)
c
c     == get info for int2e_ and set sleazy tolerance ==
      if (.not.rtdb_get(rtdb,'dft:tol2e_sleazy',
     ,     mt_dbl,1,tol2e_sleazy)) tol2e_sleazy = 1.d-4
      call scf_get_fock_param(rtdb, tol2e_sleazy)
c
c     == force sleazy SCF into "direct" mode ==
      call fock_force_direct(rtdb)
c
c     -------------------------------------------
c     calculate superposition of atomic densities
c     if ZORA or X2C functionality is requested
c     -------------------------------------------

c
c     == get atomic zora contributions ==
      do_zora = .false.
c
c
c     == determine guess ==
      if(.not.rtdb_get(rtdb,'back_to_atguess',mt_int,1,
     N     back_to_atguess)) 
     .     back_to_atguess=0
      if(.not.rtdb_get(rtdb,'forceatguess',mt_log,1,forceatguess)) 
     .     forceatguess=.false.
 1993 if(forceatguess) movecs_in='atomic'
      if(movecs_in.eq.'atomic') back_to_atguess=0
      call dft_guessin(movecs_in,ldmix,ncanorg,fon,
     ,     vecs_or_dens,
     ,     ipol,nbf_ao,g_movecs,g_gmovecs,
     ,     toll_s,svals)
      if(.not.rtdb_get(rtdb,'dft:densmat',mt_log,1,densmat))
     D     densmat=.false.
      if(densmat) goto 2022
c
c     == better have static ldb at high node counts ==
      if (odftps) call pstat_on(ps_guess)
      if (.not.rtdb_get(rtdb,'dft:staticguess',
     ,     mt_log,1,staticguess)) staticguess=.false.
      if(staticguess) llldb=util_statldb(.true.,rtdb)
      movecs_in_org=movecs_in
      call scf_vectors_guess(rtdb, tol2e_sleazy, geom, ao_bas_han, 
     &                       basis_trans, movecs_in, movecs_out, 
     &                       movecs_guess, scftype, nclosed, nopen, 
     &                       nbf, nmo, noc(1), noc(2),  k_eval, k_occ, 
     &                       k_ir, g_gmovecs, g_dens, vecs_or_dens, 
     &                       'dft', title, oskel, oadapt, 
     &                       .true.) 
c     check if movecs read has failed
      if(movecs_in_org.ne.movecs_in) then
         if(ga_nodeid().eq.0) then
            write(luout,*) ' WARNING: movecs_in_org=',
     c           movecs_in_org(1:inp_strlen(movecs_in_org)),
     c           ' not equal to movecs_in=',
     c           movecs_in(1:inp_strlen(movecs_in))
         endif
c     no movecs available, use svecs
         ncanorg=0
      endif
      if(staticguess) llldb_out=util_statldb(llldb,rtdb)
      call dft_guessout(nmo,nbf_ao,g_gmovecs,g_movecs,ipol)
      if (me.eq.0.and.oprint)
     &     write(LuOut,'(2x," Time after variat. SCF: ",f8.1)')
     &     util_wallsec()
      if (odftps) call pstat_off(ps_guess)
2022  continue
      if(.not.rtdb_get(rtdb,'quickguess',mt_log,1,quickguess)) 
     .     quickguess=.false.
      if(quickguess.and.movecs_in.eq.'atomic') then
        if(.not.rtdb_put(rtdb,'quickguess',mt_log,1,.false.)) 
     .       call errquit(' dftscf: rtbd_put failed ',0,RTDB_ERR)
        if (.not. rtdb_cput(rtdb, 'dft:input vectors', 1, movecs_out))
     $       call errquit('dft_scf: DFT MO vectors not defined',0,
     &                     RTDB_ERR)
        if (.not. rtdb_cput(rtdb, 'dft:output vectors', 1, movecs_out))
     $       call errquit('dft_scf: DFT MO vectors not defined',0,
     &                     RTDB_ERR)
      call dft_quickguess(ipol,nbf_ao,g_dens,oskel,
     N     ntotel,nclosed,
     D     densmat)
        goto 1789
      endif
      if(oskel) then
        call ga_symmetrize(g_dens(1))
        if(ipol.eq.2) call ga_symmetrize(g_dens(2))
      endif

      if(ldmix) then
        tdots=0.d0
        call dft_densm(g_dens, g_movecs, 
     &       nbf_ao, nmo, ipol, 
     &       geom, AO_bas_han,noc,ntotel,
     &       Dbl_MB(k_eval(1)),dbl_mb(k_occ),
     &       ssmear,tdots,iter,.false.,
     .       fon, nel_fon,nmo_fon,ncore_fon,
     .       spinset.or.n_rabuck.ne.0,
     &       rtdb)
c
c ... jochen: next line added '.or. fon' which
c       gives us the number of electrons created by the
c       fractional occupations:

        if(oprint_info .or. fon) then
          pstrace=ga_ddot(g_dens,g_s)
          if(ga_nodeid().eq.0) then
            write(luout,'(5x,a)') 'FON applied'
            write (luout,'(5x,a,1x,2pe15.7)')
     &         'tr(P*S): ',pstrace 
          end if
          if (.not. rtdb_put(rtdb, 'dft:pstrace', mt_dbl, 1, pstrace))
     &       call errquit('dft_scf: rtdb_put pstrace failed',
     &       1, RTDB_ERR)
        endif
c
c       Jorge Garza/Begin/Exact_pot
c
        if (exact_pot.eq.1) then
          call open_xc_exact_pot(g_rho_exact, g_tmp_exact,
     &                           g_addit_exact, g_movecs_aux,
     &                           g_dens, dif_lamda, lamda,
     &                           p_lamda, lamda_old, Ecoul_aux)
        end if
c
c       Jorge Garza/End/Exact_pot

        call dft_prevals('inter',me,oprint_eval,oprint_vecs,
     &                   g_movecs,dbl_mb(k_eval(1)),dbl_mb(k_occ))

      endif
c
c     Tidy up SCF
c
 1789 continue
      call fock_2e_tidy(rtdb)
c
c
c     set initial coulomb acc
c
c      write(LuOut,*)' movecs_guess = ',movecs_guess
      if (movecs_guess.eq.'restart') then
         ltight=.true.
         if (noc(1).gt.0) then
           homo_lumo_gap=dbl_mb(k_eval(1)+noc(1)) -
     -                   dbl_mb(k_eval(1)+noc(1)-1)
         endif
         if (noc(2).gt.0) then
           homo_lumo_gap=min(homo_lumo_gap,
     &                   dbl_mb(k_eval(2)+noc(2)) -
     -                   dbl_mb(k_eval(2)+noc(2)-1))
         endif
      endif
c     
c     May not want levelshifting initially until sure that the
c     transformed Fock matrix will be diagonally dominant, or
c     alternatively shift the piss out of it.
c     
      if (movecs_guess.eq.'restart'.or.ncanorg.eq.1.or.
     A     movecs_guess.eq.'fragment')then
         levelshifting = .not.gotsmear(ssmear)
      else
         levelshifting = .false.
      endif
      iswitc = 0
      if (ltight)then
         iAOacc = iAOacc_max
         tol_rho = tol_rho_max
         iswitc = 2
      else
         iAOacc = iAOacc_min
         tol_rho = max(tol_rho_min,tol_rho_max)
      endif
      if (use_nwxc) call nwxc_reset_rho_tolerance(tol_rho)
      if (ltight.or.(.not.direct))then
         itol2e = itol_max
      else
         itol2e = min(itol_min,itol_max)
      endif
c     
      tol2e = 10.d0**(-itol_max)
c 
c     Restore SCF parameters
c
      call scf_get_fock_param(rtdb, tol2e)
c
c     If open shell put the total density matrix in g_dens(1)
c
      if(ipol.eq.2)then             
         call ga_dadd(one,g_dens(1),one,g_dens(2),g_dens(1))
      endif
c     
c     Call to Mulliken Population Analysis for initial density
c     
      if (mulliken) call dft_mulwrap(me,g_dens,g_s)
c
c    The loop associated with the label 135 is to find the exact exchange-
c    correlation potential.
c
 135  continue
c     
c     Top of infinite SCF iteration loop
c
cc AJL/Begin/FDE
c
c     == if frozen embedding ==
c     AJL: This is already stored in commons
c      frozemb = .false.
c      if (.not.rtdb_get(rtdb, 'dft:frozemb', mt_log, 1, frozemb))
c     &    frozemb=.false.
c
c     Write prep time required
c
      if (me.eq.0.and.oprint)then
         current_wall = util_wallsec()
         write(LuOut,20)current_wall
   20    format(2x,' Time prior to 1st pass: ',f8.1)
      endif
c     
c     start DFT_SCF timer
c     
      start_wall = util_wallsec()
      start_cpu = util_cpusec()
      dft_time = -start_wall
c
      if (oprint_time)
     &     call dft_tstamp('   Before SCF iter loop. ')
c
      iter = 0
c
      if (ipol.gt.1)then
         if (.not. ga_duplicate(g_fock, g_fockt, 'fock tr'))
     &        call errquit('dft_scf: error creating ga',0,GA_ERR)
      endif

      if (ipol.gt.1)then
       if(.not.ga_duplicate(g_fock, g_focks(2), ' beta fock'))
     &    call errquit('dft_scf:failed duplicate',g_fock, GA_ERR)
      endif
 3000 continue
c
c     ----- dft scf proper -----
c
      last_time_energy = .false.
c
      if (.not. ga_duplicate(g_fock, g_denso(1), 'oldDMa'))
     &     call errquit('dft_scf: error creating ga',0,GA_ERR)
      if(ipol.eq.2) then
         if (.not. ga_duplicate(g_fock, g_denso(2), 'oldDMb'))
     &        call errquit('dft_scf: error creating ga',0,GA_ERR)
      endif
      do i=1,nuc_num
         if (.not. ga_duplicate(g_svecs_n(i),g_ndenso(i),'oldDMa:n'))
     &          call errquit('dft_scf: error creating ga',0,GA_ERR)
      enddo

      do i=1,nuc_num
         iatom = nuc_cent(i)
         g_ndens(i) = ga_create_atom_blocked(geom,bas_n(i),'n:dens')
         call ga_zero(g_ndens(i))
         call cneo_vecs_guess_nuc(rtdb,i,iatom,g_dens(1))
         g_nfock(i)=ga_create_atom_blocked(geom,bas_n(i),
     $                        'cneo:n:fock')
         call rohf_dens(g_movecs_n(i),nbf_n,1,0, g_ndens(i),g_odens)
         call ga_dscal(g_ndens(i), 5.0d-1)
      enddo
c
      if(iterations.eq.0) then
         if (.not. ga_destroy(g_denso(1))) call errquit
     &        ('dft_scf: could not destroy g_denso', 0, GA_ERR)
         if(ipol.eq.2) then
            if (.not. ga_destroy(g_denso(2))) call errquit
     &           ('dft_scf: could not destroy g_dens2', 0, GA_ERR)
            if(iter.gt.0) then
            if (.not. MA_free_heap(l_gdiis))  call errquit
     &           (' dft_scf:cannot popstack',111, MA_ERR)
            endif
         endif
         do i=1,nuc_num
            if (.not. ga_destroy(g_ndenso(i))) call errquit
     &        ('dft_scf: could not destroy g_denso', 0, GA_ERR)
         enddo
         if (odftps) call pstat_on(ps_scfend)
        !  goto 1970
      endif

c
 1000 continue !  iteration loop

c
      if (me.eq.0 .and. oprint_conv_details)
     &   write(LuOut,124)damping, levelshifting, diising
 124  format(10x,' DAMPING=',l1,' LEVELSHIFTING=',l1,
     &           ' DIISING=',l1)
c
      if (me.eq.0.and.oprint_tol)write(LuOut,3234)itol2e,iAOacc,iXCacc
 3234 format(10x,'itol2e=',i2,' iAOacc=',i2,' iXCacc=',i2)

      Ecoul  = ZERO
      Exc(1) = ZERO
      Exc(2) = ZERO
      rms(1) = 0.d0
      rms(2) = 0.d0
      rms(3) = 0.d0
c
c
c     DIM/QM JEM
c     Calculate DIM/QM dipoles
c     We use the total density to calculate the dipoles.
c
      if(ldimqm .and. dimqm_go) then
        call dimqm_main(rtdb, ao_bas_han, geom, g_dens, .false.)
      end if

c
c     Determine HOMO energy (for Zhan-Nichols-Dixon AC)
c
      e_homo = dbl_mb(k_eval(1)+noc(1)-1)
      if (ipol.eq.2) then
        if (e_homo.lt.dbl_mb(k_eval(2)+noc(2)-1))
     &  e_homo = dbl_mb(k_eval(2)+noc(2)-1)
      endif
c
c     check for fractional occupations (fon), calculate tr[P S] and print
c     if (debug_fon) call dft_pstrace(g_dens(1),ao_bas_han,nbf_ao,oskel)
c
c     build DFT Fock matrix
c  jbecca START: adding in electric field energy interaction term
      eelec = zero
c
      if(util_print('dft timings', print_high))
     &      time1=util_wallsec()  ! start fock matrix build time
c     flip spins if needed
      if (.not. rtdb_get_info(rtdb, 'nwdft:spinflip',
     C     ma_type, no_sflip, date)) no_sflip = 0
      if (no_sflip.ne.0.and.ipol.eq.2.and.iter.eq.0) then
         if (.not.MA_Push_Get(MT_int,no_sflip,'flip',
     F        l_at_flip,k_at_flip))
     &        call errquit('dft_scf: cannot allocate flip',0, MA_ERR)
         if (.not.rtdb_get(rtdb,'nwdft:spinflip',
     M        mt_int,no_sflip,int_mb(k_at_flip)))
     &        call errquit('dft_scf: read nwdft:spinflip failed ',111,
     &        UNKNOWN_ERR)
         call mull_pop(geom,ao_bas_han,g_dens(1),g_s,'DM1')
         call mull_pop(geom,ao_bas_han,g_dens(2),g_s,'DM2')
         call dft_spinflip(g_dens,g_s,ao_bas_han,1d0,-2d0,
     S        no_sflip,int_mb(k_at_flip))
         call mull_pop(geom,ao_bas_han,g_dens(1),g_s,'DM1')
         call mull_pop(geom,ao_bas_han,g_dens(2),g_s,'DM2')
         if (.not.ma_pop_stack(l_at_flip))
     &        call errquit('dft_scf: cannot pop stack',33, MA_ERR)
c     remove flipping to avoid multiple operation in geom opt.
         if(.not.rtdb_delete(rtdb, 'nwdft:spinflip'))
     &        call errquit('dft_scf: cannot rtdb_delete',3, MA_ERR)
      endif
      
      if (.not. ga_duplicate(g_nfock(1), g_tmp_n, 'tmp n matrix'))
     &  call errquit('dft_scf:failed duplicate',g_nfock(1), GA_ERR)


      call cneo_dft_fockbld(rtdb,natoms,ntotel, g_ndens,
     ,     g_fock,g_focks,g_dens,g_vxc,g_movecs,g_xcinv, 
     ,     g_tmp,g_nq,g_wght,g_xyz,
     ,     ivcoul_Opt,nexc,ivxc_opt,wght_GA,
     ,     n_batch,n3c_int,n3c_dbl,IOLGC,
     ,     k_eval, dbl_mb(irdens_atom),
     ,     dbl_mb(k_3ceri),int_mb(k_3cwhat),
     ,     int_mb(icetobfr),
     ,     iter, ecore, ecoul,exc,rho_n,ldmix,
     ,     test_sic,exact_pot,ecoul_aux,dif_lamda,totsic,
     ,     i_degen,g_rho_exact,g_tmp_exact, g_addit_exact,
     ,     n_levels, act_levels, oprint_time,
cc AJL/Begin/FDE
! c     &     g_rel,g_frozemb
     &     g_rel,g_frozemb,g_vxc_fde,g_xcinv_fde,g_dens_fde,eelec,
     &     g_nfock, ecore_n, e_en, e_nn)
cc AJL/End
      if(util_print('dft timings', print_high))
     &      time2=util_wallsec()  ! end fock matrix build time
c
c     print fock matrix build time
      if(util_print('dft timings', print_high)) then
        if (me.eq.0) then
          write(*,"(4x,'Total Fock Matrix Build Time:',F13.1,'s')") 
     &                     time2-time1
        endif
      end if

c     Calculate the total electronic energy.
c     
      if (nExc.eq.1)then
c  jbecca START: adding in eelec here
         Etnew = Ecore + Ecoul + Exc(1) + Edisp + eelec
     $                      + ecore_n - e_en - 0.5*e_nn
      else
c  jbecca START adding in eelec here
         Etnew = Ecore + Ecoul + Exc(1) + Exc(2) + Edisp + eelec
     $                      + ecore_n - e_en- 0.5*e_nn
      endif
c
      if (oprint_conv_details.and.me.eq.0)then
         write(LuOut,*)'Etnew, Ecore, Ecoul, Exc(1), Exc(2), Totsic: ',
     &                  Etnew, Ecore, Ecoul, Exc(1), Exc(2), Totsic
         write(LuOut,*)'   Ecore_n, Enn, Een: ',
     &                     ecore_n, e_nn, e_en
      endif
c
c
c     == is this the last energy evaluation ? ==
      if (last_time_energy) then
         if (.not. converged_solve) then
            write(6,*) 'CNEO:constrain not satisfied, continue scf'
            last_time_energy = .false.
            goto 100
         endif
         etnew=etnew+tdots
c     
c        If open shell put the total density matrix back in g_dens(1) and quit.
         if (ipol.eq.2)then
            call ga_dadd(one, g_dens(1), one, g_dens(2), g_dens(1))
         endif
         goto 2000
      endif  ! last_time_energy
 100  continue
c
      call dft_convpam(oprint_conv_details,iter,etnew,etold,enuc)
c
c
      if (diising)then
c     
c     DIIS step taken here.
c     
         if (oprint_time)
     &        call dft_tstamp(' calling diis driver ')
         call cneo_diis_driver(toll_s, derr,  ipol, 
     &           icall, nfock, nbf_e, nbf_n, geom, 
     &           ao_bas_han, g_focks, g_nfock,g_dens,
     &           g_ndens, g_tmp_n,
     &           g_tmp, g_svecs, svals, 
     &           g_svecs_n, svals_n, diising, 
     &           nodiis,int_mb(k_gdiis),lkeeps,
     .           g_s12m,g_s12p,.false.,
     $           nconstr,multipliers0,dbl_mb(k_mdiis))
         if (oprint_time)
     &        call dft_tstamp(' called diis driver ')
      endif

      if(util_print('Intermediate Fock Matrix', print_debug)) then
        if (ga_nodeid() .eq. 0) then
          write(6,*)
          call util_print_centered(6,
     $         'Intermediate Elec Fock Matrix',40,.true.)
          write(6,*)
          call util_flush(6)
        end if
        call ga_print(g_focks(1))
        if (ipol.eq.2) call ga_print(g_focks(2))
      endif
cc AJL/End

      do ispin = 1, ipol
         call dft_diagn(levelshifting,
     ,        nmo(ispin),ncanorg,
     .        g_focks(ispin),g_s,g_movecs(ispin),g_tmp,g_svecs,
     ,        dbl_mb(k_eval(ispin)),svals,noc(ispin),
     ,        homo,lumo,toll_s,
     ,        oprint_intermediate_fock,oprint_time,
c
c MN solvation models -->
c
     ,        oprint_conv_details,istep_cosmo_vem)
c
c <-- MN solvation models
c
c        determine homo-lumo gap 
c
         homo_lumo_gap = min(homo_lumo_gap, (lumo-homo-rlshift))
         if (me.eq.0 .and. oprint_conv_details)
     &      write(LuOut,4224)homo,lumo,rlshift, homo_lumo_gap
c     
         call ga_sync
c     
c        Save previous density for convergence check.
c     
c        call ga_copy(g_dens(ispin), g_denso(ispin))

c Constrained DFT /When lcdft, save at the beginning
      if((.not.lcdft) .or. dlold.eq.zero) then
         call ga_copy(g_dens(ispin), g_denso(ispin))
         ! call ga_copy(g_ndens, g_ndenso)
      endif

c Constrained DFT / end

c
c        symmetry adapt vectors?
c
         if (oadapt)then
            if(nmo(ispin).ne.0) then
               if (.not. ga_create(mt_dbl, nbf_ao, nmo(ispin),
     +           'movad', nbf_ao, 0, g_movad)) call errquit(
     +              'dft_scf: gacreate failed ',0, GA_ERR)
               call ga_copy_patch('n', g_movecs(ispin), 1, nbf_ao, 
     .            1, nmo(ispin),
     ,            g_movad, 1,nbf_ao,1,nmo(ispin))

            else
               g_movad=g_movecs(ispin)
            endif
            if(ispin.eq.1)name = '- alpha'
            if(ispin.eq.2)name = '- beta'
            call scf_movecs_sym_adapt(ao_bas_han, g_movad,
     &           oprint_syma, nmo, name,
     &           .true., 
     &                                int_mb(k_ir+nbf_ao*(ispin-1)))
            if(nmo(ispin).ne.0) then
               call ga_copy_patch('n', 
     ,              g_movad, 1,nbf_ao,1,nmo(ispin),
     .              g_movecs(ispin), 1, nbf_ao, 1, nmo(ispin))
               if (.not. ga_destroy(g_movad)) call errquit(
     &              'dft_scf: could not destroy g_movad', 0, GA_ERR)
            endif
         endif      
      enddo                     ! end big loop over ispin
      call cneo_set_cdft('n')
      do i=1,nuc_num
         
         sval_start = (i-1)*nbf_n+1
         call dft_diagn(levelshifting,
     ,         nmo_n,ncanorg,
     .         g_nfock(i),g_s_n(i),g_movecs_n(i),g_tmp_n,g_svecs_n(i),
     ,         dbl_mb(k_eval_n(i)),svals_n(sval_start),1,
     ,         homo,lumo,toll_s,
     ,         oprint_intermediate_fock,oprint_time,
c
c MN solvation models -->
c
     ,        oprint_conv_details,istep_cosmo_vem)

         call ga_copy(g_ndens(i), g_ndenso(i))
      enddo
      call cneo_set_cdft('e')

      call ga_sync
c     
c     Form a new density matrix.
c     
      tdots=0.d0
      if (oprint_time)
     &     call dft_tstamp(' calling densm ')
      call dft_densm(g_dens,g_movecs, 
     &               nbf_ao,nmo,ipol, 
     &               geom,AO_bas_han,noc,ntotel,
     &               Dbl_MB(k_eval(1)),dbl_mb(k_occ),
     &               ssmear,tdots,iter,.true.,
     .               fon, nel_fon,nmo_fon,ncore_fon,
     .               spinset.or.n_rabuck.ne.0,           
     &               rtdb)

      do i=1,nuc_num
         call rohf_dens(g_movecs_n(i),nbf_n,1,0, g_ndens(i),g_odens)
         call ga_dscal(g_ndens(i), 5.0d-1)
      enddo
c
      if(oprint_info) then
         pstrace=ga_ddot(g_dens,g_s)
         if(ga_nodeid().eq.0) write (luout,'(5x,a,1x,e15.7)')
     &      'tr(P*S): ',pstrace 
         if(ipol.eq.2) then
            pstrace=ga_ddot(g_dens(2),g_s)
            if(ga_nodeid().eq.0) write (luout,'(5x,a,1x,e15.7)')
     &           'tr(P_beta*S): ',pstrace 
         endif
      endif
c
      if (oprint_time)
     &     call dft_tstamp(' called densm ')
      

      etnew=etnew+tdots
c
      delta = -etold+etnew

      do ispin = 1, ipol
c     
c     
c     Check convergence on Density.
c     
         rms(ispin) = dft_dencvg(g_dens(ispin), g_denso(ispin), nbf_ao)
c     
c        Damping implemented here.
c     
         if (damping)then
            pp = dble(ndamp)*1.d-2
            call ga_dadd(pp, g_denso(ispin),(1d0-pp), g_dens(ispin), 
     &                   g_dens(ispin))
         else
            ndamp = 0
         endif
      enddo                     ! end loop over ispin
      do i=1,nuc_num
         rms(ipol+i) = dft_dencvg(g_ndens(i), g_ndenso(i), nbf_n)
      enddo
      call ga_sync
c     
      if (oprint_conv.and.iter.eq.1.and.me.eq.0)then
c
         nheap = MA_Inquire_Heap(MT_Dbl)
         nstack = MA_Inquire_Stack(MT_Dbl)
         write(LuOut,21)
         write(LuOut,'(10x,a,f10.2,i20)')
     &        ' Heap Space remaining (MW):  ',dble(nheap)*1.D-06,nheap
         write(LuOut,'(10x,a,f10.2,i20)')
     &        'Stack Space remaining (MW):  ',dble(nstack)*1.D-06,nstack
         call util_flush(LuOut)
         write(LuOut,1)
      endif

      if (oprint_conv.and.me.eq.0)then
         current_wall = util_wallsec()
         if (diising)then
            write(LuOut,2)ndamp,rlshift,
     &           iter, Etnew+Enuc,
     &           delta,sqrt(rms(1)),derr(1),current_wall
            if (ipol.eq.2)write(LuOut,3)sqrt(rms(2)),derr(2)
         else
            write(LuOut,22)ndamp,rlshift,
     &           iter, Etnew+Enuc,
     &           delta,sqrt(rms(1)), current_wall
            if (ipol.eq.2)write(LuOut,23)sqrt(rms(2))
         endif
         call util_flush(LuOut)
      endif
c
c     ecce ouput
c
      call ecce_print1 ('iteration counter', mt_int, iter, 1)
      call ecce_print1 ('iterative total energy difference', 
     &                  mt_dbl, delta, 1)
      call ecce_print1 ('iterative total density difference', 
     &                  mt_dbl, sqrt(rms(1)), 1)
c
      call ga_sync
c     
c     save eigenvectors to movecs file
c     
      if (.not.movecs_write(rtdb, ao_bas_han, movecs_out, 'dft', title,
     &                      nbf_ao, ipol, nmo, dbl_mb(k_occ), nbf_ao, 
     &                      dbl_mb(k_eval(1)), nbf_ao, g_movecs))
     &   call errquit('dft_scf: movec_write failed', 0, DISK_ERR)
! c
! c     print out the eigenvalues
! c
   !    call dft_prevals('inter',me,oprint_eval,oprint_vecs,
   !   ,     g_movecs,dbl_mb(k_eval(1)),dbl_mb(k_occ))

! c     
! c     If open shell compute overlap of alpha orbitals with beta 
! c     orbitals.
! c     
!       if ((ipol.gt.1).and.(oprint_interm_overlap)) then
!          call dft_mxspin_ovlp(nbf_ao,nmo,ao_bas_han, g_movecs(1), 
!      &       g_movecs(2),g_tmp)
!       endif
! c     
! c     computation of <S2> for open shell
! c     
!       if ((ipol.gt.1).and.(oprint_interm_S2)) then

!          call dft_s2_value(geom, AO_bas_han, .false., noc(1), noc(2),
!      &                     nbf_ao, g_dens(1), g_dens(2))
!       endif
! c     
c     
c     Form the total density matrix.
c     
      if (ipol.eq.2)then
         call ga_dadd(one, g_dens(1), one, g_dens(2), g_dens(1))
      endif
c     
c     Check for SCF convergence.
c     
      call ga_sync
      if (oprint_time)
     &     call dft_tstamp(' calling scfcvg ')

c
      call cneo_dft_scfcvg(rms, derr, Etold, Etnew,
     &            e_conv, d_conv, g_conv, ipol+nuc_num, 
     &            iter, iterations, idone, rtdb,
     &     converged, diising,
     B     back_to_atguess)

c
      if(movecs_in.ne.'atomic'.and.back_to_atguess.lt.0) then
         movecs_in='atomic'
         if(ga_nodeid().eq.0)
     W        write(luout,*) ' restart from atomic guess'
         if (.not. ga_destroy(g_denso(1))) call errquit
     &        ('dft_scf: could not destroy g_denso', 0, GA_ERR)
         if(ipol.eq.2) then
            if (.not. ga_destroy(g_denso(2))) call errquit
     &           ('dft_scf: could not destroy g_dens2', 0, GA_ERR)
         endif
         goto 1993
      endif
      if (oprint_time)
     &     call dft_tstamp(' called scfcvg ')
c
c     Check for remaining time to exit "gracefully"
c
      current_wall = util_wallsec()
      if ((iter-1).gt.1)then
         elapsed_wall = current_wall - save_wall
         save_wall = current_wall
      else
         elapsed_wall = current_wall - start_wall
         save_wall = current_wall
      endif
c
      if (converged)then
c
c        If converged probably need a few seconds to clean things up 
c        and calculate a few properties.
c
         wall_time_reqd = 5.0
      else
c
c        If not converged probably need at least the amount time
c        required for previous iteration (multiply by 1.2 to be on the safe side).
c
         wall_time_reqd = elapsed_wall*1.2d0
      endif
      int_wall_time_reqd = wall_time_reqd
      if (.not.util_test_time_remaining(rtdb, int_wall_time_reqd))then
         if (me.eq.0)then
            write(LuOut,*)
            call util_print_centered(LuOut,
     &           'Exiting due to time limitations.', 20, .true.)
            write(LuOut,*)
         endif
         goto 2000
      endif
      if (idone.eq.0.or.(iswitc.lt.2.and.iter.lt.iterations).or.
     .     (iterations.eq.0))
     &   go to 1000 ! begin new iteration
      if (idone.eq.1.and.(.not.last_time_energy))then
         last_time_energy = .true.
         ! if(rlshift.eq.0.and..not.lcdft) goto 2000
         go to 1000             ! build final total energies
      endif

c     
 2000 continue
c
      if (nExc.eq.1)then
c  jbecca START: adding in eelec here for field int. 
        Etnew = Ecore + Ecoul + Exc(1) + tdots + Edisp + eelec
     $                   + ecore_n - e_en - 0.5*e_nn
      else
        Etnew = Ecore + Ecoul + Exc(1) + Exc(2) + tdots + Edisp + eelec
     $                   + ecore_n - e_en - 0.5*e_nn
      endif
c
c
c    
      if (.not. ga_destroy(g_denso(1))) call errquit
     &   ('dft_scf: could not destroy g_denso', 0, GA_ERR)
      if(ipol.eq.2) then
          if (.not. ga_destroy(g_denso(2))) call errquit
     &     ('dft_scf: could not destroy g_denso', 0, GA_ERR)
      end if
      do i=1,nuc_num
         if (.not. ga_destroy(g_ndenso(i))) call errquit
     &      ('dft_scf: could not destroy g_denso', 0, GA_ERR)
      enddo
c
c
      call ga_sync()
c
c <-- MN solvation models
c
      if (me.eq.0.and.oprint)then
         if (.not.converged)then
            write(LuOut,*)
            call util_print_centered(LuOut,
     &           'Calculation failed to converge', 20, .true.)
            write(LuOut,*)
         endif
c
c        Tally up last energy 
c
c        DIM/QM JEM
        
         dft_time = dft_time+util_wallsec()
         if (nexc.le.1)then
          write(LuOut,222)etnew+enuc, 
     &                      ecore, 
     &                      ecoul, 
     &                      exc(1),
     &                      ecore_n,
     &                      e_en, 
     &                      enuc
         else
          write(LuOut,223)etnew+enuc, 
     &                      ecore, 
     &                      ecoul, 
     &                      exc(1), 
     &                      exc(2), 
     &                      ecore_n,
     &                      e_en, 
     &                      enuc
         end if
         write(luout,2222) rho_n
         write(luout,2223) dft_time
c
      endif
c
c MN solvation models -->
c
 1970 continue
c
c     do analysis of MO occupations per irrep
c
      call print_irrep_occp(ao_bas_han, nbf_ao, int_mb(k_ir),
     &                      dbl_mb(k_occ), ipol)
c
c     DIM/QM JEM
c     Final DIM dipoles
c
      if (ldimqm) then
        call dimqm_main(rtdb, ao_bas_han, geom, g_dens, .true.)
      end if

c
c     Print all the eigenvalues
c
      nprint = nmo(1)
      if (ga_nodeid() .eq. 0) then
        if (util_print('final evals', print_high)) then
         do ispin = 1, ipol
            if (ipol.eq.1) then
              call util_print_centered(6,'Final eigenvalues',
     &           20,.true.)
            else if (ipol.eq.2.and.ispin.eq.1) then
              call util_print_centered(6,'Final alpha eigenvalues',
     &           20,.true.)
            else
              call util_print_centered(6,'Final beta eigenvalues',
     &           20,.true.)
            end if ! ipol
            call output(dbl_mb(k_eval(ispin)),1,nprint,1,1,nmo,1,1)
            write(6,*)
            call util_flush(6)
         end do ! ispin
        end if ! final evals
      end if ! ga_nodeid()  
c
c     Vector analysis stolen from rohf.F
c
      do ispin = 1, ipol
         call movecs_fix_phase(g_movecs(ispin))
         if (util_print('final vectors analysis', print_default)) then
           ilo = 1
           ihi = nmo(ispin)
           if (ipol.eq.1) then
              blob='DFT Final Molecular Orbital Analysis'
           else if(ipol.eq.2.and.ispin.eq.1) then
              blob='DFT Final Alpha Molecular Orbital Analysis' 
           else
              blob='DFT Final Beta Molecular Orbital Analysis' 
           endif
           
           call movecs_print_anal(ao_bas_han, ilo, ihi, tanalyze, 
     &          g_movecs(ispin), 
     &          blob, 
     &          .true., dbl_mb(k_eval(ispin)), oadapt, 
     &          int_mb(k_ir+(ispin-1)*nbf_ao),
     &          .true., dbl_mb(k_occ+(ispin-1)*nbf_ao))
         endif
      enddo

      if (.not. rtdb_get(rtdb, 'dft:molden', mt_log, 1, molden))
     &   molden=.false.
      if(molden) call dft_nwmolden(rtdb,geom,
     I     ipol,dbl_mb(k_eval(1)),
     O     dbl_mb(k_occ),
     I     int_mb(k_ir),
     N     nbf_ao, nmo,
     I     g_movecs,scftype,'write')

      if (oprint_final_vecs) then
         if (me .eq. 0) then
            call util_print_centered(6,'Final MO vectors',40,.true.)
            write(6,*)
            call util_flush(6)
         end if
         call ga_sync()
         do ispin = 1, ipol
            call ga_print(g_movecs(ispin))
            call util_flush(6)
         end do
         do i = 1, nuc_num
            call ga_print(g_movecs_n(i))
            call util_flush(6)
         enddo
      end if
! c
c Jorge Garza/Begin
c
      if (condfukui.eq.1) then
         call fukui(g_movecs, k_eval, tol2e, rtdb, nExc, iVxc_opt, 
     &              g_xcinv, IOLGC, g_wght, g_xyz, g_nq, wght_GA,
     &              rho_n, irdens_atom,
     &              icetobfr, natoms)
      endif
c
c Jorge Garza/End
c
      if (mulliken) then
         call dft_mulwrap(me,g_dens,g_s)
c Qin / print Lowdin populations.
        if (.not. ga_duplicate(g_fock, g_tmp2, 'ga_temp'))
     &       call errquit('dft_scf: error creating ga',0,GA_ERR)
        if (me.eq.0) call dft_header
     &        (' Total Density - Lowdin Population Analysis')
        call diis_bld12(toll_s, svals, g_svecs,  g_tmp2, g_tmp, 3)
        call lowd_pop(geom, ao_bas_han, g_dens(1), g_tmp2, g_tmp)
        if (ipol.eq.2)then
         if (me.eq.0)call dft_header
     &        (' Spin Density - Lowdin Population Analysis')
         call ga_dadd(1d0,g_dens(1),-2.d0,g_dens(2),g_dens(2))
         call lowd_pop(geom,ao_bas_han,g_dens(2), g_tmp2, g_tmp)
         call ga_dadd(1d0,g_dens(1),-1.d0,g_dens(2),g_dens(2))
         call ga_dscal(g_dens(2),0.5d0)
        endif
        if (.not. ga_destroy(g_tmp2)) call errquit
     &    ('dft_scf: could not destroy g_tmp2', 0, GA_ERR)
      endif
c     
c     end infinite loop for SCF iterations
c     
c     Store energy and convergence status ... must store before
c     write movecs since date of insertion is used.
c     
      if (.not. rtdb_put(rtdb,'uhf:coulomb', mt_dbl, 1, Ecoul)) call
     $     errquit('uhf: writing ecoul failed', 0, RTDB_ERR)
c
      if (.not. rtdb_put(rtdb, 'dft:energy', MT_DBL, 1, (Etnew+Enuc)))
     &   call errquit('dft_scf: failed to store energy in rtdb', 0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:converged', MT_LOG, 1, converged))
     &   call errquit('dft_scf: failed to store converged in rtdb', 0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:alpha irreps', MT_INT, nbf_ao, 
     &   int_mb(k_ir)))
     &   call errquit('dft_scf: failed to MO irreps in rtdb', 0,
     &       RTDB_ERR)
      if (ipol.eq.2) then
      if (.not. rtdb_put(rtdb, 'dft:beta irreps', MT_INT, nbf_ao, 
     &   int_mb(k_ir+nbf_ao)))
     &   call errquit('dft_scf: failed to MO irreps in rtdb', 0,
     &       RTDB_ERR)
      endif
c
c     if present, force noscf=f to avoid having it in future restarts
      if (noscf) then 
        noscf = .false.
        if (.not.rtdb_put(rtdb, 'dft:noscf', mt_log, 1, noscf))
     &   call errquit('dft_scf: failed to store nscf rtdb', 0,
     &       RTDB_ERR)
      endif ! noscf
c
c     output energies and eigenvectors to disk
c     
      if (.not.movecs_write(rtdb, ao_bas_han, movecs_out, 'dft', title,
     &                      nbf_ao, ipol, nmo, dbl_mb(k_occ), nbf_ao, 
     &                      dbl_mb(k_eval(1)), nbf_ao, g_movecs))
     &                      call errquit('dft_scf: movec_write failed',
     &                      0, DISK_ERR)

      do i=1,nuc_num
         movecs_out_n = ' '
         write(str_i, '(I0)') i
         call util_file_name('movecs_n'//trim(str_i), .false.,
     $                                          .false.,movecs_out_n)
      !    print *, 'debug:movecs:',i,movecs_out_n
         if (.not.movecs_write(rtdb, bas_n(i), movecs_out_n, 'dft', 
     &                      title, nbf_n, 1, nmo_n, 
     &                      dbl_mb(k_occ_n(i)), nbf_n, 
     &                      dbl_mb(k_eval_n(i)), nbf_n, g_movecs_n(i)))
     &                      call errquit('dft_scf: movec_write failed',
     &                      0, DISK_ERR)
         if (.not. rtdb_cput(rtdb,'cneo:output vectors:n'//trim(str_i),
     $                      1, movecs_out_n))
     $       call errquit('dft_scf: CNEO nuc MO vectors not saved',0,
     &                     RTDB_ERR)
      enddo
c
c     Write out the alpha/beta split eigenvalues
c
      if (.not. rtdb_get(rtdb,'dft:save_evals',mt_log,1,save_evals))
     c     save_evals=.false.
      if(save_evals) then
         call dft_print_evals_ab(
     c        nbf_ao, dbl_mb(k_occ), dbl_mb(k_eval(1)))
         save_evals=.false.
c     unset save evals
         if (.not. rtdb_put(rtdb,'dft:save_evals',mt_log,1,save_evals))
     &        call errquit('dft_scf: rtdbput failed',1234,RTDB_ERR)
      endif
c     
c     Shut down DIIS.
c     
      if (icall(1).gt.0 .or. noscf)then
         call cneo_diis_driver(toll_s, derr,  ipol, 
     &           icall, nfock, nbf_e, nbf_n, geom, 
     &           ao_bas_han, g_focks, g_nfock,g_dens,
     &           g_ndens, g_tmp_n,
     &           g_tmp, g_svecs, svals, 
     &           g_svecs_n, svals_n, diising, 
     &           nodiis,int_mb(k_gdiis),lkeeps,
     .           g_s12m,g_s12p,.true.,
     $           nconstr,multipliers0,dbl_mb(k_mdiis))
      endif
      if(nconstr.ne.0) then
         if(k_hess.ne.-1) then
            if (.not. MA_free_heap(l_hess))
     .           call errquit(' dft_scf:cannot popstack',111, MA_ERR)
         endif
         if(k_mdiis.ne.-1) then
            if (.not. MA_free_heap(l_mdiis))
     .           call errquit(' dft_scf:cannot popstack',111, MA_ERR)
         endif
      endif
      if(k_gdiis.ne.-1) then
         if (.not. MA_free_heap(l_gdiis))
     .        call errquit(' dft_scf:cannot popstack',111, MA_ERR)
      endif
c     
c     If open shell compute overlap of alpha orbitals with beta orbitals.
c     
      if (ipol.gt.1 .and. util_print('final overlap',print_default))then
         call dft_mxspin_ovlp(nbf_ao,nmo,ao_bas_han, g_movecs(1), 
     &        g_movecs(2),g_tmp)
      endif
c
      if (wght_GA)then
         if (.not. ga_destroy(g_wght)) call errquit
     &      ('dft_scf: could not destroy g_wght', 0, GA_ERR)
         if (.not. ga_destroy(g_xyz)) call errquit
     &      ('dft_scf: could not destroy g_xyz', 0, GA_ERR)
         if (.not. ga_destroy(g_nq)) call errquit
     &      ('dft_scf: could not destroy g_nq', 0, GA_ERR)
      endif
c     
c     Restore alpha and beta densities.
c
      if (ipol .gt. 1)
     &   call ga_dadd(one,g_dens(1),onem,g_dens(2),g_dens(1))
c     
c     computation of <S2> for open shell
c     
      if (ipol.gt.1 .and. util_print('final s2', print_low))then

         call dft_s2_value(geom,AO_bas_han,.false.,noc(1),noc(2),
     &        nbf_ao,g_dens(1),g_dens(2))

      endif
c
c     write to rtdb that fractional occupancy exists
c
      if(fon.or.gotsmear(ssmear)) then
      if (.not.rtdb_put(rtdb, 'dft:fractional_occup', mt_log, 1,
     &   .true.))call errquit('dftscf: cannot write rtdb',1, DISK_ERR)
      endif
c     
c     computation of moments
c
      if (natoms .gt. 1) then
         if(oprintinertia)call geom_momint(geom)
         call dft_mpole(rtdb, ao_bas_han, ipol, g_dens(1), g_dens(2))
      endif
c     
c     print stolen for uhf.F
c     
      if (util_print('schwarz',print_high).and.(.not.CDFIT))then
         call schwarz_print(natoms, nshells_ao)
      endif
c     
      call dft_prevals('final',me,oprint_eval,oprint_vecs,
     ,     g_movecs,dbl_mb(k_eval(1)),dbl_mb(k_occ))
c
c     calculate virtual spectrum
c
      lvspec = .false.
      if (.not.rtdb_get(rtdb,'dft:lvspec',mt_log,1,lvspec)) 
     &         lvspec=.false.
      if (.not. rtdb_get(rtdb,'dft:owstart',mt_int,2,owstart)) 
     &         lvspec=.false.
      if (.not. rtdb_get(rtdb,'dft:owend',mt_int,2,owend)) 
     &         lvspec=.false.
      if (lvspec) then
       call util_file_name('vspec',.false.,.false.,vspecfilename)
       call dft_vspec(rtdb, ao_bas_han, vspecfilename, ipol, nmo, 
     &    nbf_ao, noc, g_movecs, k_eval, owstart, owend)
      endif
c
c     ECCE printout
c     
      call movecs_ecce(nbf_ao, nmo, 1, nmo(1), dbl_mb(k_eval(1)),
     &                 dbl_mb(k_occ), int_mb(k_ir), 
     &                 g_movecs(1), 'dft', 'alpha')
      if (ipol.eq.2)then ! spin-unrestricted
         call movecs_ecce(nbf_ao, nmo, 1, nmo(2), dbl_mb(k_eval(2)),
     &                    dbl_mb(k_occ+nbf_ao), int_mb(k_ir+nbf_ao), 
     &                    g_movecs(2), 'dft', 'beta')
      endif
      call ecce_print1 ('total energy', mt_dbl, (Etold+Enuc), 1)
      call ecce_print1 ('nuclear repulsion energy', mt_dbl, Enuc, 1)
      call ecce_print1 ('coulomb energy', mt_dbl, Ecoul, 1)
      call ecce_print1 ('exchange energy', mt_dbl, Exc(1), 1)
      if (nexc.gt. 1)then
         call ecce_print1 ('correlation energy', mt_dbl, Exc(2), 1)
      endif
c
      if (.not.ma_chop_stack(l_ir))
     &   call errquit('dft_scf: cannot chop stack',98, MA_ERR)
c     
      if (ipol.gt.1)then
         if (.not. ga_destroy(g_fockt)) call errquit
     &      ('dft_scf: could not destroy g_fockt', 0, GA_ERR)
      endif
c
      if(lkeeps) then
      if (.not. ga_destroy(g_s12p)) call errquit
     &   ('dft_scf: could not destroy g_s12p', 0, GA_ERR)
      if (.not. ga_destroy(g_s12m)) call errquit
     &   ('dft_scf: could not destroy g_s12m', 0, GA_ERR)
      endif
c
      call fock_2e_tidy(rtdb)
c     
      if (converged)then
         call ecce_print_module_exit('dft', 'ok')
      else
         call ecce_print_module_exit('dft', 'failed')
      endif
c     
c     eval deallocation moved here from inside iteration loop
c     
      if (.not.ma_chop_stack(lcntoce))
     &   call errquit('dft_scf: cannot chop stack',99, MA_ERR)
c
      cneo_dft_scf = converged
c
c !!! BGJ
      if (.not. rtdb_get(rtdb, 'bgj:poliz', mt_log,
     &     1, do_poliz)) then
         do_poliz = .false.
      endif
      if (do_poliz) then
         write(LuOut,*)'*** dft_scf: calling cphf_poliz'
         if (.not. cphf_poliz(rtdb)) ! Never executed.
     $        call errquit(' cphf_poliz: failed from dft_scf !',0,
     &       CALC_ERR)
      endif
c !!! BGJ
c

C
C     Print out matrices to stdout (useful for development and debugging)
C
      if (.not.rtdb_get(rtdb, 'dft:print_mats', mt_log, 1, lprint_mats))
     $     lprint_mats = .false.
      if (lprint_mats) call dft_print_mats (g_s, svals, g_svecs,
     $     nmo, g_dens, g_focks)

      
C
C     Save select matrices to file
C     
      if (.not.rtdb_get(rtdb, 'dft:save_mats', mt_log, 1, lsave_mats))
     $     lsave_mats = .false.
      if (lsave_mats) call dft_save_mats (g_dens)

      
      if (.not. ga_destroy(g_tmp)) call errquit
     &   ('dft_scf: could not destroy g_tmp', 0, GA_ERR)
      if (.not. ga_destroy(g_focks(1))) call errquit
     &        ('dft_scf: could not destroy g_focks1', 0, GA_ERR)
      if (ipol.gt.1) then
         if (.not. ga_destroy(g_focks(2))) call errquit
     &           ('dft_scf: could not destroy g_focks2', 0, GA_ERR)
      end if
c
c     == if frozen embedding ==
      if (frozemb) then
        if (.not. ga_destroy(g_frozemb)) call errquit
     &     ('dft_scf: could not destroy g_frozemb', 0, GA_ERR)
cc AJL/Begin/FDE
        if (frozemb_fde) 
c Destroy all allocations and reset variables that have been altered
     &  call dft_frozemb_fde_end(rtdb,g_vxc_fde,g_dens_fde, g_xcinv_fde)
      end if !frozemb
cc AJL/End
c
cdft...Deallocate cdft arrays
      if(lcdft) then
        do i = 1, nconstr
         if (.not. ga_destroy(g_constr(i))) call errquit
     &      ('dft_scf: could not destroy g_constr', 0, GA_ERR)
        enddo

        if (.not. ga_destroy(g_fockv(1))) call errquit
     &      ('dft_scf: could not destroy g_fockv', 0, GA_ERR)

        if(ipol.gt.1) then
         if (.not. ga_destroy(g_fockv(2))) call errquit
     &      ('dft_scf: could not destroy g_fockv', 0, GA_ERR)
        endif
      endif
c
c     Deallocate zora related arrays
      if (do_zora) then
       if (.not. ga_destroy(g_rel(1))) call errquit(
     &          'dft_scf: ga_destroy failed ',0, GA_ERR)
       if (.not. ga_destroy(g_zora_scale_sf(1))) call errquit(
     &          'dft_scf: ga_destroy failed ',0, GA_ERR)
       if(ipol.gt.1) then
         if (.not. ga_destroy(g_rel(2))) call errquit(
     &        'dft_scf: ga_destroy failed ',0, GA_ERR)
         if (.not. ga_destroy(g_zora_scale_sf(2))) call errquit(
     &        'dft_scf: ga_destroy failed ',0, GA_ERR)
       end if
      end if
c
c     deallocate x2c arrays if necessary
      if (do_x2c) then
         call x2c_scf_terminate(rtdb,ipol,g_rel,
     g     g_x2c_u,g_x2c_u2c,g_x2c_c2u)
      end if ! do_x2c
c
c     AOR begin
c     deallocate xdm arrays
      if (lxdm.eq.1) then
         call xc_xdm_cleanup(rtdb)
      endif
c     AOR end

      if (odftps) call pstat_off(ps_scfend)
      return
c     
 21   format(/,10x,' Memory utilization after 1st SCF pass: ')
    1 format(/,
     &     1x,'  convergence    iter        energy       DeltaE   ',
     &     'RMS-Dens  Diis-err    time'/
     &     1x,'---------------- ----- ----------------- --------- ',
     &     '--------- ---------  ------')
    2 format(1x,'d=',i2,',ls=',f3.1,',diis',1x,i5,f18.10,
     &     1p,3d10.2,0p,f8.1)
    3 format(51x,1p,2d10.2)
 22   format(1x,'d=',i2,',ls=',f3.1,6x,i5,f18.10,
     &     1p,2d10.2,10x,0p,f8.1)
 23   format(51x,1p,1d10.2)
 1111 format(15x,'Core Energy:              ',f22.10)
c
 222  format(//
     &     '       Total CNEO-DFT energy =', f22.12/
     &     '         One electron energy =', f22.12/
     &     '              Coulomb energy =', f22.12/
     &     '       Exchange-Corr. energy =', f22.12/
     &     '          Quantum nuc energy =', f22.12/
     &     'Quantum nuc-elec rep. energy =', f22.12/
     &     '    Nuclear repulsion energy =', f22.12/)
c
 223  format(//
     &     '       Total CNEO-DFT energy =', f22.12/
     &     '         One electron energy =', f22.12/
     &     '              Coulomb energy =', f22.12/
     &     '             Exchange energy =', f22.12/
     &     '          Correlation energy =', f22.12/
     &     '          Quantum nuc energy =', f22.12/
     &     'Quantum nuc-elec rep. energy =', f22.12/
     &     '    Nuclear repulsion energy =', f22.12/)
 224  format('    Dispersion correction =', f22.12/)
 225  format('             COSMO energy =', f22.12/)
c
 2221 format('       Scaling correction =', f22.12/)
 2222 format(' Numeric. integr. density =', f22.12/)
 2223 format('     Total iterative time =', f9.1,'s'//)
! c
 4224 format(10x,' HOMO = ',f7.3,' LUMO = ',f7.3,
     &              ' RLSHIFT = ',f7.3,' HL_GAP = ',f7.3)
c     
      end
