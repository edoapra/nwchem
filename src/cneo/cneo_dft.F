      logical function cneo_dft_energy(rtdb)
      implicit none 
#include "errquit.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "ccneo.fh"
#include "geom.fh"
#include "bas.fh"
      integer rtdb
      logical status
      logical cneo_nwdft
      external cneo_nwdft

      integer geom, geom_cneo, basis
      integer bases(max_quantum_nuc+1)
      integer bases_nuc(max_quantum_nuc)
      integer i
      logical int_normalize
      external int_normalize

      cneo_dft_energy = .false.

      status = geom_create(geom,'geometry')
      status = geom_rtdb_load(rtdb,geom,'geometry')
      
      call cneo_nuc_init(rtdb,geom)
      call cneo_geom_init(rtdb,geom_cneo)

      status = bas_create(basis,'ao basis')
      status = bas_rtdb_load(rtdb,geom_cneo,basis,'ao basis')
      status = int_normalize(rtdb,basis)
      bases(1) = basis
      call cneo_bas_init(rtdb,geom_cneo,bases_nuc)
      do i=1,nuc_num
         status = int_normalize(rtdb,bases_nuc(i))
         bases(i+1) = bases_nuc(i)
      enddo

      call int_init(rtdb,nuc_num+1,bases(1:nuc_num+1))
      call cneo_int1er_init(rtdb,geom_cneo,bases_nuc)
      call cneo_dft_init(rtdb,geom_cneo,basis,bases_nuc)
      call cneo_schwarz_init(geom_cneo,nuc_num+1,bases)

      if(.not. rtdb_cput(rtdb,'dft:theory', 1, 'dft')) call errquit
     +      ('cneo_dft_energy: could not put dft:theory',0,0)
      cneo_dft_energy = cneo_nwdft(rtdb)
      call cneo_int1er_tidy
      call cneo_schwarz_tidy(nuc_num+1,bases)
      call cneo_dft_tidy
      call grid_cleanup(.true.)
      if (.not. geom_destroy(geom)) call errquit
     $   ('cneo_dft:geom_destroy?',0,0)

      end

      logical function cneo_nwdft(rtdb)
c     
C$Id$
c
      implicit none
c     
c     == Include standard input and memory checking common structures ==
c     
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "cdft.fh"
#include "util.fh"
#include "cgridfile.fh"
#include "cosmo.fh"
#include "modelpotential.fh"
#include "zora.fh"
#include "ccneo.fh"
#include "ccneo_dft.fh"
cinclude "crohf.fh"
c
      integer rtdb
c     
c     == Local declarations ==
      double precision energy
      integer iproc,itype
cc AJL/Begin/FDE
c      integer bases(3), i3
      integer bases(nw_max_basis_sets), i3
cc AJL/End
      logical LResult, oprint, converged
      logical dyall_mod_dir  ! do modified Dirac calculation
c
      logical cneo_dft_main0d,movecs_converged,grid_reopen,xc_gotxc  
      external cneo_dft_main0d,movecs_converged,grid_reopen,xc_gotxc  
      logical grid_ok,l1ecache
      integer igok
      character*80 theory
      logical lrt_tddft
      character*32 pname
      integer nclosed, nopen
c
      integer natoms
      logical dftmp2
c
cc AJL/Begin/FDE
c      logical frozemb_fde
cc AJL/End

      logical rt_tddft
      external rt_tddft
      logical xc_hascam,int_forcenwints
      external xc_hascam,int_forcenwints
      character*24 rtdb_errmsg

      integer geom0
      integer i
      double precision tmp_charge

      logical status
c     
      pname = "cneo_dft:"
      cneo_nwdft = .false.
      iproc  = ga_nodeid()
c
c     == Make sure database is parallel, push context down to DFT ==
      LResult = rtdb_parallel(.true.)
      call util_print_push
      call util_print_rtdb_load(rtdb, 'dft')
c
c     == Init pstat ==
      call dft_pstat_init(rtdb)
c
c     == Init fdist ==
      call dft_fdist_init(rtdb)
c
c     == If Dyall-mod-dir is set, exit with error ==
      if (rtdb_get(rtdb,'dyall_mod_dir',mt_log,1,dyall_mod_dir)) then
        if (dyall_mod_dir) call errquit(pname//
     &      'Cannot do DFT calculation with Dyall-modified-Dirac method'
     &      ,99, RTDB_ERR)
      end if
c
c     == Print header ==
      oprint = util_print('information', print_low)
      lrt_tddft = .false.
      dftmp2 = .false.
c
      if (oprint) then
        if (.not. rtdb_cget(rtdb, 'title', 1, title))
     &       title = ' '
        if(iproc.eq.0) then
         write(LuOut,*)
         call util_print_centered(
     &         LuOut, 'NWChem CNEO-DFT Module', 40, .true.)
          write(LuOut,*)
          write(LuOut,*)
          if (title .ne. ' ') then
            call util_print_centered(LuOut, title, 40, .false.)
            write(LuOut,*)
            write(LuOut,*)
          endif
        endif
      endif
      if (.not. rtdb_put(rtdb,'dft:skeleton',mt_log,1,.false.))
     +      call errquit(pname//'block symmetry',0,0)
      if (.not. rtdb_put(rtdb,'dft:adapt',mt_log,1,.false.))
     +      call errquit(pname//'block adapt',0,0)
c
c     == Gather input from all available resources e.g. input unit, ==
c     == old checkpoint files, etc. ==
      call dft_rdinput(rtdb)
      tol2e_cneo = min(1.0d-8,g_conv*1d-2)
      if (.not. bas_destroy(AO_bas_han))
     $      call errquit('cneo_dft:bas_destroy?',0,0)

      AO_bas_han = bas_e
c
      call dft_printinfo()
c
      cosmo_last = .false.
      do_zora = .false.
c
      use_modelpotential = .false.
c
c     == Check to see if system is already converged? ==
      if (movecs_in .eq. movecs_out) then
         if (.not. rtdb_cget(rtdb, 'dft:theory', 1, theory))
     $        call errquit(pname//'Theory not specified',0, RTDB_ERR)
c
c
         if(theory .eq. 'dft')then 
            if(.not.xc_gotxc()) then
               grid_ok=.true.
            else
               grid_ok=grid_reopen(geom)
            endif
            igok=0
            if(grid_ok) igok=1
            call ga_igop(999,igok,1,'+')
            grid_ok=igok.eq.ga_nnodes()
            if(.not.grid_ok) call grid_cleanup(.true.)
            if(.not.dftmp2) then
            if (movecs_converged(rtdb, ao_bas_han,theory, movecs_in).
     &      and.grid_ok)     then
c
c     check if we can get the grid from a file
c         
               if(.not.rtdb_get(rtdb,'cneo:energy', mt_dbl,
     &             1,energy)) call errquit(
     &             pname//'Failed getting converged energy',0,RTDB_ERR)
               if (iproc.eq.0 .and. oprint)then
                  write(LuOut,1011)energy
                  call util_flush(Luout)
               endif
               converged = .true.
               cneo_nwdft = .true.
               if (.not. lrt_tddft) goto 1101
            endif
         endif
      endif
      endif
 1011 format(/'  The CNEO-DFT is already converged '//,
     &     '         Total DFT energy =', f20.12/)
c     
c     == Initialize the integral object ==
c 
cc AJL/Begin/FDE
c
      i3 = 1
      bases(1) = AO_bas_han
      do i=1,nuc_num
         bases(i3+i) = bas_n(i)
      enddo
      i3 = i3+nuc_num
c
cc AJL/End
c
      if(xc_hascam(rtdb)) then
c
c     Force use of nw ints for CAM calculations
         if(.not.int_forcenwints(rtdb,rtdb_errmsg))
     &        call errquit(rtdb_errmsg,0, RTDB_ERR)
      endif
      
c
      call print_integrals(bases, oskel) ! print integrals, if needed
c
c     == Get 1-el in cache ==
      call cneo_schwarz_set(ao_bas_han)
c
      l1ecache=.false.
c
c     == Build the grid ==
      grid_written=.false.
      call grid_setquad(rtdb)
c
c     == Analyze set parameters for consistency ==
      geom0 = geom
      geom = geom_cneo

      call dft_inpana(rtdb)
c
c     == Pump the settings into cscf.fh ==
      call dft_cscf_init(rtdb,g_conv,rcharge,oskel,oadapt,iterations,
     &                   geom,ao_bas_han,mult-1,nbf_ao,noc,
     &                   movecs_in,movecs_out,movecs_guess,title)

c
c     == Set stuff in crohf.fh ==
c     if (mult.eq.1) then
c       nclosed = noc(1)
c       nopen   = 0
c     else
c       nclosed = min(noc(1),noc(2))
c       nopen   = mult - 1
c     endif
c     crohf_vlen = (nbf_ao-nclosed-nopen)*(nclosed+nopen)
c    &           + (nclosed*nopen)
c     crohf_init_flag = 1
c     
c     == Determine system type (molecule, polymer,surface,crystal) ==
c     == if no system type exists in db --> default = molecular case ==    
      if(.not. geom_systype_get(geom,itype)) itype = 0

c
      if(itype.eq.0) then
#ifdef GA_TRACE
        call trace_init(100000)                       ! initialize trace
#endif
         cneo_nwdft = cneo_dft_main0d(rtdb)

#ifdef GA_TRACE
        call trace_end(ga_nodeid())                   ! end trace
#endif
      else
        call errquit(pname//'Problem with system type',0, INPUT_ERR)
      endif

c     
c     == Terminate integrals ==
      if(l1ecache) call int_1e_uncache_ga()
      ! call schwarz_tidy
      call int_terminate
c
c     == Go here if already converged ==
c
 1101 continue
c     
c     == Destroy basis set handles and geom object ==
      XCFIT    = .FALSE.
      CDFIT    = .FALSE.
c
cc AJL/BEGIN
c      frozemb_fde  = .true.
c      if( FD_bas_han.eq.-99999 ) frozemb_fde = .false.

cc AJL/END
c
      status = .true.
      status = status .and. bas_destroy(ao_bas_han)
      status = status .and. geom_destroy(geom)
      status = status .and. geom_destroy(geom0)
      do i=1,nuc_num
         status = status .and. bas_destroy(bas_n(i))
      enddo
      if (.not. status) call errquit
     $      ('cneo_nwdft:destroying geom and basis handles')

c     == Restore charge ==
      if (.not. rtdb_get(rtdb,'charge_real',MT_DBL,1,tmp_charge))
     $   call errquit('cneo_nwdft:read tot charge',0,0)
      if (.not. rtdb_put(rtdb,'charge',MT_DBL,1,tmp_charge))
     $   call errquit('cneo_nwdft:restore tot charge',0,0)
c
c     == Flush output buffer ==
      if(iproc.eq.0) call util_flush(LuOut)
c
      call dft_pstat_print
      call util_print_pop
c
      return
      end