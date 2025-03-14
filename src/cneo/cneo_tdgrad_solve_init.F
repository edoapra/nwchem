      subroutine cneo_tdgrad_solve_init(rtdb,geom,ihdl_bfao,ipol,
     +      nroots,nao,nocc,naoc,nav,nfc,nfv,nmo,g_mo,iptr_mo_e,
     +      tda,oskel,lhashf,otriplet,lhascd,tol2e,kfac,
     +      ihdl_bfn,nao_n,naoc_n,nav_n,nmo_n,g_mo_n,iptr_mo_e_n,
     +      g_int1e_r)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "ccneo_tdgrad_solve.fh"
#include "rtdb.fh"
! #include "cconstrain.fh"

      integer rtdb
      integer geom
      integer ihdl_bfao,ipol
      integer nroots
      integer nao,nocc(2),naoc(2),nav(2),nfc(2),nfv(2),nmo(2)
      integer g_mo(2)
      integer iptr_mo_e
      logical tda,oskel,lhashf,otriplet,lhascd
      double precision tol2e,kfac
      
      integer ihdl_bfn(nuc_num)
      integer nao_n(nuc_num)
      integer naoc_n(nuc_num)
      integer nav_n(nuc_num)
      integer nmo_n(nuc_num)
      integer g_mo_n(nuc_num)
      integer iptr_mo_e_n
      integer g_int1e_r(3,nuc_num)

      character*32 pname
      integer i

      pname= "cneo_tdgrad_solve_init: "

      ctdg_tol2e = tol2e
      ctdg_kfac = kfac

      ctdg_rtdb = rtdb
      ctdg_geom = geom
      ctdg_bfao = ihdl_bfao
      ctdg_ipol = ipol
      ctdg_nroots = nroots
      ctdg_nao = nao
      ctdg_nocc(1) = nocc(1)
      ctdg_nocc(2) = nocc(2)
      ctdg_naoc(1) = naoc(1)
      ctdg_naoc(2) = naoc(2)
      ctdg_nav(1)  = nav(1)
      ctdg_nav(2)  = nav(2)
      ctdg_nfc(1)  = nfc(1)
      ctdg_nfc(2)  = nfc(2)
      ctdg_nfv(1)  = nfv(1)
      ctdg_nfv(2)  = nfv(2)
      ctdg_nmo(1)  = nfc(1)+naoc(1)+nav(1)+nfv(1)
      ctdg_nmo(2)  = nfc(2)+naoc(2)+nav(2)+nfv(2)
      ctdg_g_mo(1) = g_mo(1)
      ctdg_g_mo(2) = g_mo(2)
      ctdg_e_mo    = iptr_mo_e
      ctdg_tda     = tda
      ctdg_oskel   = oskel
      ctdg_lhashf  = lhashf
      ctdg_triplet = .false.
      ctdg_lhascd  = lhascd

      if (otriplet) then
        if (.not.rtdb_put(rtdb,'fock_xc:triplet',mt_log,1,.false.))
     1    call errquit(pname//'failed to set triplet',0,RTDB_ERR)
      endif

      ctdg_tot_size = 0
      do i=1,nuc_num
         ctdg_bfn(i)    = ihdl_bfn(i)
         ctdg_nao_n(i)  = nao_n(i)
         ctdg_naoc_n(i) = naoc_n(i)
         ctdg_nav_n(i)  = nav_n(i)
         ctdg_nmo_n(i)  = nmo_n(i)
         ctdg_g_mo_n(i) = g_mo_n(i)
         ctdg_n_size(i) = naoc_n(i)*nav_n(i)
         ctdg_tot_size  = ctdg_tot_size+ctdg_n_size(i)
      enddo
      ctdg_e_mo_n = iptr_mo_e_n
      ctdg_e_size = naoc(1)*nav(1)+(ipol-1)*naoc(2)*nav(2)
      ctdg_tot_size = ctdg_tot_size+ctdg_e_size+3*nuc_num
      print *, 'tdgrad:debug:size info:',ctdg_e_size,ctdg_tot_size,
     +         ctdg_n_size(1),ctdg_n_size(2),nao,ipol

      call cneo_tdgrad_int1er_ov(nuc_num,g_int1e_r,g_mo_n,
     +            nao_n,naoc_n,nav_n,ctdg_g_int1er_ao,
     +            ctdg_g_int1er_ov)

      end

      subroutine cneo_tdgrad_solve_precond(g_x,shift)
c
c     This subroutine applies a preconditioner to the electronic part of
c     cneo_tdgrad_solve g_x array. Currently, the nuclear component and 
c     the constrained component remain unchanged.
c
c     This subroutine is modified from subroutine tddft_grad_solve_precond.
c     The modifications are:
c     1. decompose g_x into g_x_e,g_x_n(nuc_num),g_tdg(nuc_num)
c     2. apply tddft_grad_solve_precond to g_x_e
c     3. assemble updated g_x_e with old g_x_n and g_tdg to form the updated g_x
c
      implicit none 
#include "mafdecls.fh"
#include "global.fh"
#include "ccneo_tdgrad_solve.fh"
c
c     Input:
c
      double precision shift ! a level shift
c
c     Input/Output:
c
      integer g_x ! the new trial vectors
c
c     Local:
c
      integer ir  ! counter running over roots
      integer ioc ! counter over active occupied orbitals
      integer ivt ! counter over active virtual orbitals
      integer iov ! counter over occupied-virtual orbital pairs
c
      integer ilo,  ihi  ! limits on CPHF vector block
      integer irlo, irhi ! limits on roots block
      integer k_x        ! index for the local CPHF block
      integer ld_x       ! leading dimension for the local CPHF block
c
      logical oalpha ! .true. if there is an alpha orbital part
      logical obeta  ! .true. if there is an beta orbital part
c
      integer ialo, iahi ! limits on alpha CPHF vector block
      integer iblo, ibhi ! limits on beta CPHF vector block
c
      integer iav_start ! the first local alpha virtual orbital
      integer iao_start ! the first local alpha occupied orbital
      integer ibv_start ! the first local beta virtual orbital
      integer ibo_start ! the first local beta occupied orbital
c
      integer iproc     ! my processor rank
      logical doit

      integer g_x_e
      integer g_x_n(nuc_num)
      integer g_tdg(nuc_num)
      integer i
      logical status
      character*255 pname

      pname = 'cneo_tdgrad_precond: '
c
c     1. create local arrays and disassemble g_x
      if (.not. ga_create(mt_dbl,ctdg_e_size,ctdg_nroots,
     $         pname//'g_x_e',32,32,g_x_e)) call errquit
     $         (pname//'ga_create?',0,0)
      do i=1,nuc_num
         if (.not. ga_create(mt_dbl,ctdg_n_size(i),ctdg_nroots,
     $         pname//'g_x_n',32,32,g_x_n(i))) call errquit
     $         (pname//'ga_create?',0,0)

         if (.not. ga_create(mt_dbl,3,ctdg_nroots,
     $         pname//'g_tdg',32,32,g_tdg(i))) call errquit
     $         (pname//'ga_create?',0,0)
         
      enddo

      call cneo_disassamble(ctdg_nroots,nuc_num,g_x,ctdg_e_size,
     $         ctdg_n_size,g_x_e,g_x_n,g_tdg)

c     2. (copied) apply tddft_grad_solve_precond to g_x_e
      call ga_sync()
      iproc = ga_nodeid()
      call ga_distribution(g_x_e,iproc,ilo,ihi,irlo,irhi)
      doit = .not. (ilo.eq.0.and.ihi.eq.-1.and.
     &     irlo.eq.0.and.irhi.eq.-1)
      if (doit) then
       call ga_access(g_x_e,ilo,ihi,irlo,irhi,k_x,ld_x)
c
c      Do we have alpha or beta orbitals or both?
c
       oalpha = ilo.le.ctdg_naoc(1)*ctdg_nav(1)
       obeta  = ihi.gt.ctdg_naoc(1)*ctdg_nav(1)
c
c      If we have alpha orbitals then work on those elements
c
       ialo =  0
       iahi = -1
       if (oalpha) then
c
c       First we work out the limits on the loops
c
        ialo      = ilo
        iahi      = min(ctdg_naoc(1)*ctdg_nav(1),ihi)
        iav_start = (ialo+ctdg_naoc(1)-1)/ctdg_naoc(1)
        iao_start = mod(ialo-1,ctdg_naoc(1))+1
c
c       Loop over the roots
c
        do ir = irlo, irhi
          ioc = iao_start
          ivt = iav_start
          iov = 0
c
c         Loop over occupied (fastest) and virtual orbitals
c
          do while (iov.le.iahi-ialo)
            do while (iov.le.iahi-ialo.and.ioc.le.ctdg_naoc(1))
              dbl_mb(k_x+(ir-1)*ld_x+iov) 
     &        = dbl_mb(k_x+(ir-1)*ld_x+iov)/
     &         (dbl_mb(ctdg_e_mo+ctdg_nocc(1)+ivt-1)-
     &          dbl_mb(ctdg_e_mo+ctdg_nfc(1)+ioc-1))
              ioc = ioc + 1
              iov = iov + 1
            enddo
            ioc = 1
            ivt = ivt + 1
          enddo
        enddo
      endif
      if (obeta) then
        iblo      = max(ilo-ctdg_naoc(1)*ctdg_nav(1),1)
        ibhi      = ihi-ctdg_naoc(1)*ctdg_nav(1)
        ibv_start = (iblo+ctdg_naoc(2)-1)/ctdg_naoc(2)
        ibo_start = mod(iblo-1,ctdg_naoc(2))+1
        do ir = irlo, irhi
          ioc = ibo_start
          ivt = ibv_start
          iov = iahi-ialo+1  ! the first after the alpha spin elements
          do while (iov.le.iahi-ialo+1+ibhi-iblo)
            do while (iov.le.iahi-ialo+1+ibhi-iblo.and.
     &                ioc.le.ctdg_naoc(2))
              dbl_mb(k_x+(ir-1)*ld_x+iov) 
     &        = dbl_mb(k_x+(ir-1)*ld_x+iov)/
     &         (dbl_mb(ctdg_e_mo+ctdg_nmo(1)+ctdg_nocc(2)+ivt-1)-
     &          dbl_mb(ctdg_e_mo+ctdg_nmo(1)+ctdg_nfc(2)+ioc-1))
              ioc = ioc + 1
              iov = iov + 1
            enddo
            ioc = 1
            ivt = ivt + 1
          enddo
        enddo
      endif
      call ga_release(g_x_e,ilo,ihi,irlo,irhi)
      endif ! doit
c
c     3. assemble updated g_x_e with old g_x_n and g_tdg
c
      call ga_zero(g_x)
      call cneo_assamble(ctdg_nroots,nuc_num,g_x_e,g_x_n,g_tdg,
     +      ctdg_e_size,ctdg_n_size,g_x)

      status = .true.
      status = status .and. ga_destroy(g_x_e)
      do i=1,nuc_num
         status = status .and. ga_destroy(g_x_n(i))
         status = status .and. ga_destroy(g_tdg(i))
      enddo
      if (.not. status) call errquit(pname//'ga_destroy?',0,0)
c
      call ga_sync
      return
      end      
