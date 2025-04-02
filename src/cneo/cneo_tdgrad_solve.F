      subroutine cneo_tdgrad_solve(rtdb,ihdl_geom,ihdl_bfao,ipol,
     +      nroots,nao,nocc,naoc,nav,nfc,nfv,nmo,g_mo,iptr_e_mo,
     +      tda,oskel,lhashf,otriplet,lhascd,tol2e,kfac,
     +      solve_thresh,solve_maxsub,solve_maxit,
     +      ihdl_bfn,nao_n,naoc_n,nav_n,nmo_n,g_mo_n,iptr_e_mo_n,
     +      g_r_e,g_r_n,g_z,g_z_n,tol2e_cneo)
      implicit none 
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "util.fh"
#include "ccneo_tdgrad_solve.fh"
#include "ccneo_constrain.fh"

      integer rtdb
      integer ihdl_geom
      integer ihdl_bfao,ipol
      integer nroots
      integer nao(2),nocc(2),naoc(2),nav(2),nfc(2),nfv(2),nmo(2)
      integer g_mo(2)
      integer iptr_e_mo
      logical tda,oskel,lhashf,otriplet,lhascd
      double precision tol2e,kfac
      double precision tol2e_cneo

      double precision solve_thresh ! convergence threshold for solver
      integer solve_maxit  ! the maximum number of iterations in solver
      integer solve_maxsub ! the maximum number expansion vectors
                           ! per root
      
      integer ihdl_bfn(nuc_num)
      integer nao_n(nuc_num)
      integer naoc_n(nuc_num)
      integer nav_n(nuc_num)
      integer nmo_n(nuc_num)
      integer g_mo_n(nuc_num)
      integer iptr_e_mo_n

      integer g_r_e(2)        ![input] global arrays for elec RHS
      integer g_r_n(nuc_num)  ![input] global arrays for nuc RHS
      integer g_z(2)          ![output] global arrays for Ze vectors
      integer g_z_n(nuc_num)  ![output] global arrays for Zn vectors

      external cneo_tdgrad_solve_product
      external cneo_tdgrad_solve_precond

      integer i
      character*255 pname

      integer idim,jdim   ! GA dimensions
      integer ichk,jchk   ! GA chunk sizes
      integer g_b         ! 2D GA for right-hand-sides
      integer g_x         ! 2D GA for solutions
      integer klo         ! location in vectors in 2D GA
      integer col         ! number of columns in 2D GA
      integer alo(3)      ! block lower limits in 3D GA
      integer ahi(3)      ! block upper limits in 3D GA
      integer ip          ! counter over spin components
      integer ir          ! counter over roots
      logical oroot
      logical oprint, oprint_debug

      integer g_x_e
      integer g_x_n(nuc_num)
      integer g_x_g(nuc_num)
      logical status

      pname = 'cneo_tdgrad_solve: '
      oprint=util_print('cneo_tdgrad_solve',print_default)
      oroot = ga_nodeid().eq.0

      call cneo_tdgrad_solve_init(rtdb,ihdl_geom,ihdl_bfao,ipol,
     +      nroots,nao,nocc,naoc,nav,nfc,nfv,nmo,g_mo,iptr_e_mo,
     +      tda,oskel,lhashf,otriplet,lhascd,tol2e,kfac,
     +      ihdl_bfn,nao_n,naoc_n,nav_n,nmo_n,g_mo_n,iptr_e_mo_n,
     +      g_int1e_r,tol2e_cneo)

      if (.not.ga_create(mt_dbl,ctdg_tot_size,nroots,'g_b',-1,-1,g_b))
     +   call errquit(pname//'failed to create g_b',0,GA_ERR)
      if (.not.ga_create(mt_dbl,ctdg_tot_size,nroots,"g_x",-1,-1,g_x)) 
     +  call errquit(pname//'failed to create g_x',0, GA_ERR)
      call ga_zero(g_b)

      if (.not. ga_create(mt_dbl,ctdg_e_size,ctdg_nroots,
     $         pname//'g_x_e',-1,-1,g_x_e)) call errquit
     $         (pname//'ga_create?',0,0)
      do i=1,nuc_num
         if (.not. ga_create(mt_dbl,ctdg_n_size(i),ctdg_nroots,
     $         pname//'g_x_n',-1,-1,g_x_n(i))) call errquit
     $         (pname//'ga_create?',0,0)

         if (.not. ga_create(mt_dbl,3,ctdg_nroots,
     $         pname//'g_x_g',-1,-1,g_x_g(i))) call errquit
     $         (pname//'ga_create?',0,0)
         call ga_zero(g_x_g(i))      
      enddo

      do ip = 1, ctdg_ipol
        do ir = 1, ctdg_nroots
          alo(1) = ir
          ahi(1) = ir
          alo(2) = 1
          ahi(2) = ctdg_naoc(ip)
          alo(3) = 1
          ahi(3) = ctdg_nav(ip)
          klo = (ip-1)*(ctdg_naoc(1)*ctdg_nav(1))+1
          col = ir
          call tddft_grad_mat_to_vec(g_r_e(ip),alo,ahi,g_b,klo,col,
     +                               -1.0d0,"=")
        enddo
      enddo

      klo = ctdg_e_size+1
      do i = 1, nuc_num
         do ir = 1, ctdg_nroots
            alo(1) = ir
            ahi(1) = ir
            alo(2) = 1
            ahi(2) = ctdg_naoc_n(i)
            alo(3) = 1
            ahi(3) = ctdg_nav_n(i)
            col = ir
            call tddft_grad_mat_to_vec(g_r_n(i),alo,ahi,g_b,klo,col,
     +                                 -1.0d0,'=')
         enddo
         klo = klo+ctdg_n_size(i)
      enddo


      call ga_copy(g_b,g_x)
      call cneo_tdgrad_solve_precond(g_x,0.0d0)

      call ga_lkain(rtdb,g_x,g_b,cneo_tdgrad_solve_product,
     +              cneo_tdgrad_solve_precond,solve_thresh,
     +              solve_maxsub*nroots,solve_maxit,.true.,oprint)

      call cneo_disassamble(ctdg_nroots,nuc_num,g_x,ctdg_e_size,
     $      ctdg_n_size,g_x_e,g_x_n,g_x_g)

      do ip = 1, ctdg_ipol
        do ir = 1, ctdg_nroots
          alo(1) = ir
          ahi(1) = ir
          alo(2) = 1
          ahi(2) = ctdg_naoc(ip)
          alo(3) = 1
          ahi(3) = ctdg_nav(ip)
          klo = (ip-1)*(ctdg_naoc(1)*ctdg_nav(1))+1
          col = ir
          call tddft_grad_vec_to_mat(g_z(ip),alo,ahi,g_x_e,klo,col)
        enddo
      enddo

      do i = 1, nuc_num
         do ir = 1, ctdg_nroots
            alo(1) = ir
            ahi(1) = ir
            alo(2) = 1
            ahi(2) = ctdg_naoc_n(i)
            alo(3) = 1
            ahi(3) = ctdg_nav_n(i)
            call tddft_grad_vec_to_mat(g_z_n(i),alo,ahi,g_x_n(i),1,ir)
         enddo
      enddo

      if (util_print('cneo tdgrad z',print_debug)) then
         if (oroot) write(LuOut,*)'DEBUG: '//pname//'Ze'
         call tddft_grad_print_array(ipol,nroots,g_z,dble(ipol))
         if (oroot) write(LuOut,*)'DEBUG: '//pname//'Zn'
         call tddft_grad_print_array(nuc_num,nroots,g_z_n,1.0d0)
       endif

      status = .true.
      status = status .and. ga_destroy(g_x)
      status = status .and. ga_destroy(g_b)
      status = status .and. ga_destroy(g_x_e)
      do i=1,nuc_num
         status = status .and. ga_destroy(g_x_n(i))
         status = status .and. ga_destroy(g_x_g(i))
      enddo
      if(.not. status) call errquit
     $            (pname//'ga_destroy failed',0,GA_ERR)
      call cneo_tdgrad_solve_tidy

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