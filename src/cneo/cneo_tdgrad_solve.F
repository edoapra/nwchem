      subroutine cneo_tdgrad_solve(rtdb,ihdl_geom,ihdl_bfao,ipol,
     +      nroots,nao,nocc,naoc,nav,nfc,nfv,nmo,g_mo,iptr_e_mo,
     +      tda,oskel,lhashf,otriplet,lhascd,tol2e,kfac,
     +      ihdl_bfn,nao_n,naoc_n,nav_n,nmo_n,g_mo_n,iptr_e_mo_n,
     +      g_r_e,g_r_n,g_z,g_z_n)
      implicit none 
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "ccneo_tdgrad_solve.fh"
#include "cconstrain.fh"

      integer rtdb
      integer ihdl_geom
      integer ihdl_bfao,ipol
      integer nroots
      integer nao(2),nocc(2),naoc(2),nav(2),nfc(2),nfv(2),nmo(2)
      integer g_mo(2)
      integer iptr_e_mo
      logical tda,oskel,lhashf,otriplet,lhascd
      double precision tol2e,kfac
      
      integer ihdl_bfn(nuc_num)
      integer nao_n(nuc_num)
      integer naoc_n(nuc_num)
      integer nav_n(nuc_num)
      integer nmo_n(nuc_num)
      integer g_mo_n(nuc_num)
      integer iptr_e_mo_n

      integer g_r_e(2)
      integer g_r_n(nuc_num)
      integer g_z(2)
      integer g_z_n(nuc_num)

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

      integer g_x_e
      integer g_x_n(nuc_num)
      integer g_x_g(nuc_num)
      logical status

      pname = 'cneo_tdgrad_solve: '

      call cneo_tdgrad_solve_init(rtdb,ihdl_geom,ihdl_bfao,ipol,
     +      nroots,nao,nocc,naoc,nav,nfc,nfv,nmo,g_mo,iptr_e_mo,
     +      tda,oskel,lhashf,otriplet,lhascd,tol2e,kfac,
     +      ihdl_bfn,nao_n,naoc_n,nav_n,nmo_n,g_mo_n,iptr_e_mo_n,
     +      g_int1e_r)

      if (.not.ga_create(mt_dbl,ctdg_tot_size,nroots,'g_b',-1,-1,g_b))
     +   call errquit(pname//'failed to create g_b',0,GA_ERR)
      if (.not.ga_create(mt_dbl,ctdg_tot_size,nroots,"g_x",-1,-1,g_x)) 
     +  call errquit(pname//'failed to create g_x',0, GA_ERR)
      call ga_zero(g_b)

      ! call ga_fill(g_x,1.0d0)
      ! call cneo_tdgrad_solve_product(0.0d0,g_x,g_b)
      ! print *,'tdgrad:debug:product:test:'
      ! call ga_print(g_b)
      ! call ga_print(g_x)
      ! call ga_zero(g_b)
      ! return
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
      print *, 'tdgrad:debug:g_b:'
      call ga_print(g_b)
      call cneo_tdgrad_solve_precond(g_x,0.0d0)

      call ga_lkain(rtdb,g_x,g_b,cneo_tdgrad_solve_product,
     +              cneo_tdgrad_solve_precond,1.0d-8,
     +              250,40,.true.,.true.)

      print *, 'tdgrad:debug:solve g_x:'
      call ga_print(g_x)

   !    do i=1,nuc_num
   !       if (.not. ga_create(mt_dbl,3,ctdg_nroots,
   !   $         pname//'g_tdg_tmp',-1,-1,g_x_g(i))) call errquit
   !   $         (pname//'ga_create?',0,0)
   !    enddo

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

      print *, 'tdgrad:debug:g_z'
      call ga_print(g_z(1))
      call ga_print(g_z_n(1))
      ! call ga_print(g_z_n(2))

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

      end