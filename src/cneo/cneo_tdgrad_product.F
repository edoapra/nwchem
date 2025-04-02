subroutine cneo_tdgrad_solve_product(acc,g_x,g_Ax)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "ccneo_tdgrad_solve.fh"
      double precision acc
      integer g_x
      integer g_Ax

      integer i, ip, ir, j
      character*255 pname

      integer g_x_e
      integer g_x_n(nuc_num)
      integer g_tdg(nuc_num)

      integer g_x_3d_mo(2)
      integer g_x_3d_ao
      integer g_Ax_3d_ao
      integer g_Ax_3d_mo(2)
      integer g_x_3d_n(nuc_num)
      integer g_Ax_3d_n(nuc_num)
      integer g_Ax_n_mo(nuc_num)
      integer g_ambt
      ! integer g_t

      integer idim(3)
      integer ichnk(3)
      integer alo(3)
      integer ahi(3)
      integer blo(2)
      integer bhi(2)
      integer clo(3)
      integer chi(3)
      integer klo
      integer col
      integer iwhich

      integer offset_e  ! offset for nuclear orbital energy pointer

      integer g_Ax_e
      integer g_Ax_n(nuc_num)
      integer g_Ax_g(nuc_num)
      double precision tmp_dbl

      logical status

      pname = 'cneo_tdgrad_product: '
      iwhich =0
      if (.not. ga_create(mt_dbl,ctdg_e_size,ctdg_nroots,
     $         pname//'g_x_e',32,32,g_x_e)) call errquit
     $         (pname//'ga_create?',0,0)
      if (.not. ga_duplicate(g_x_e,g_Ax_e,pname//'g_Ax_e'))
     $      call errquit(pname//'create g_Ax_e',0,0)
      do i=1,nuc_num
         if (.not. ga_create(mt_dbl,ctdg_n_size(i),ctdg_nroots,
     $         pname//'g_x_n',32,32,g_x_n(i))) call errquit
     $         (pname//'ga_create?',0,0)
         if (.not. ga_duplicate(g_x_n(i),g_Ax_n(i),pname//'g_Ax_n'))
     $      call errquit(pname//'create g_Ax_n',0,0)

         if (.not. ga_create(mt_dbl,3,ctdg_nroots,
     $         pname//'g_tdg',32,32,g_tdg(i))) call errquit
     $         (pname//'ga_create?',0,0)
         if (.not. ga_duplicate(g_tdg(i),g_Ax_g(i),pname//'g_Ax_g'))
     $      call errquit(pname//'create g_Ax_g',0,0)
      enddo

      call cneo_disassamble(ctdg_nroots,nuc_num,g_x,ctdg_e_size,
     $         ctdg_n_size,g_x_e,g_x_n,g_tdg)

c  ========= Start elec part:H^{+}[Ze]+orbital energy+H_en[Zn]=======
c     The first part of electronic part is copied from 
c     nwdft/lr_tddft_grad/tddft_grad_solve_product.F
c     In principle, this can be done by calling subroutine 
c     tddft_grad_solve_product with g_x_e as input.
c     But the common file tddft_grad_solve_common.fh is a 
c     local common file and it's not straightforward to call
c
c     1. g_x 2D MO -> g_x_3d_mo (both e and n)
c
      call tddft_grad_create_r(ctdg_ipol,ctdg_naoc,ctdg_nav,ctdg_nroots,
     +                         g_x_3d_mo)
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
          call tddft_grad_vec_to_mat(g_x_3d_mo(ip),alo,ahi,
     $               g_x_e,klo,col)
        enddo
      enddo

      call cneo_tdgrad_create_r_n(nuc_num,ctdg_naoc_n,ctdg_nav_n,
     $         ctdg_nroots,g_x_3D_n)
      do ir=1, ctdg_nroots
         do i=1,nuc_num
            alo(1)=ir
            ahi(1)=ir
            alo(2)=1
            ahi(2)=ctdg_naoc_n(i)
            alo(3)=1
            ahi(3)=ctdg_nav_n(i)
            call tddft_grad_vec_to_mat(g_x_3D_n(i),alo,ahi,
     $                  g_x_n(i),1,ir)
         enddo
      enddo
c
c     2. g_x_3d_mo -> g_x_3d_ao
c
      idim(1) = ctdg_nroots*ctdg_ipol
      idim(2) = ctdg_nao
      idim(3) = ctdg_nao
      ichnk(1) = ctdg_nroots*ctdg_ipol
      ichnk(2) = -1
      ichnk(3) = -1
      if (.not.nga_create(mt_dbl,3,idim,'g_x_3d_ao',ichnk,g_x_3d_ao))
     +    call errquit(pname//'failed to create g_x_3d_ao',0,GA_ERR)
      call ga_zero(g_x_3d_ao)
      call tddft_grad_trans_mo2ao(ctdg_ipol,ctdg_nao,ctdg_nfc,ctdg_naoc,
     +     ctdg_nocc,ctdg_nav,ctdg_nfv,ctdg_nroots,
     +     1.0d0,0.0d0,"ib",ctdg_g_mo,g_x_3d_mo,"ib",g_x_3d_ao)
c
c     3. compute (A+B)g_x currently we have to compute (A-B)g_x as well
c        although we do not want it.
c
      if (.not.nga_create(mt_dbl,3,idim,'vec (A+B)X',ichnk,g_Ax_3d_ao))
     +    call errquit(pname//'failed to create g_Ax_3d_ao',0,GA_ERR)
c
      if (.not.ctdg_tda) then
        if (.not.nga_create(mt_dbl,3,idim,'vec (A-B)X',ichnk,g_ambt))
     +      call errquit(pname//'failed to create g_ambt',0,GA_ERR)
      endif
      call ga_zero(g_Ax_3d_ao)
      call tddft_nga_cont(ctdg_rtdb,ctdg_geom,ctdg_bfao,g_x_3d_ao,
     +     g_Ax_3d_ao,g_ambt,ctdg_nao,ctdg_ipol,ctdg_tol2e,ctdg_tda,
     +     ctdg_oskel,ctdg_kfac,ctdg_lhashf,ctdg_triplet,ctdg_nroots,
     +     iwhich,ctdg_lhascd)
      if (.not.ctdg_tda) then
        if (.not.ga_destroy(g_ambt))
     +      call errquit(pname//'failed to destroy g_ambt',0,GA_ERR)
      endif
c
c     4. Destroy g_x_3d_ao
c
      if (.not.ga_destroy(g_x_3d_ao)) 
     +   call errquit(pname//'failed to destroy g_x_3d_ao',0,GA_ERR)
c
c     5. Transform g_Ax_3d_ao -> g_Ax_3d_mo
c
      call tddft_grad_create_r(ctdg_ipol,ctdg_naoc,ctdg_nav,ctdg_nroots,
     +     g_Ax_3d_mo)
      ! call ga_zero(g_Ax_3d_mo)
      call tddft_grad_trans_ao2mo(ctdg_ipol,ctdg_nao,ctdg_nfc,ctdg_naoc,
     +     ctdg_nocc,ctdg_nav,ctdg_nfv,ctdg_nroots,1.0d0,0.0d0,"ib",
     +     ctdg_g_mo,g_Ax_3d_ao,g_Ax_3d_mo,"ib")
c
c     6. Destroy g_Ax_3d_ao
c
      if (.not.ga_destroy(g_Ax_3d_ao)) 
     +   call errquit(pname//'failed to destroy g_Ax_3d_ao',0,GA_ERR)
c
c     7. Add orbital energy difference term
c
      call tddft_grad_epsilon(g_x_3d_mo,g_Ax_3d_mo,ctdg_ipol,
     +     ctdg_nroots,ctdg_nao,ctdg_nocc,ctdg_nfc,dbl_mb(ctdg_e_mo))
c     
c     8. add en term: H_en[Z_n]
c  
      do ip=1,ctdg_ipol
         do i=1,nuc_num
            call cneo_tdgrad_response_ov(i,0,'enints',ctdg_bfao,
     $            ctdg_bfn(i),-4.0d0,g_x_3d_n(i),ctdg_g_mo(ip),   ! check later
     $            ctdg_g_mo_n(i),ctdg_nao,ctdg_nao_n(i),ctdg_naoc(ip),
     $            ctdg_naoc_n(i),ctdg_nav(ip),ctdg_nav_n(i),
     $            ctdg_nfc(ip),0,ctdg_nfv(ip),0,ctdg_nroots,
     $            g_Ax_3d_mo(ip),ctdg_geom,
     $            ctdg_nocc(ip),ctdg_naoc_n(i),'ib',ctdg_tol2e_cneo)
         enddo
      enddo
c
c     9. Transform g_Ax_3d_mo -> g_Ax 2D MO
c
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
          call tddft_grad_mat_to_vec(g_Ax_3d_mo(ip),alo,ahi,g_Ax_e,
     +                               klo,col,1.0d0,"=")
        enddo
      enddo

c  ========= Start nuc part:orbital energy+H_ne[Ze]+H_nn[Zn]=======
c
c     10. calculate g_Ax_n: add orbital contribution
      call cneo_tdgrad_create_r_n(nuc_num,ctdg_naoc_n,ctdg_nav_n,
     $         ctdg_nroots,g_Ax_n_mo)
      offset_e = 0
      call ga_zero(g_Ax_n_mo(1))
      do i=1,nuc_num
         call tddft_grad_epsilon(g_x_3d_n(i),g_Ax_n_mo(i),1,
     $            ctdg_nroots,ctdg_nao_n(i),ctdg_naoc_n(i),0,
     $            dbl_mb(ctdg_e_mo_n+offset_e))
         offset_e = offset_e+ctdg_nao_n(i)
      enddo
c     
c     11. add ne and nn response to g_Ax_n_mo
c
      do i=1,nuc_num
         do ip=1,ctdg_ipol
            call cneo_tdgrad_response_ov(i,0,'neints',ctdg_bfn(i),
     $            ctdg_bfao,-2.0d0,g_x_3d_mo(ip),ctdg_g_mo_n(i),
     $            ctdg_g_mo(ip),ctdg_nao_n(i),ctdg_nao,ctdg_naoc_n(i),
     $            ctdg_naoc(ip),ctdg_nav_n(i),ctdg_nav(ip),0,
     $            ctdg_nfc(ip),0,ctdg_nfv(ip),ctdg_nroots,g_Ax_n_mo(i),
     $            ctdg_geom,ctdg_naoc_n(i),ctdg_nocc(ip),'ib',
     $            ctdg_tol2e_cneo)
         enddo
         do j=1,nuc_num
            if (j.eq.i) goto 10
            call cneo_tdgrad_response_ov(i,j,'nnints',ctdg_bfn(i),
     $            ctdg_bfn(j),2.0d0,g_x_3d_n(j),ctdg_g_mo_n(i),
     $            ctdg_g_mo_n(j),ctdg_nao_n(i),ctdg_nao_n(j),
     $            ctdg_naoc_n(i),ctdg_naoc_n(j),
     $            ctdg_nav_n(i),ctdg_nav_n(j),0,0,0,0,
     $            ctdg_nroots,g_Ax_n_mo(i),
     $            ctdg_geom,ctdg_naoc_n(i),ctdg_naoc_n(j),'ib',
     $            ctdg_tol2e_cneo)

 10         continue
         enddo
      enddo

c     12. add 2g.r to g_Ax_n_mo
c    
      do i=1,nuc_num
         do j=1,ctdg_nav_n(i)
            alo(1)=1
            ahi(1)=3
            alo(2)=1
            ahi(2)=ctdg_naoc_n(i)
            alo(3)=j
            ahi(3)=j
            blo(2)=1
            bhi(2)=3
            blo(1)=1
            bhi(1)=ctdg_nroots
            clo(1)=1
            chi(1)=ctdg_nroots
            clo(2)=1
            chi(2)=ctdg_naoc_n(i)
            clo(3)=j
            chi(3)=j
            call nga_matmul_patch('t','n',2.0d0,1.0d0,
     $                     g_tdg(i),blo,bhi,
     $                     ctdg_g_int1er_ov(i),alo,ahi,
     $                     g_Ax_n_mo(i),clo,chi)
         enddo
      enddo
c
c     13. Transform g_Ax_n_mo -> g_Ax_n 2D MO
c
      do ip = 1, nuc_num
        do ir = 1, ctdg_nroots
          alo(1) = ir
          ahi(1) = ir
          alo(2) = 1
          ahi(2) = ctdg_naoc_n(ip)
          alo(3) = 1
          ahi(3) = ctdg_nav_n(ip)
          klo = 1
          col = ir
          call tddft_grad_mat_to_vec(g_Ax_n_mo(ip),alo,ahi,g_Ax_n(ip),
     +                               klo,col,1.0d0,"=")
        enddo
      enddo

c     14. Calculate g_Ax_g
c         Yg       = \sum{ia} R(i)_{ia}      Zn(i)_{ia}
c        g_Ax_g(i)        g_int1e_r_vo(i)     g_x_n_3D(i) 
c        (3,nroot)        (3,nocc,nvir)     (nroot,nocc,nvir)
      do i=1,nuc_num
         do ip=1,3
            do ir=1,ctdg_nroots
               alo(1)=ip
               ahi(1)=ip
               alo(2)=1
               ahi(2)=ctdg_naoc_n(i)
               alo(3)=1
               ahi(3)=ctdg_nav_n(i)
               clo(1)=ir
               chi(1)=ir
               clo(2)=1
               chi(2)=ctdg_naoc_n(i)
               clo(3)=1
               chi(3)=ctdg_nav_n(i)
               tmp_dbl=nga_ddot_patch(ctdg_g_int1er_ov(i),'n',alo,ahi,
     $                                g_x_3d_n(i),'n',clo,chi)
               call ga_put(g_Ax_g(i),ip,ip,ir,ir,tmp_dbl,1)
            enddo ! end loop nroot
         enddo ! end loop x,y,z
      enddo ! end loop nuc_num
c
c     15. Assemble g_Ax_e,g_Ax_n,g_Ax_g to g_Ax
c
      call cneo_assamble(ctdg_nroots,nuc_num,g_Ax_e,g_Ax_n,g_Ax_g,
     +      ctdg_e_size,ctdg_n_size,g_Ax)

c     16. Destroy local arrays
      call tddft_grad_destroy_r(ctdg_ipol,g_Ax_3d_mo)
      call tddft_grad_destroy_r(ctdg_ipol,g_x_3d_mo)
      call cneo_tdgrad_destroy_r_n(nuc_num,g_x_3d_n)
      call cneo_tdgrad_destroy_r_n(nuc_num,g_Ax_n_mo)
      status = .true.
      status = status .and. ga_destroy(g_x_e)
      status = status .and. ga_destroy(g_Ax_e)
      do i=1,nuc_num
         status = status .and. ga_destroy(g_x_n(i))
         status = status .and. ga_destroy(g_Ax_n(i))
         status = status .and. ga_destroy(g_tdg(i))
         status = status .and. ga_destroy(g_Ax_g(i))
      enddo
      if (.not. status) call errquit(pname//'ga_destroy?',0,0)

      end

      subroutine cneo_disassamble(nroot,nuc_num,g_x,e_size,n_size,
     +      g_x_e,g_x_n,g_tdg)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer nroot
      integer nuc_num
      integer g_x    ![input] concatenated array; (tot_size,nroot)
      integer e_size ![input] elec array size
      integer n_size(nuc_num) ![input] nuclear array size
      integer g_x_e  ![output]  electron part of the concatenated array (e_size,nroot)
      integer g_x_n(nuc_num)  ![output] nuclear part of the concatenated array g_x_n(i):(n_size(i),nroot)
      integer g_tdg(nuc_num)  ![output] constrained multiplier for cneo tdgrad g_tdf(i):(3,nroot)

      integer i
      integer current

      call ga_copy_patch('n',g_x,1,e_size,1,nroot,
     $            g_x_e,1,e_size,1,nroot)
      current=e_size
      do i=1,nuc_num
         call ga_copy_patch('n',g_x,current+1,current+n_size(i),
     $            1,nroot,g_x_n(i),1,n_size(i),1,nroot)
         current = current+n_size(i)
      enddo
      do i=1,nuc_num
         call ga_copy_patch('n',g_x,current+1,current+3,1,nroot,
     $            g_tdg(i),1,3,1,nroot)
         current = current+3
      enddo

      end

      subroutine cneo_assamble(nroot,nuc_num,g_y_e,g_y_n,g_y_g,
     +      e_size,n_size,g_y)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer nroot
      integer nuc_num
      integer g_y_e  ![input]  electron array (e_size,nroot)
      integer g_y_n(nuc_num)  ![input] nuclear array g_x_n(i):(n_size(i),nroot)
      integer g_y_g(nuc_num)  ![input] constrained multiplier array:(3,nroot)
      integer e_size ![input] elec array size
      integer n_size(nuc_num) ![input] nuclear array size
      integer g_y    ![output] concatenated array; (tot_size,nroot)

      integer i
      integer current

      call ga_copy_patch('n',g_y_e,1,e_size,1,nroot,
     $         g_y,1,e_size,1,nroot)

      current=e_size
      do i=1,nuc_num
         call ga_copy_patch('n',g_y_n(i),1,n_size(i),1,nroot,
     $            g_y,current+1,current+n_size(i),1,nroot)
         current = current+n_size(i)
      enddo

      do i=1,nuc_num
         call ga_copy_patch('n',g_y_g(i),1,3,1,nroot,
     $            g_y,current+1,current+3,1,nroot)
         current = current+3
      enddo
      
      end

      subroutine cneo_tdgrad_response_ov(idx1,idx2,inttype,bas1,bas2,
     +            jfac,g_moin,g_mo1,
     +            g_mo2,nao1,nao2,naoc1,naoc2,nav1,nav2,nfc1,nfc2,
     +            nfv1,nfv2,nroot,g_r,geom,nocc1,nocc2,blockout,
     +            tol2e)
      implicit none 
c
c     This subroutine is written because the input g_moin for 
c     cneo_tdgrad_en_response needs to be a symmetric full matrix (nao,nao),
c     but the g_x in cneo_tdgrad_solve_product only has asymmetric ov block.
c     Therefore, the input g_x (nocc, nvir) needs to be transfered to
c
c           | 0       Xia/2 |
c           |                 |  where t(A) is the transpose of A.
c           | t(Xia/2)  0   |
c     
c     cneo_tdgrad_en_response is then called to get the Coulomb response
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer idx1,idx2    ! quantum_nuc index; to read int_2e file
      character*(*) inttype
      integer bas1, bas2,geom
      double precision jfac
      integer g_moin ! [input] input (ov) matrix in MO2 basis: (nroots,nocc2,nvir2)
      integer g_mo1
      integer g_mo2
      integer nao1, naoc1, nav1, nfc1, nfv1, nocc1
      integer nao2, naoc2, nav2, nfc2, nfv2, nocc2
      integer nroot
      integer g_r ! [output] output (ov) matrix in MO1 basis: (nroots,nocc1,nvir1)
      character*(*) blockout
      double precision tol2e
      character*32 pname
      integer g_moin_tmp
      integer idim(3)
      integer ichnk(3)
      integer alo(3),ahi(3)
      integer blo(3),bhi(3)

      pname = 'cneo_tdgrad_response_ov: '

      idim(1) = nroot
      idim(2) = nao2
      idim(3) = nao2
      ichnk(1)= nroot
      ichnk(2)= -1
      ichnk(3)= -1
      if (.not. nga_create(mt_dbl,3,idim,'vectors full mo',ichnk,
     +   g_moin_tmp))
     +   call errquit(pname//'failed to create g_tuv',0,GA_ERR)

      call ga_zero(g_moin_tmp)
      alo(1) = 1
      ahi(1) = nroot
      alo(2) = 1
      ahi(2) = naoc2
      alo(3) = 1
      ahi(3) = nav2
      blo(1) = 1
      bhi(1) = nroot
      blo(2) = 1
      bhi(2) = naoc2
      blo(3) = naoc2+nfc2+1
      bhi(3) = blo(3)+nav2-1
      call nga_copy_patch('n', g_moin, alo, ahi, g_moin_tmp, blo, bhi)
      blo(2) = naoc2+nfc2+1
      bhi(2) = blo(2)+nav2-1
      blo(3) = 1
      bhi(3) = naoc2
      call nga_copy_patch('t', g_moin, alo, ahi, g_moin_tmp, blo, bhi)
      blo(2) = 1
      bhi(2) = nao2
      blo(3) = 1
      bhi(3) = nao2
      call nga_scale_patch(g_moin_tmp, blo, bhi, 5.0d-1)

      call cneo_tdgrad_en_response(idx1,idx2, inttype,bas1,bas2,
     +            jfac,g_moin_tmp,1,'pq',g_mo1,
     +            g_mo2,nao1,nao2,naoc1,naoc2,nav1,nav2,nfc1,nfc2,
     +            nfv1,nfv2,nroot,g_r,geom,nocc1,nocc2,blockout,
     +            tol2e)

      if (.not. ga_destroy(g_moin_tmp)) call errquit
     $      (pname//'destroy temp array',0,GA_ERR)

      end
       