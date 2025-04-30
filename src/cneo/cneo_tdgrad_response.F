      subroutine cneo_tdgrad_r_n(nuc_num,geom,bas_e,bas_n,nao,naoc,nav,
     +            nocc,nfc,nfv,nao_n,naoc_n,nav_n,g_mo,g_mo_n,g_tp,
     +            nroot,ipol,g_r,tol2e)
      implicit none 
#include "errquit.fh"
#include "global.fh"
#include "util.fh"
      integer nuc_num,geom,bas_e
      integer bas_n(nuc_num)
      integer nao(2),naoc(2),nav(2),nfc(2),nfv(2),nocc(2)
      integer nao_n(nuc_num),naoc_n(nuc_num),nav_n(nuc_num)
      integer g_mo(2),g_mo_n(nuc_num),g_tp
      integer nroot,ipol,g_r(nuc_num)
      double precision tol2e

      integer i
      character*2 blocksin(2)

      call cneo_tdgrad_create_r_n(nuc_num,naoc_n,nav_n,nroot,g_r)
      blocksin(1)='ij'
      blocksin(2)='ab'
      call ga_zero(g_r)
      do i=1,nuc_num
         call cneo_tdgrad_en_response(i,0,'neints',bas_n(i),bas_e,
     $        -2.0d0,g_tp,2,blocksin,g_mo_n(i),g_mo(1),nao_n(i),nao,
     $            naoc_n(i),naoc(1),nav_n(i),nav(1),0,nfc(1),0,nfv(1),
     $            nroot,g_r(i),geom,naoc_n(i),nocc(1),'ib',tol2e)
      enddo

      if (util_print('tddft_grad_r', print_debug)) then
         do i=1,nuc_num
            call ga_print(g_r(i))
         enddo
         call util_flush(6)
      endif

      end
      
      
      subroutine cneo_tdgrad_create_r_n(nuc_num,naoc,nav,nroot,g_r)
c     To create R_n=Q^n_{ia}
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

      integer nuc_num
      integer naoc(nuc_num)
      integer nav(nuc_num)
      integer nroot
      integer g_r(nuc_num)

      integer i

      do i=1,nuc_num
         call tddft_grad_create_r(1,naoc(i),nav(i),nroot,g_r(i))
      enddo

      end
c
c     Coulomb interaction response matrix in MO1 basis from 
c     input matrix in MO2 basis
c     R_{PQ} = \sum_{pq} (PQ|pq) Zpq
c
c     Notes: 
c     1. this subroutine requires the input g_moin is symmetric
c        (this is actually the requirement for get_j_dm)
c     2. The final result is added to the original g_r;
c        zero g_r in advance if needed
c
      subroutine cneo_tdgrad_en_response(idx1,idx2, inttype,bas1,bas2,
     +            jfac,g_moin,nblocksin,blocksin,g_mo1,
     +            g_mo2,nao1,nao2,naoc1,naoc2,nav1,nav2,nfc1,nfc2,
     +            nfv1,nfv2,nroot,g_r,geom,nocc1,nocc2,blockout,tol2e)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer idx1,idx2       ! [input] quantum_nuc index; to read int_2e file
      character*(*) inttype   ! [input] interaction type: 'neints'/'enints'/'nnints'
      integer bas1, bas2,geom
      double precision jfac   ! [input] scaling factor
      integer g_moin          ! [input] input matrix in MO basis
      integer nblocksin 
      character*(*) blocksin(nblocksin)   ![input] blocks for g_moin
      integer g_mo1           ! [input] MO coeffs for bas1
      integer g_mo2           ! [input] MO coeffs for bas2
      integer nao1, naoc1, nav1, nfc1, nfv1, nocc1
      integer nao2, naoc2, nav2, nfc2, nfv2, nocc2
      integer nroot
      integer g_r             ! [output] response matrix in MO2 basis
      character*2 blockout    ! [input] output block for g_r

c     Local variables
      character*32 pname

      integer g_tuv
      integer idim(3), ichnk(3)
      integer idim2(2), ichnk2(2)
      integer g_tuv2,g_tuv2_t  ! 2D matrix
      integer alo(3),ahi(3)
      integer blo(3),bhi(3)
      integer clo(2),chi(2)

      integer ir,ip,ib
      integer g_r_ao_2D, g_r_ao, g_r_ao_2d_t
      integer g_list(15),ig

      double precision tol2e

      pname="cneo_tdgrad_compute_r_n: "

      idim(1) = nroot
      idim(2) = nao2
      idim(3) = nao2
      ichnk(1) = nroot
      ichnk(2) = -1
      ichnk(3) = -1
      if (.not.nga_create(mt_dbl,3,idim,'vectors Tuv',ichnk,
     +   g_tuv))
     +   call errquit(pname//'failed to create g_tuv',0,GA_ERR)
      ig=1
      g_list(ig) = g_tuv
      call ga_zero(g_tuv)

      do ib=1,nblocksin
         call tddft_grad_trans_mo2ao(1,nao2,nfc2,naoc2,nocc2,nav2,
     +     nfv2,nroot,1.0d0,1.0d0,blocksin(ib),g_mo2,g_moin,"pq",g_tuv)
      enddo
      
      if (.not.ga_create(mt_dbl,nao2,nao2,'vectors Tuv_2D',32, 32,
     +   g_tuv2))
     +   call errquit(pname//'failed to create g_tuv',0,GA_ERR)
      ig=ig+1
      g_list(ig) = g_tuv2

      if (.not.ga_create(mt_dbl,nao1,nao1,'vectors R_ao2',32, 32,
     +   g_r_ao_2D))
     +   call errquit(pname//'failed to create g_r_ao',0,GA_ERR)
      idim(1) = nroot
      idim(2) = nao1
      idim(3) = nao1
      ichnk(1) = nroot
      ichnk(2) = -1
      ichnk(3) = -1
      if (.not.nga_create(mt_dbl,3,idim,'vectors R_ao',ichnk,
     +   g_r_ao))
     +   call errquit(pname//'failed to create g_r_ao',0,GA_ERR)

      ig=ig+1
      g_list(ig) = g_r_ao_2D
      ig=ig+1
      g_list(ig) = g_r_ao

      call ga_zero(g_r_ao)
      do ir=1,nroot
         call ga_zero(g_tuv2)
         call ga_zero(g_r_ao_2D)
         alo(1)=ir
         alo(2)=1
         alo(3)=1
         ahi(1)=ir
         ahi(2)=nao2
         ahi(3)=nao2
         clo(1)=1
         clo(2)=1
         chi(1)=nao2
         chi(2)=nao2
         call nga_copy_patch('n',g_tuv,alo,ahi,g_tuv2,clo,chi)
         call get_j_dm(geom,bas1,bas2,jfac,tol2e,.false.,g_tuv2,
     $            g_r_ao_2D,inttype,idx1,idx2,.true.)
         clo(1)=1
         clo(2)=1
         chi(1)=nao1
         chi(2)=nao1
         alo(1)=ir
         alo(2)=1
         alo(3)=1
         ahi(1)=ir
         ahi(2)=nao1
         ahi(3)=nao1
         call nga_copy_patch('n',g_r_ao_2D,clo,chi,g_r_ao,alo,ahi)
      enddo

      call tddft_grad_trans_ao2mo(1,nao1,nfc1,naoc1,nocc1,nav1,nfv1,
     $         nroot,1.0d0,1.0d0,blockout,g_mo1,g_r_ao,g_r,blockout)

      do ir=1,ig
         if (.not.ga_destroy(g_list(ir))) call errquit
     $         (pname//'ga_destroy?',0,0)
      enddo

      end