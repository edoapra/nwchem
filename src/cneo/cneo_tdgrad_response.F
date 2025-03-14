      subroutine cneo_tdgrad_r_n(nuc_num,geom,bas_e,bas_n,nao,naoc,nav,
     +            nocc,nfc,nfv,nao_n,naoc_n,nav_n,g_mo,g_mo_n,g_tp,
     +            nroot,ipol,g_r)
      implicit none 
#include "errquit.fh"
      integer nuc_num,geom,bas_e
      integer bas_n(nuc_num)
      integer nao(2),naoc(2),nav(2),nfc(2),nfv(2),nocc(2)
      integer nao_n(nuc_num),naoc_n(nuc_num),nav_n(nuc_num)
      integer g_mo(2),g_mo_n(nuc_num),g_tp
      integer nroot,ipol,g_r(nuc_num)

      integer i
      character*5 blocksin(2)

      call cneo_tdgrad_create_r_n(nuc_num,naoc_n,nav_n,nroot,g_r)
      print *,'tdgrad:r_n: done creating'
      ! blocksin(1)=' '
      ! blocksin(2)=' '
      blocksin(1)='ij'
      blocksin(2)='ab'
      print *,'tdgrad:r_n:blocksin:',blocksin(1),blocksin(2)
      call ga_zero(g_r)
      call ga_print(g_tp)
      do i=1,nuc_num
   !       call cneo_tdgrad_compute_r_n(i,bas_e,bas_n(i),g_tp,g_mo,
   !   $            g_mo_n(i),nao,nao_n(i),naoc,naoc_n(i),nav,nav_n(i),
   !   $            nfc,nfv,nroot,ipol,g_r(i),geom,nocc)
         call cneo_tdgrad_en_response(i,0,'neints',bas_n(i),bas_e,
     $        -2.0d0,g_tp,2,blocksin,g_mo_n(i),g_mo(1),nao_n(i),nao,
     $            naoc_n(i),naoc(1),nav_n(i),nav(1),0,nfc(1),0,nfv(1),
     $            nroot,g_r(i),geom,naoc_n(i),nocc(1),'ib')
         if (ipol.gt.1) then
            call cneo_tdgrad_en_response(i,0,'neints',bas_n(i),bas_e,
     $        -2.0d0,g_tp,2,blocksin,g_mo_n(i),g_mo(2),nao_n(i),nao,
     $            naoc_n(i),naoc(2),nav_n(i),nav(2),0,nfc(2),0,nfv(2),
     $            nroot,g_r(i),geom,naoc_n(i),nocc(2),'ib')
         endif

      enddo

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
      integer idim(3),ichunk(3)

      ! idim(1)=nroot
      do i=1,nuc_num
         print *,'tdgrad:r_n:create:',i,nuc_num,nroot,naoc(i),nav(i)
      !    idim(2)=naoc(i)
      !    idim(3)=nav(i)
      !    ichunk(1) = nroot
      !   ichunk(2) = -1
      !   ichunk(3) = -1
   !      if (.not.nga_create(mt_dbl,3,idim,'vectors R_n',ichunk,
   !   +                      g_r(i)))
   !   +    call errquit('tddft_grad_create_r: failed to create g_r',
   !   +                 0, GA_ERR)
         call tddft_grad_create_r(1,naoc(i),nav(i),nroot,g_r(i))
         ! g_r(i) = g_tmp(1)
      enddo

      end

!       subroutine cneo_tdgrad_compute_r_n(idx, bas_e,bas_n,g_tp,g_mo,
!      +            g_mo_n,nao,nao_n,naoc,naoc_n,nav,nav_n,nfc,nfv,
!      +            nroot,ipol,g_r,geom,nocc)
!       implicit none 
! #include "errquit.fh"
! #include "mafdecls.fh"
! #include "global.fh"
!       integer idx    ! quantum_nuc index; to read int_2e file
!       integer bas_e, bas_n,geom
!       integer g_tp(2)
!       integer g_mo(2)
!       integer g_mo_n
!       integer nao
!       integer nao_n, naoc_n,nav_n
!       integer naoc(2), nav(2), nfc(2), nfv(2), nocc(2)
!       integer nroot, ipol
!       integer g_r
!       character*32 pname

!       integer g_tuv_pol, g_tuv_tot
!       integer idim(3), ichnk(3)
!       integer idim2(2), ichnk2(2)
!       integer g_tuv_tmp
!       integer g_tuv2  ! 2D matrix
!       integer alo(3),ahi(3)
!       integer blo(3),bhi(3)
!       integer clo(2),chi(2)

!       integer ir,ip
!       integer g_r_ao_2D, g_r_ao
!       integer g_list(10),ig

!       pname="cneo_tdgrad_compute_r_n: "

!       idim(1) = nroot*ipol
!       idim(2) = nao
!       idim(3) = nao
!       ichnk(1) = nroot*ipol
!       ichnk(2) = -1
!       ichnk(3) = -1
!       if (.not.nga_create(mt_dbl,3,idim,'vectors Tuv_pol',ichnk,
!      +   g_tuv_pol))
!      +   call errquit(pname//'failed to create g_tuv',0,GA_ERR)
!       ig=1
!       g_list(ig) = g_tuv_pol
! c
! c     Tuv = sum_ij Cui*Tij*Cvj
! c
!       call tddft_grad_trans_mo2ao(ipol,nao,nfc,naoc,nocc,nav,nfv,
!      +     nroot,1.0d0,0.0d0,"ij",g_mo,g_tp,"pq",g_tuv_pol)
! c
! c     Tuv = Tuv + sum_ab Cua*Tab*Cvb
! c
!       call tddft_grad_trans_mo2ao(ipol,nao,nfc,naoc,nocc,nav,nfv,
!      +     nroot,1.0d0,1.0d0,"ab",g_mo,g_tp,"pq",g_tuv_pol)
! c     sum spin-up and spin-down
!       idim(1) = nroot
!       idim(2) = nao
!       idim(3) = nao
!       ichnk(1) = nroot
!       ichnk(2) = -1
!       ichnk(3) = -1
!       if (.not.nga_create(mt_dbl,3,idim,'vectors Tuv_tot',ichnk,
!      +   g_tuv_tot))
!      +   call errquit(pname//'failed to create g_tuv',0,GA_ERR)
!       if (.not.ga_create(mt_dbl,nao,nao,'vectors Tuv_2D',32, 32,
!      +   g_tuv2))
!      +   call errquit(pname//'failed to create g_tuv',0,GA_ERR)
!       ig=ig+1
!       g_list(ig) = g_tuv_tot
!       ig=ig+1
!       g_list(ig) = g_tuv2

!       call ga_zero(g_tuv_tot)
!       do ir=1,nroot
!          do ip=1,ipol
!             ! call ga_zero(g_tuv_tmp)
!             alo(1)=(ip-1)*nroot+ir
!             alo(2)=1
!             alo(3)=1
!             ahi(1)=alo(1)
!             ahi(2)=nao
!             ahi(3)=nao
!             blo(1)=ir
!             blo(2)=1
!             blo(3)=1
!             bhi(1)=ir
!             bhi(2)=nao
!             bhi(3)=nao
!             call nga_add_patch(1d0,g_tuv_pol,alo,ahi,
!      $            1d0,g_tuv_tot,blo,bhi,g_tuv_tot,blo,bhi)
!          enddo
!       enddo

!       if (.not.ga_create(mt_dbl,nao_n,nao_n,'vectors R_ao2',32, 32,
!      +   g_r_ao_2D))
!      +   call errquit(pname//'failed to create g_r_ao',0,GA_ERR)
!       idim(1) = nroot
!       idim(2) = nao_n
!       idim(3) = nao_n
!       ichnk(1) = nroot
!       ichnk(2) = -1
!       ichnk(3) = -1
!       if (.not.nga_create(mt_dbl,3,idim,'vectors R_ao',ichnk,
!      +   g_r_ao))
!      +   call errquit(pname//'failed to create g_r_ao',0,GA_ERR)

!       ig=ig+1
!       g_list(ig) = g_r_ao_2D
!       ig=ig+1
!       g_list(ig) = g_r_ao

!       call ga_zero(g_r_ao)
!       do ir=1,nroot
!          call ga_zero(g_tuv2)
!          call ga_zero(g_r_ao_2D)
!          alo(1)=ir
!          alo(2)=1
!          alo(3)=1
!          ahi(1)=ir
!          ahi(2)=nao
!          ahi(3)=nao
!          clo(1)=1
!          clo(2)=1
!          chi(1)=nao
!          chi(2)=nao
!          call nga_copy_patch('n',g_tuv_tot,alo,ahi,g_tuv2,clo,chi)
!          call get_j_dm(geom,bas_n,bas_e,-1.0d0,1.0d-9,.false.,g_tuv2,
!      $            g_r_ao_2D,'neints',idx,0,.true.) ! here, jfac is 1.0 or -1.0? check later
!          clo(1)=1
!          clo(2)=1
!          chi(1)=nao_n
!          chi(2)=nao_n
!          alo(1)=ir
!          alo(2)=1
!          alo(3)=1
!          ahi(1)=ir
!          ahi(2)=nao_n
!          ahi(3)=nao_n
!          call nga_copy_patch('n',g_r_ao_2D,clo,chi,g_r_ao,alo,ahi)
!       enddo

!       call tddft_grad_trans_ao2mo(1,nao_n,0,naoc_n,1,nav_n,0,nroot,
!      $         1.0d0,0.0d0,'ib',g_mo_n,g_r_ao,g_r,'ib')

!       do ir=1,ig
!          if (.not.ga_destroy(g_list(ir))) call errquit
!      $         (pname//'ga_destroy?',0,0)
!       enddo

!       print *,'tdgrad:r_n:',idx
!       call ga_print(g_r)

!       end

      subroutine cneo_tdgrad_en_response(idx1,idx2, inttype,bas1,bas2,
     +            jfac,g_moin,nblocksin,blocksin,g_mo1,
     +            g_mo2,nao1,nao2,naoc1,naoc2,nav1,nav2,nfc1,nfc2,
     +            nfv1,nfv2,nroot,g_r,geom,nocc1,nocc2,blockout)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer idx1,idx2    ! quantum_nuc index; to read int_2e file
      character*(*) inttype
      integer bas1, bas2,geom
      double precision jfac
      integer g_moin ! [input] input matrix in MO basis
      integer nblocksin
      character*(*) blocksin(nblocksin)
      integer g_mo1
      integer g_mo2
      integer nao1, naoc1, nav1, nfc1, nfv1, nocc1
      integer nao2, naoc2, nav2, nfc2, nfv2, nocc2
      integer nroot
      integer g_r
      character*2 blockout
      character*32 pname

      integer g_tuv
      integer idim(3), ichnk(3)
      integer idim2(2), ichnk2(2)
      ! integer g_tuv_tmp
      integer g_tuv2,g_tuv2_t  ! 2D matrix
      integer alo(3),ahi(3)
      integer blo(3),bhi(3)
      integer clo(2),chi(2)

      integer ir,ip,ib
      integer g_r_ao_2D, g_r_ao, g_r_ao_2d_t
      integer g_list(15),ig

      pname="cneo_tdgrad_compute_r_n: "
      print *,'tdgrad:en response:entering...',nblocksin
      ! print *,nblocksin,blocksin(1)

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

      ! print *,pname,'line 286'
      ! print *,nao2,nfc2,naoc2,nocc2,nav2,nfv2,nroot
      ! print *,g_mo2,g_moin,g_tuv
      ! call ga_print(g_moin)
      ! call ga_print(g_mo2)
      
      call tddft_grad_trans_mo2ao(1,nao2,nfc2,naoc2,nocc2,nav2,
     +     nfv2,nroot,1.0d0,0.0d0,blocksin(1),g_mo2,g_moin,"pq",g_tuv)
      ! print *,pname,'line 291'
      ! call ga_print(g_tuv)
      if (nblocksin .gt.1) then
         do ib=2,nblocksin
c
c     Tuv = Tuv + sum_ab Cua*Tab*Cvb
c
            call tddft_grad_trans_mo2ao(1,nao2,nfc2,naoc2,nocc2,nav2,
     +     nfv2,nroot,1.0d0,1.0d0,blocksin(ib),g_mo2,g_moin,"pq",g_tuv)
         enddo
      endif
      ! print *,pname,'line 294'
! c     sum spin-up and spin-down
!       idim(1) = nroot
!       idim(2) = nao
!       idim(3) = nao
!       ichnk(1) = nroot
!       ichnk(2) = -1
!       ichnk(3) = -1
!       if (.not.nga_create(mt_dbl,3,idim,'vectors Tuv_tot',ichnk,
!      +   g_tuv_tot))
!      +   call errquit(pname//'failed to create g_tuv',0,GA_ERR)
      if (.not.ga_create(mt_dbl,nao2,nao2,'vectors Tuv_2D',32, 32,
     +   g_tuv2))
     +   call errquit(pname//'failed to create g_tuv',0,GA_ERR)
      if (.not.ga_create(mt_dbl,nao2,nao2,'vectors Tuv_2D',32, 32,
     +   g_tuv2_t))
     +   call errquit(pname//'failed to create g_tuv',0,GA_ERR)
      ! ig=ig+1
      ! g_list(ig) = g_tuv_tot
      ig=ig+1
      g_list(ig) = g_tuv2

      ! call ga_zero(g_tuv)
   !    do ir=1,nroot
   !       do ip=1,ipol
   !          ! call ga_zero(g_tuv_tmp)
   !          alo(1)=(ip-1)*nroot+ir
   !          alo(2)=1
   !          alo(3)=1
   !          ahi(1)=alo(1)
   !          ahi(2)=nao
   !          ahi(3)=nao
   !          blo(1)=ir
   !          blo(2)=1
   !          blo(3)=1
   !          bhi(1)=ir
   !          bhi(2)=nao
   !          bhi(3)=nao
   !          call nga_add_patch(1d0,g_tuv_pol,alo,ahi,
   !   $            1d0,g_tuv_tot,blo,bhi,g_tuv_tot,blo,bhi)
   !       enddo
   !    enddo

      if (.not.ga_create(mt_dbl,nao1,nao1,'vectors R_ao2',32, 32,
     +   g_r_ao_2D))
     +   call errquit(pname//'failed to create g_r_ao',0,GA_ERR)
      if (.not.ga_create(mt_dbl,nao1,nao1,'vectors R_ao2',32, 32,
     +   g_r_ao_2D_t))
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
      ig=ig+1
      g_list(ig) = g_r_ao_2d_t
      ig=ig+1
      g_list(ig) = g_tuv2_t
      ! print *,pname,'line 351'

      call ga_zero(g_r_ao)
      do ir=1,nroot
         call ga_zero(g_tuv2)
         call ga_zero(g_r_ao_2D)
         call ga_zero(g_r_ao_2d_t)
         call ga_zero(g_tuv2_t)
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
         ! print *,'tdgrad:debug:test:response:',idx1,idx2,inttype
         ! print *,'input:',bas1,bas2,jfac
         call get_j_dm_asymm(geom,bas1,bas2,jfac,1.0d-9,.false.,g_tuv2,
     $            g_r_ao_2D,inttype,idx1,idx2,.true.)
   !       call ga_transpose(g_tuv2,g_tuv2_t)
   !       call get_j_dm(geom,bas1,bas2,jfac,1.0d-9,.false.,g_tuv2_t,
   !   $            g_r_ao_2D_t,inttype,idx1,idx2,.true.)
   !       call ga_add(5.0d-1,g_r_ao_2D,5.0d-1,g_r_ao_2d_t,g_r_ao_2D)
   !       print *,'tdgrad:debug:asymm'
         ! call ga_print(g_r_ao_2D)
         ! call ga_print(g_tuv2)
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

      ! print *,pname,'line 383',nfc1,naoc1,nocc1,nav1,nfv1
      ! Note that facmo is set to be 1, so that original g_r is added
      ! to the final result. Zero g_r in advance if needed
      call tddft_grad_trans_ao2mo(1,nao1,nfc1,naoc1,nocc1,nav1,nfv1,
     $         nroot,1.0d0,1.0d0,blockout,g_mo1,g_r_ao,g_r,blockout)

      do ir=1,ig
         if (.not.ga_destroy(g_list(ir))) call errquit
     $         (pname//'ga_destroy?',0,0)
      enddo

      ! print *,'tdgrad:r_n:',idx1, idx2
      ! call ga_print(g_r)

      end