      subroutine cneo_tdgrad_solve_init(rtdb,geom,ihdl_bfao,ipol,
     +      nroots,nao,nocc,naoc,nav,nfc,nfv,nmo,g_mo,iptr_mo_e,
     +      tda,oskel,lhashf,otriplet,lhascd,tol2e,kfac,
     +      ihdl_bfn,nao_n,naoc_n,nav_n,nmo_n,g_mo_n,iptr_mo_e_n,
     +      g_int1e_r,tol2e_cneo)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "ccneo_tdgrad_solve.fh"
#include "rtdb.fh"

      integer rtdb
      integer geom
      integer ihdl_bfao,ipol
      integer nroots
      integer nao,nocc(2),naoc(2),nav(2),nfc(2),nfv(2),nmo(2)
      integer g_mo(2)
      integer iptr_mo_e
      logical tda,oskel,lhashf,otriplet,lhascd
      double precision tol2e,kfac
      double precision tol2e_cneo
      
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
      ctdg_tol2e_cneo = tol2e_cneo
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

      call cneo_tdgrad_int1er_ov(nuc_num,g_int1e_r,g_mo_n,
     +            nao_n,naoc_n,nav_n,ctdg_g_int1er_ao,
     +            ctdg_g_int1er_ov)

      end

      subroutine cneo_tdgrad_int1er_ov(nuc_num,g_int1e_r,g_mo,
     +            nao,naoc,nav,g_int1e_r_ao,g_int1e_r_mo)
      implicit none 
#include "errquit.fh"
#include "global.fh"
      integer nuc_num
      integer g_int1e_r(3,nuc_num)
      integer g_mo(nuc_num)
      integer nao(nuc_num)
      integer naoc(nuc_num)
      integer nav(nuc_num)
      integer g_int1e_r_mo(nuc_num) ![output]

      integer i,ii
      integer g_int1e_r_ao(nuc_num)

      integer alo(3),ahi(3)
      integer blo(2),bhi(2)

      call cneo_tdgrad_create_r_n(nuc_num,nao,nao,3,g_int1e_r_ao)
      call cneo_tdgrad_create_r_n(nuc_num,naoc,nav,3,g_int1e_r_mo)

      do i=1,nuc_num
         do ii=1,3
            alo(1)=ii
            ahi(1)=ii
            alo(2)=1
            ahi(2)=nao(i)
            alo(3)=1
            ahi(3)=nao(i)
            blo(1)=1
            bhi(1)=nao(i)
            blo(2)=1
            bhi(2)=nao(i)
            call nga_copy_patch('n',g_int1e_r(ii,i),blo,bhi,
     $                              g_int1e_r_ao(i),alo,ahi)
         enddo

         call tddft_grad_trans_ao2mo(1,nao(i),0,naoc(i),naoc(i),nav(i),
     $                     0,3,1.0d0,0d0,'ib',g_mo(i),g_int1e_r_ao(i),
     $                     g_int1e_r_mo(i),'ib')
      enddo

      end

      subroutine cneo_tdgrad_solve_tidy
      implicit none
#include "global.fh"
#include "errquit.fh"
#include "ccneo_tdgrad_solve.fh"

      integer i
      character*32 pname

      pname = 'cneo_tdgrad_solve_tidy: '

      do i = 1, nuc_num
         if (.not. ga_destroy(ctdg_g_int1er_ao(i))) call errquit
     $         (pname//'failed to destroy g_int1er_ao',0,0)
         if (.not. ga_destroy(ctdg_g_int1er_ov(i))) call errquit
     $         (pname//'failed to destroy g_int1er_ov',0,0)
      enddo

      end
