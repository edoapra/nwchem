      subroutine cneo_dft_init(rtdb, geomin, basisin, basin_cneo)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "nwc_const.fh"
#include "util.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "cscf.fh"
#include "cfock.fh"
#include "ccneo.fh"
#include "ccneo_constrain.fh"
#include "ccneo_dft.fh"

      integer rtdb
      integer geomin
      integer basisin, basin_cneo(max_quantum_nuc)
      double precision nuclear_charge
      integer len_occ
      integer i

      bas_e = basisin
      bas_n = basin_cneo
      geom_cneo = geomin

      if (.not. bas_numbf(bas_e,nbf_e)) call errquit
     $      ('cneo_dft_init:bas_e',0,0)
      if (.not. bas_numbf(bas_n(1),nbf_n)) call errquit
     $      ('cneo_dft_init:bas_n',0,0)
      nmo_n = nbf_n

      do i=1,nuc_num
         if (.not. ga_create(MT_DBL, nbf_n, nmo_n, 'cneo_dft: n_MOs',
     $        32, 32, g_movecs_n(i))) call errquit('cneo_init: MOs', 0,
     &       GA_ERR)
         if (.not. ga_create(MT_DBL, nbf_n, nmo_n, 'cneo_dft: fock0',
     $        32, 32, g_fock0(i))) call errquit('cneo_init: fock0', 0,
     &       GA_ERR)
         if (.not. ga_create(MT_DBL, nbf_n,nmo_n,'cneo_dft:fock_add',
     $        32, 32, g_fock_add(i))) call errquit('cneo_init: fock0',
     &        0,GA_ERR)
      enddo

      len_occ = nmo_n
      do i=1,nuc_num
         if (.not. ma_push_get(mt_dbl,len_occ,'cneo_dft: mo evals n',
     $     l_eval_n(i), k_eval_n(i))) call errquit
     $     ('cneo_init: insufficient memory?', len_occ, MA_ERR)
c     
         if (.not. ma_push_get(mt_dbl,len_occ,'cneo_dft: mo occ n',
     $     l_occ_n(i), k_occ_n(i))) call errquit
     $     ('cneo_init: insufficient memory?', len_occ, MA_ERR)
      enddo

      end 

      subroutine cneo_dft_tidy
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "nwc_const.fh"
#include "ccneo.fh"
#include "ccneo_dft.fh"
#include "ccneo_constrain.fh"
      integer i
      character*32 pname

      pname = 'cneo_dft_tidy: '

      do i=nuc_num,1,-1
         if (.not. ga_destroy(g_movecs_n(i))) call errquit
     +      (pname//'could not destroy n MOs',0,GA_ERR)
         if (.not. ga_destroy(g_fock0(i))) call errquit
     +      (pname//'could not destroy n fock0',0,GA_ERR)
         if (.not. ga_destroy(g_fock_add(i))) call errquit
     +      (pname//'could not destroy n fock_add',0,GA_ERR)

         if (.not. ma_pop_stack(l_occ_n(i))) call errquit
     +      (pname//'could not deallocate n occ',0,MA_ERR)
         if (.not. ma_pop_stack(l_eval_n(i))) call errquit
     +      (pname//'could not deallocate n eval',0,MA_ERR)
      enddo

      end