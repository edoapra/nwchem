      subroutine cneo_scf_init(rtdb, geomin, basisin, basin_cneo)
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
#include "cconstrain.fh"
#include "cscf_cneo.fh"

      integer rtdb
      integer geomin
      integer basisin, basin_cneo(max_quantum_nuc)
      double precision nuclear_charge
      integer len_occ
      integer i

      bas_e = basisin
      bas_n = basin_cneo
      geom_cneo = geomin

      if (ccneo_init_flag.gt.0)
     $     call errquit('cneo internals already initialised?',0,
     &       UNKNOWN_ERR)

      ccneo_init_flag = 1
      if (.not. bas_numbf(bas_e,nbf_e)) call errquit
     $      ('cneo:fock:bas_e',0,0)
      if (.not. bas_numbf(bas_n(1),nbf_n)) call errquit
     $      ('cneo:fock:bas_n',0,0)
      ! print *, 'cneo:nbf:',nbf_e,nbf_n
      nmo_e = nbf_e
      nmo_n = nbf_n

      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge_cneo)) then
         call errquit('cneo_init:no quantum nuc set up',0,0)
      end if
      if (.not. geom_nuc_charge(geomin, nuclear_charge))
     $     call errquit('scf: geom_nuc_charge failed', 0, GEOM_ERR)
      nelec_e = nint(nuclear_charge - charge_cneo)
      if (nelec_e .le. 0) call errquit
     $     ('scf: negative no. of electrons ?', nelec, INPUT_ERR)

      nopen_e = 0
      nclosed_e = (nelec_e-nopen_e) / 2

      ccneo_vlen_e = (nmo_e-nclosed_e-nopen_e)*
     $         (nclosed_e+nopen_e)+(nclosed_e*nopen_e)
      ccneo_vlen_n = nmo_n-1
      ccneo_vlen = ccneo_vlen_e + ccneo_vlen_n*nuc_num

      if (ccneo_vlen .eq. 0) ccneo_vlen = 1

      if (.not. ga_create(MT_DBL,nmo_e,nmo_e,'cneo_vec_init: FCV_e',
     $     32, 32, ccneo_g_fcv_e)) call errquit('cneo_vec_init: FCV',
     &        0,GA_ERR)
c     
      if (nopen_e .gt. 0) then
         if (.not. ga_create(MT_DBL, nmo_e, nmo_e, 'cneo_vec_init: FPV',
     $        32, 32, ccneo_g_fpv)) call errquit('cneo_init: FPV', 0,
     &       GA_ERR)
c     
         if (.not. ga_create(MT_DBL, nmo_e, nmo_e, 'cneo_vec_init: FCP',
     $        32, 32, ccneo_g_fcp)) call errquit('cneo_init: FCP', 0,
     &       GA_ERR)
      endif

      do i=1,nuc_num
         if (.not. ga_create(MT_DBL,nmo_n,nmo_n,'cneo_vec_init: FCV_n',
     $     32, 32, ccneo_g_fcv_n(i))) call errquit('cneo_vec_init: FCV',
     &        0,GA_ERR)
      enddo

      if (util_print('cneo_scf_init', print_debug)) then
         write(6,*) 'ccneo_vlen,e/n:', ccneo_vlen, ccneo_vlen_e,
     $                        ccneo_vlen_n
         write(6,*) 'nelec, nclosed=', nelec_e, nclosed_e
      end if

      if (.not. ga_create(MT_DBL, nbf_e, nmo_e, 'cneo_init: e_MOs',
     $        32, 32, g_movecs_e)) call errquit('cneo_init: MOs', 0,
     &       GA_ERR)

      do i=1,nuc_num
         if (.not. ga_create(MT_DBL, nbf_n, nmo_n, 'cneo_init: n_MOs',
     $        32, 32, g_movecs_n(i))) call errquit('cneo_init: MOs', 0,
     &       GA_ERR)
         if (.not. ga_create(MT_DBL, nbf_n, nmo_n, 'cneo_init: fock0',
     $        32, 32, g_fock0(i))) call errquit('cneo_init: fock0', 0,
     &       GA_ERR)
         if (.not. ga_create(MT_DBL, nbf_n,nmo_n,'cneo_init:fock_add',
     $        32, 32, g_fock_add(i))) call errquit('cneo_init: fock0',
     &        0,GA_ERR)
      enddo


      len_occ = nmo_e
      if (.not. ma_push_get(mt_dbl, len_occ, 'cneo_init: mo evals e',
     $     l_eval_e, k_eval_e)) call errquit
     $     ('cneo_init: insufficient memory?', len_occ, MA_ERR)
c     
      if (.not. ma_push_get(mt_dbl, len_occ, 'cneo_init: mo occ e',
     $     l_occ_e, k_occ_e)) call errquit
     $     ('cneo_init: insufficient memory?', len_occ, MA_ERR)
c
      if (.not. ma_push_get(mt_int, len_occ, 'cneo_init: mo irs',
     $     l_irs, k_irs)) call errquit
     $     ('cneo_init: insufficient memory?', len_occ, MA_ERR)
c
      call ifill(len_occ, 1, int_mb(k_irs), 1) ! In case not adapting

      len_occ = nmo_n
      do i=1,nuc_num
         if (.not. ma_push_get(mt_dbl,len_occ,'cneo_init: mo evals n',
     $     l_eval_n(i), k_eval_n(i))) call errquit
     $     ('cneo_init: insufficient memory?', len_occ, MA_ERR)
c     
         if (.not. ma_push_get(mt_dbl,len_occ,'cneo_init: mo occ n',
     $     l_occ_n(i), k_occ_n(i))) call errquit
     $     ('cneo_init: insufficient memory?', len_occ, MA_ERR)
      enddo

      call scf_get_conv_info(rtdb)

      do i=1,6
         shifts_e(i) = shifts(i)
      enddo

      call dfill(6, -1.0d0, shifts_n, 1)
      if (shifts_n(1) .eq. -1.0d0) shifts_n(1) = 5.0d0
      if (shifts_n(2) .eq. -1.0d0) shifts_n(2) = 0.5d0
      if (shifts_n(3) .eq. -1.0d0) shifts_n(3) = 0.0d0
      if (shifts_n(4) .eq. -1.0d0) shifts_n(4) = 0.0d0
      if (shifts_n(5) .eq. -1.0d0) shifts_n(5) = 0.0d0
      if (shifts_n(6) .eq. -1.0d0) shifts_n(6) = 0.0d0

      gnorm_tol_cneo = 1.0d-5
      maxiter_cneo = 30
      tol2e_cneo = min(1.0d-9,gnorm_tol*1d-2)

      if (.not. rtdb_cput(rtdb,'task:theory',1,'scf'))
     $      call errquit('cneo:fock:init theory',0,0)
      if (.not. rtdb_put(rtdb, 'bgj:scf_type', MT_INT, 1, 1))
     $     call errquit('cneo:fock:put bgj:scf_type failed',0, RTDB_ERR)

      maxquartet = 10000
      memsize=-1
      maxeri = 1296*100
      maxscr = 0
      oreplicated = .true.
      incmax = -1
      incmin = -1
      incomplete = -1

      end 

      subroutine cneo_scf_set(partical_type, idx)
      implicit none 
#include "errquit.fh"
#include "nwc_const.fh"
#include "ccneo.fh"
#include "cscf_cneo.fh"
#include "cscf.fh"
#include "crohf.fh"

      character partical_type
      integer idx

      geom = geom_cneo
      charge = charge_cneo
      crohf_init_flag = ccneo_init_flag
      nelec = nelec_e

      gnorm_tol = gnorm_tol_cneo
      iter = iter_cneo
      maxiter = maxiter_cneo
      tol2e = tol2e_cneo

      if (partical_type .eq. 'e') then
         call cneo_scf_set_elec
      else if (partical_type .eq. 'n') then
         call cneo_scf_set_nuc(idx)
      else
         call errquit('cneo_scf_set: invalid input',0,0)
      end if

      end 

      subroutine cneo_scf_set_elec
      implicit none 
#include "errquit.fh"
#include "nwc_const.fh"
#include "cscf.fh"
#include "crohf.fh"
#include "ccneo.fh"
#include "cscf_cneo.fh"

      integer i

      basis = bas_e
      nopen = nopen_e
      nclosed = nclosed_e
      nbf = nbf_e
      nmo = nmo_e

      nalpha = nclosed + nopen
      nbeta = nclosed

      do i=1,6
         shifts(i) = shifts_e(i)
      enddo

      g_movecs = g_movecs_e
      l_eval = l_eval_e
      k_eval = k_eval_e
      l_occ = l_occ_e
      k_occ = k_occ_e

      crohf_g_fcv = ccneo_g_fcv_e

      end

      subroutine cneo_scf_set_nuc(idx)
      implicit none 
#include "errquit.fh"
#include "nwc_const.fh"
#include "cscf.fh"
#include "crohf.fh"
#include "ccneo.fh"
#include "cscf_cneo.fh"

      integer i, idx

      basis = bas_n(idx)
      nopen = 0
      nclosed = 1
      nbf = nbf_n
      nmo = nmo_n

      nalpha = nclosed + nopen
      nbeta = nclosed

      do i=1,6
         shifts(i) = shifts_n(i)
      enddo

      g_movecs = g_movecs_n(idx)
      l_eval = l_eval_n(idx)
      k_eval = k_eval_n(idx)
      l_occ = l_occ_n(idx)
      k_occ = k_occ_n(idx)

      crohf_g_fcv = ccneo_g_fcv_n(idx)

      end