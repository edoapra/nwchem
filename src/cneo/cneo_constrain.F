      subroutine cneo_solve_constrain(idx)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "nwc_const.fh"
#include "ccneo.fh"
#include "ccneo_constrain.fh"
#include "ccneo_dft.fh"

      integer rtdb
      integer g_cfock, g_nfock
      double precision eone, etwo, eone_n, e_en
      integer g_f
      logical status
      integer idx

      external position_analysis

      double precision One
      parameter (One=1.0d0)

!DMR: This is only a 3D-vector, there is no need to use a global array

!      if (.not. ga_create(MT_DBL,3,1,'cneo_constrain:f_multiplier',
!     $        32, 32, g_f)) call errquit('cneo_constrain:ga_create',
!     &       0,GA_ERR)
!      call ga_zero(g_f)
!      call ga_put(g_f, 1,3,1,1,f_cneo(1,idx),3)

      f_cneo(1:3, idx) = 0d0    
      call pos_ana_set(idx, nbf_n, bas_n(idx))
      call util_Newton_solve(f_cneo(1,idx), position_analysis,
     $            1.0d-9, 10, converged_solve)

      call cneo_fock_constrain(g_f,g_int1e_r(1,idx),g_fock_add(idx))

      call ga_get(g_f,1,3,1,1,f_cneo(1,idx),3)

      if (.not. ga_destroy(g_f)) call errquit
     &       ('cneo_constrain:ga_destroy',0,GA_ERR)

      end
      
      
      subroutine position_analysis(acc, g_f, g_r)
      implicit none 
#include "errquit.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "ccneo.fh"
#include "ccneo_constrain.fh"

      double precision acc
      integer g_f
      integer g_r

      integer g_fock_add_local
      integer g_fock_local
      integer type, nrow, ncol
      double precision One
      integer g_mo, g_over_n
      double precision eval(nbf_solve)
      integer g_ndens, g_odens, g_rD
      double precision rD(nbf_solve, nbf_solve)
      double precision exp_pos_local(3)
      integer i,j
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      parameter (One = 1.0d0)

      if (.not. ga_duplicate(g_fock_add(1),g_fock_add_local,'tmp'))
     $      call errquit('pos_ana:ga_duplicate?',0,0)
      call ga_zero(g_fock_add_local)

      call cneo_fock_constrain(g_f, g_int1e_r_solve, g_fock_add_local)

      if (.not. ga_create(mt_dbl, nbf_solve,nbf_solve,'pos_ana: fock',
     $         32, 32, g_fock_local)) call errquit
     &         ('pos_ana:ga_create?',0,0)

      call ga_dadd(One, g_fock_add_local, One, g_fock0_solve, 
     $                                             g_fock_local)

      if (.not. ga_create(mt_dbl, nbf_solve, nbf_solve, 'pos_ana: mo',
     $         32, 32, g_mo)) call errquit
     &         ('pos_ana:ga_create?',0,0)

      if (.not. ga_create(MT_DBL,nbf_solve,nbf_solve,'pos_ana:over',
     $     0, 0, g_over_n)) call
     $     errquit('pos_ana: ga failed for over', 0, GA_ERR)

      if (.not. ga_create(MT_DBL,nbf_solve,nbf_solve,'pos_ana:over',
     $     0, 0, g_rD)) call
     $     errquit('pos_ana: ga create?', 0, GA_ERR)
     
      call ga_zero(g_over_n)
      call int_1e_ga(bas_solve, bas_solve, g_over_n, 
     $                                 'overlap', .false. )

      call ga_diag_seq(g_fock_local, g_over_n, g_mo, eval)

      if (.not. ga_create(MT_DBL,nbf_solve,nbf_solve,'pos_ana:dens',
     $     0, 0, g_ndens)) call
     $     errquit('pos_ana: ga failed for over', 0, GA_ERR)

      call rohf_dens(g_mo, nbf_solve, 1, 0, g_ndens, g_odens)
      call ga_dscal(g_ndens, 5.0d-1)

      do i=1,3
         call ga_zero(g_rD)
         call ga_dgemm('n','n',nbf_solve,nbf_solve,nbf_solve,
     $            One,g_int1e_r_solve(i),g_ndens,0.0d0,g_rD)

         call ga_get(g_rD, 1, nbf_solve, 1, nbf_solve, rD, nbf_solve)

         exp_pos_local(i) = 0.0d0
         do j=1, nbf_solve
            exp_pos_local(i) = exp_pos_local(i)+rD(j,j)
         enddo
      enddo

      call ga_put(g_r, 1,3,1,1,exp_pos_local,1)

      if (.not. ga_destroy(g_fock_local)) call errquit
     $      ('position_analysis:ga_destroy?',0,0)
      if (.not. ga_destroy(g_mo)) call errquit
     $      ('position_analysis:ga_destroy?',0,0)
      if (.not. ga_destroy(g_over_n)) call errquit
     $      ('position_analysis:ga_destroy?',0,0)
      if (.not. ga_destroy(g_rD)) call errquit
     $      ('position_analysis:ga_destroy?',0,0)
      if (.not. ga_destroy(g_ndens)) call errquit
     $      ('position_analysis:ga_destroy?',0,0)
      if (.not. ga_destroy(g_fock_add_local)) call errquit
     $      ('position_analysis:ga_destroy?',0,0)

      end

      subroutine pos_ana_set(idx, nbf, basis)
      implicit none 
#include "nwc_const.fh"
#include "ccneo.fh"
#include "ccneo_constrain.fh"

      integer idx
      integer nbf, basis
      integer i

      do i=1,3
         g_int1e_r_solve(i) = g_int1e_r(i,idx)
      enddo
      g_fock0_solve = g_fock0(idx)
      nbf_solve = nbf
      bas_solve = basis

      end



