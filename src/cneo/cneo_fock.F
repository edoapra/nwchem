      subroutine hcore_nuc(geom, bas_nuc, g_nhcore,iatom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "geomP.fh"

      integer geom, bas_nuc, g_nhcore
      double precision masselec, avogadro, mass_au
      double precision scale
      integer iatom

      parameter(masselec=9.1093826d-31)
      parameter(avogadro=6.022140857d23)

      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      mass_au = 1.0d-3/avogadro

      g_nhcore = ga_create_atom_blocked(geom,bas_nuc,'cneo:n:hcore')
      scale = geom_mass(iatom,geom)*mass_au/masselec
      call ga_zero(g_nhcore)
      call int_1e_ga(bas_nuc, bas_nuc, g_nhcore, 'kinetic',.false.)
      call ga_dscal(g_nhcore,-1.0d0/scale)
      call int_1e_ga(bas_nuc, bas_nuc,g_nhcore, 'potential',.false.)
      call ga_dscal(g_nhcore, -1.0d0)

      end

      subroutine cneo_fock_constrain(g_f, g_int1e_r, g_fock)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

      integer g_f    ! array handle for multiplier: f
      double precision f(3)
      integer g_int1e_r(3)
      integer g_fock
      integer nrow, ncol
      integer type
      integer i
      double precision One
      parameter (One=1.0d0)
      
      call ga_zero(g_fock)

      call ga_get(g_f, 1, 3, 1, 1, f, 1)

      do i=1,3
         call ga_dadd(f(i),g_int1e_r(i),One,g_fock,g_fock)
      enddo
      
      end

      subroutine cneo_vecs_guess_nuc(rtdb, i,iatom,g_edens)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "ccneo.fh"
#include "ccneo_dft.fh"

      integer rtdb,i, iatom
      integer g_edens
      integer me

      integer g_nhcore, g_nfock, g_over, g_tmp
      integer ga_create_atom_blocked
      external ga_create_atom_blocked

!DMR: Get the node id
      me = ga_nodeid()

!DMR: If this is a single atom, this array will be allocated
!     on a single rank
      g_nfock = ga_create_atom_blocked(geom_cneo, bas_n(i), 
     $            'Temp:n:fock')

!DMR: This might result in too much fragmentation
!      if (.not. ga_create(MT_DBL,nbf_n,nbf_n,'cneo_vecs_guess_n:over',
!     $     0, 0, g_over)) call
!     $     errquit('cneo_vecs_guess_n: ga failed for over', 0, GA_ERR)
      if(.not.ga_duplicate(g_nfock, g_over, "cneo_vecs_guess_n: over"))
     $  call errquit("cneo_vecs_guess_nuc: failed to create array",89,
     $                GA_ERR) 

!      if (.not. ga_create(MT_DBL,nbf_n,nbf_n,'cneo_vecs_guess_n:tmp',
!     $     0, 0, g_tmp)) call
!     $     errquit('cneo_vecs_guess_n: ga failed for tmp', 0, GA_ERR)
      if(.not.ga_duplicate(g_nfock, g_tmp, "cneo_vecs_guess_n: over"))
     $  call errquit("cneo_vecs_guess_nuc: failed to create array",96,
     $                GA_ERR) 

      call ga_zero(g_nfock)
      call get_j_dm(geom_cneo, bas_n(i), bas_e, -1.0d0, tol2e_cneo,
     $       .false., g_edens, g_nfock, 'neints', i,0,.false.)
      call hcore_nuc(geom_cneo, bas_n(i), g_nhcore, iatom)
      call ga_dadd(1.0d0, g_nhcore, 1.0d0, g_nfock, g_nfock)
      if (util_print('cneo_vec_nuc', print_debug)) 
     $         call ga_print(g_nfock)

      call ga_zero(g_over)
      call int_1e_ga(bas_n(i), bas_n(i), g_over, 'overlap', .false.)

!DMR: This function is deprecated
!      call ga_diag(g_nfock, g_over, g_tmp, dbl_mb(k_eval_n(i)))

!DMR: Since each nulcear Fock matrix should be rather
!     small, there is no need to use parallel eigensolvers.
      call ga_diag_seq(g_nfock, g_over, g_tmp, dbl_mb(k_eval_n(i)))
         
!
      call ga_copy(g_tmp, g_movecs_n(i))

      if (.not. ga_destroy(g_tmp)) call errquit
     $      ('cneo_vecs_guess_nuc: destroy ga',0,GA_ERR)
      if (.not. ga_destroy(g_over)) call errquit
     $      ('cneo_vecs_guess_nuc: destroy ga',0,GA_ERR)
      if (.not. ga_destroy(g_nfock)) call errquit
     $      ('cneo_vecs_guess_nuc: destroy ga',0,GA_ERR)
      if (.not. ga_destroy(g_nhcore)) call errquit
     $      ('cneo_vecs_guess_nuc: destroy ga',0,GA_ERR)

      if (util_print('cneo_vec_guess_nuc', print_debug)) 
     $         call ga_print(g_movecs_n(i))

      end
