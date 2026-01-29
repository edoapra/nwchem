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

      subroutine cneo_vecs_guess_nuc(rtdb, i, iatom, g_edens, work, n)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "ccneo.fh"
#include "ccneo_dft.fh"

      integer rtdb,i, iatom, n
      integer g_edens
      integer me

      integer g_nhcore, g_nfock, g_over, g_tmp
      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      double precision work(n**2, 3)

!DMR: Get the node id
      me = ga_nodeid()

      g_nfock = ga_create_atom_blocked(geom_cneo, bas_n(i), 
     $            'Temp:n:fock')

      call ga_zero(g_nfock)
      call get_j_dm(geom_cneo, bas_n(i), bas_e, -1.0d0, tol2e_cneo,
     $       .false., g_edens, g_nfock, 'neints', i,0,.false.)
      call hcore_nuc(geom_cneo, bas_n(i), g_core_n(i), iatom)
      call ga_dadd(1.0d0, g_core_n(i), 1.0d0, g_nfock, g_nfock)
      if (util_print('cneo_vec_nuc', print_debug)) 
     $         call ga_print(g_nfock)

      call cneo_diag_wrapper(g_nfock, g_ortho_n(i), g_movecs_n(i), work,
     $                       dbl_mb(k_eval_n(i)), n)

      if (.not. ga_destroy(g_nfock)) call errquit
     $      ('cneo_vecs_guess_nuc: destroy ga',0,GA_ERR)
c      if (.not. ga_destroy(g_nhcore)) call errquit
c     $      ('cneo_vecs_guess_nuc: destroy ga',0,GA_ERR)

      if (util_print('cneo_vec_guess_nuc', print_debug)) 
     $         call ga_print(g_movecs_n(i))

      end
