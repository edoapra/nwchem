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

      subroutine get_fock_add_cneo(g_f, g_int1e_r, g_fock)
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

      subroutine cneo_get_fock(nuc_num,g_fock_e,g_fock_n)
      implicit none
      integer nuc_num
      integer g_fock_e
      integer g_fock_n(nuc_num)
      integer i


      call cneo_scf_set('e',0)
      call rohf_get_fock(g_fock_e)

      do i=1,nuc_num
         call cneo_scf_set('n',i)
         call rohf_get_fock(g_fock_n(i))
      enddo

      end