      subroutine cneo_schwarz_init(geom, nbas, bases)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "cschwarz.fh"
#include "cschwarz_cneo.fh"
#include "schwarz.fh"
#include "basP.fh"

      integer geom
      integer nbas
      integer bases(nbas)

      integer basis, basisin
      integer i

      do i=1,nbas
         basisin = bases(i)
         call schwarz_init(geom, basisin)
         basis = basisin + Basis_Handle_Offset
         sch_maxs(basis) = sch_max
         nb_per_is(basis) = nb_per_i
         l_cshs(basis) = l_csh
         k_cshs(basis) = k_csh
         l_cats(basis) = l_cat
         k_cats(basis) = k_cat
      enddo

      end


      subroutine cneo_schwarz_set(basisin)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "cschwarz.fh"
#include "cschwarz_cneo.fh"
#include "basP.fh"

      integer basisin, basis
      basis = basisin + Basis_Handle_Offset

      sch_max = sch_maxs(basis)
      nb_per_i = nb_per_is(basis)
      l_csh = l_cshs(basis)
      k_csh = k_cshs(basis)
      l_cat = l_cats(basis)
      k_cat = k_cats(basis)

      end

      subroutine cneo_schwarz_tidy(nbas,bases)
C$Id$
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "cschwarz_cneo.fh"
#include "mafdecls.fh"
#include "basP.fh"
      integer nbas
      integer bases(nbas)

      integer i

      do i=1,nbas
         call cneo_schwarz_set(bases(i))
         call schwarz_tidy
      enddo
      
      end