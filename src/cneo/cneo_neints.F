      subroutine cneo_neints(geom, bas1, bas2, jfac, g_dens, g_fock)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "cfock.fh"
#include "util.fh"
#include "stdio.fh"

      integer geom, bas1, bas2
      integer g_dens
      integer g_fock
      double precision jfac

      integer cnmax
      parameter(cnmax = nw_max_shells)
      integer,dimension(cnmax)   :: cntoce1, cntoce2
      integer,dimension(2,cnmax) :: cntobfr1, cntobfr2

      logical status

      ! Load basis1 info (Fock matrix)
      status = bas_numbf(bas1, nbf1)
      status = status.and.bas_numcont(bas1, nsh1)
      if (.not.status)
     $  call errquit("cneo_neints: failed to load bas1", bas1, BASIS_ERR)

     ! Load basis2 info (Density matrix)
      status = status.and.bas_numbf(bas2, nbf2)
      status = status.and.bas_numcont(bas2, nsh2)
      if (.not.status)
     $  call errquit("cneo_neints: failed to load bas2", bas2, BASIS_ERR)

      status = status.and.geom_ncent(geom, natoms)

      ! Check available memory
      call intb_mem_2e4c(maxg, lenscr)
      mem_rep = nbf1*nbf1*2 + nbf2*nbf2*2 + nsh2*nsh2 + maxquartet
     $        + ma_sizeof(MT_INT,4*maxquartet+nbf1*(nbf1+1),MT_DBL)
     $        + nbf1*(nbf1+1)/2
     $        + ma_sizeof(MT_INT,2*tdim,MT_DBL) 
     $        + 8*ma_sizeof_overhead(MT_DBL)      

      ! Allocate Fock matrix
      if(.not.ma_push_get(mt_dbl, nbf1*nbf1, "fock", l_fock, k_fock))
     $  call errquit("cneo_neints: failed to allocate fock matrix", 
     $                nbf1, MA_ERR)
      call dcopy(nbf1*nbf1, 0d0, 0, dbl_mb(k_fock),1)

      ! Allocate Density matrix
      if(.not.ma_push_get(mt_dbl, nbf2*nbf2, "dens", l_dens, k_dens))
     $  call errquit("cneo_neints: failed to allocate dens matrix",
     $                nbf2, MA_ERR)  
      call ga_get(g_dens, 1, nbf2, 1, nbf2, dbl_mb(k_dens), nbf2) 

      ! Allocate reduced density matrix
      if(.not.ma_push_get(mt_dbl, nsh2*nsh2, "rdens", l_rdens, k_rdens))
     $  call errquit("cneo_neints: failed to allocate rdens",
     $                nsh2*nsh2, MA_ERR) 
      call dcopy(nsh2*nsh2, 0d0, 0, dbl_mb(k_rdens),1)

      ! Build maps and construct density over atoms
      do ish = 1, nsh1
         if(.not.bas_cn2bfr(bas1,ish,cntobfr1(1,ish),cntobfr1(2,ish)))
     $        call errquit('cneo_neints: bad basis', 0, BASIS_ERR)
         if (.not. bas_cn2ce(bas1, ish, cntoce1(ish)))
     $        call errquit('cneo_neints: bad basis', 0, BASIS_ERR)
      end do
      do ish = 1, nsh2
         if (.not. bas_cn2bfr(bas2,ish,cntobfr2(1,ish),cntobfr2(2,ish)))
     $        call errquit('cneo_neints: bad basis', 0, BASIS_ERR)
         if (.not. bas_cn2ce(bas2, ish, cntoce2(ish)))
     $        call errquit('cneo_neints: bad basis', 0, BASIS_ERR)
      end do
      call util_mat_reduce(nbf2, nsh2, cntobfr2, dbl_mb(k_dens), 1, 
     $                      dbl_mb(k_rdens), 'rms')
     

      end