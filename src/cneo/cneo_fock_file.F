      subroutine fock_en_rep_from_file(geom,bas1,bas2,
     $		nbf1, nbf2, jfac, tol2e, dens, fock, inttype,
     $      idx1, idx2)
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "cfock.fh"
#include "global.fh"
      integer geom, bas1, bas2
      integer nbf1, nbf2
      double precision jfac, tol2e
      double precision dens(nbf2*nbf2)
      double precision fock(nbf1*nbf1)
      character*(*) inttype
      integer idx1, idx2

      double precision den_tol, denmax
      integer ilo, jlo, klo, llo
      integer ihi, jhi, khi, lhi, i, j
c
      logical int2e_get_bf_range, int2e_file_read
      external int2e_get_bf_range, int2e_file_read
      integer idamax
      external idamax
      logical  oreverse

      oreverse = ((inttype .eq. 'neints') .or. 
     $      ((inttype .eq. 'nnints') .and. (idx1 .gt. idx2)))
      denmax = 0.0d0
      i = idamax(nbf2*nbf2, dens,1)
      denmax = max(denmax, abs(dens(i)))
      if(denmax.lt.1d-12) return
      den_tol = min(tol2e/denmax,tol2e/denmax**2) ! Threshold to screen integs only
c
      if (ga_nodeid().eq.0 .and. util_print('fockfile',print_debug))
     $     write(6,*) 'fockfile: tols ',tol2e, denmax,den_tol

 10   if (int2e_get_bf_range(ilo,ihi,jlo,jhi,klo,khi,llo,lhi)) then
         if (ilo.ne.1 .or. jlo.ne.1 .or. klo.ne.1 .or. llo.ne.1)
     $        call errquit('not a replicated file?',0, UNKNOWN_ERR)
         if (.not. oreverse) then
            if (ihi.ne.nbf1 .or. jhi.ne.nbf1 .or. 
     $                  khi.ne.nbf2.or. lhi.ne.nbf2)
     $         call errquit('not a replicated file?',0, UNKNOWN_ERR)
         else
            if (ihi.ne.nbf2 .or. jhi.ne.nbf2 .or. 
     $                  khi.ne.nbf1.or. lhi.ne.nbf1)
     $         call errquit('not a replicated file?',0, UNKNOWN_ERR)
         end if
c
         call inten_file_rep_fock(nbf1, nbf2,
     $        jfac, den_tol, dens, fock, inttype,oreverse)
c
         goto 10
c     
      end if

      end

      subroutine inten_file_rep_fock(nbf1, nbf2,
     $      jfac, tol2e, dens, fock, inttype, oreverse)
      implicit none
#include "errquit.fh"
#include "cint2efile.fh"
      integer nbf1, nbf2
      double precision tol2e, jfac(*)
      double precision dens(nbf2*nbf2), fock(nbf1*nbf1)
      character*(*) inttype
c
      integer neri
      logical int2e_buf_read
      external int2e_buf_read
      logical oreverse
c
c     While (integrals available) read buffers
c
 20   if (next_value .gt. n_in_rec) then
         if (.not. int2e_buf_read()) goto 1000 ! EOF
      end if
      if (next_value .gt. n_in_rec) call errquit('i2effb: uh?', 0,
     &       UNKNOWN_ERR)
      if (nleft_in_range .le. 0) goto 1000 ! end of range
c
      neri = nleft_in_range
c
      call fock_en_rep_mod_1_label(nbf1, nbf2, jfac, tol2e,
     $     neri, labels(1,next_value), values(next_value), 
     $     dens, fock, oreverse)
c
      nleft_in_range = nleft_in_range - neri
      next_value = next_value + neri
c     
      goto 20
c
 1000 continue
c
      end

      subroutine fock_en_rep_mod_1_label(nbf1, nbf2, jfac, tol2e,
     $     neri, labels, eri, dens, fock, oreverse)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"

      integer nbf1, nbf2
      double precision jfac, tol2e
      integer neri
      integer labels(4,neri)
      double precision eri(neri)
      double precision fock(*), dens(*)
      logical oreverse

      integer i, j, k, l, ind
      integer ij_fock, ik, jk, lk_dens, li, lj
      double precision g, gk, gj, fij, fkl, fik, fil, fjl, fjk

      integer maxim
      parameter (maxim = nw_max_nbf)
      integer im1(1:maxim), im2(1:maxim)

      do i=1, nbf1
         im1(i) = (i-1)*nbf1
      enddo

      do i=1, nbf2
         im2(i) = (i-1)*nbf2
      enddo

      do ind = 1, neri
         g = eri(ind)
         if (abs(g) .gt. tol2e) then
            gj = g*jfac

            if (.not. oreverse) then
               i = labels(1,ind)+1
               j = labels(2,ind)+1
               k = labels(3,ind)+1
               l = labels(4,ind)+1
            else
               i = labels(3,ind)+1
               j = labels(4,ind)+1
               k = labels(1,ind)+1
               l = labels(2,ind)+1
            end if
c
            ij_fock = i + im1(j)
            lk_dens = l + im2(k)
c
            fock(ij_fock) = fock(ij_fock) + gj*dens(lk_dens)
         end if
      end do
c     
      end