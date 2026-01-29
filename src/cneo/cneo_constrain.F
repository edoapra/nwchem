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

      double precision threshold
      parameter( threshold = 1.0d-9 )

      integer maxiter
      parameter( maxiter = 10 )

      call pos_ana_set(idx, bas_n(idx))

      ! Allocate memory
      if(.not.ma_push_get(mt_dbl, 8*nbf_solve**2, 'cneo_solve: work',
     $                    l_work, k_work))
     $   call errquit("allocation failed", 0, MA_ERR)
      k_r1       = k_work
      k_r2       = k_r1 + nbf_solve**2
      k_r3       = k_r2 + nbf_solve**2
      k_fock0    = k_r3 + nbf_solve**2
      k_fock     = k_fock0 + nbf_solve**2
      k_mo       = k_fock + nbf_solve**2
      k_ortho    = k_mo + nbf_solve**2
      k_dens     = k_ortho + nbf_solve**2

      ! Prefetch necessary matrices
      call ga_get(g_int1e_r_solve(1),1,nbf_solve,1,nbf_solve,
     $            dbl_mb(k_r1),nbf_solve)
      call ga_get(g_int1e_r_solve(2),1,nbf_solve,1,nbf_solve,
     $            dbl_mb(k_r2),nbf_solve)
      call ga_get(g_int1e_r_solve(3),1,nbf_solve,1,nbf_solve,
     $            dbl_mb(k_r3),nbf_solve)
      call ga_get(g_fock0_solve,1,nbf_solve,1,nbf_solve,
     $            dbl_mb(k_fock0),nbf_solve)
      call ga_get(g_ortho_n(idx), 1, nbf_solve, 1, nbf_solve, 
     $            dbl_mb(k_ortho), nbf_solve)

 
      !f_cneo(1:3, idx) = 0d0
      call util_Newton_solve(f_cneo(1,idx), 3, position_analysis,
     $            threshold, maxiter, converged_solve)

      call cneo_fock_constrain(f_cneo(1,idx), dbl_mb(k_r1), 
     $                         dbl_mb(k_fock), nbf_solve)

      call ga_put(g_fock_add(idx),1,nbf_solve,1,nbf_solve,
     $            dbl_mb(k_fock), nbf_solve)

      if(.not.ma_chop_stack(l_work))
     $  call errquit('cneo_solve: deallocation failed', 0, MA_ERR)

      end
      
      
      subroutine position_analysis(f, r, flen)
      implicit none 
#include "errquit.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "ccneo.fh"
#include "ccneo_constrain.fh"

      integer flen
      double precision f(flen), r(flen)

      integer type, nrow, ncol
      double precision One
      integer g_mo, g_over_n
      double precision eval(nbf_solve)
      integer g_ndens, g_odens, g_rD
      double precision exp_pos_local(3)
      double precision vec1(nbf_solve), vec2(nbf_solve)
      integer i,j
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      parameter (One = 1.0d0)

      call cneo_fock_constrain(f, dbl_mb(k_r1),dbl_mb(k_fock),nbf_solve)
      call daxpy(nbf_solve**2, 1d0, dbl_mb(k_fock0),1,dbl_mb(k_fock),1)

      ! Diagonalization
      call cneo_diag(dbl_mb(k_fock), dbl_mb(k_ortho),dbl_mb(k_mo),eval, 
     $               nbf_solve)

      ! Get Density Matrix
      call dcopy(nbf_solve**2, 0d0, 0, dbl_mb(k_dens), 1)
      call dsyr('u', nbf_solve, 1d0, dbl_mb(k_mo), 1, dbl_mb(k_dens), 
     $           nbf_solve)
      call cneo_symm('u', dbl_mb(k_dens), nbf_solve)

      ! Get new position/residual
      call dgemv('t', nbf_solve**2, 3, 1d0, dbl_mb(k_r1), nbf_solve**2, 
     $            dbl_mb(k_dens),1,0d0,r,1)

      end

      subroutine pos_ana_set(idx, basis)
      implicit none 
#include "nwc_const.fh"
#include "ccneo.fh"
#include "ccneo_constrain.fh"
#include "errquit.fh"
#include "bas.fh"

      integer idx
      integer basis
      integer i

      do i=1,3
         g_int1e_r_solve(i) = g_int1e_r(i,idx)
      enddo
      g_fock0_solve = g_fock0(idx)
      if (.not. bas_numbf(basis, nbf_solve)) 
     $   call errquit('pos_ana_set:bas_n',0,0)
      bas_solve = basis

      end

      subroutine cneo_diag_wrapper(g_fock, g_ortho, g_vectors, work, 
     $                             values, n)
      implicit none
#include "global.fh"
      integer n, g_fock, g_ortho, g_vectors
      double precision work(n**2, 3), values(n)

      call ga_get(g_fock, 1, n, 1, n, work(1,1), n)
      call ga_get(g_ortho, 1, n, 1, n, work(1,2), n)
      call cneo_diag(work(1,1), work(1,2), work(1,3), values, n)
      call ga_put(g_vectors, 1, n, 1, n, work(1,3), n)
      end

      subroutine cneo_diag(matrix, metric, vectors, values, n)
      implicit none
      integer n 
      double precision, dimension(n,n) :: matrix, metric, vectors
      double precision, dimension(n) :: values
      integer idum, info

      ! Transform to standard eigenvalue problem
      call dgemm('n', 'n', n, n, n, 1d0, matrix, n, metric, n, 0d0, 
     $            vectors, n)
      call dgemm('t', 'n', n, n, n, 1d0, metric, n, vectors, n, 0d0,
     $            matrix, n)

      ! Diagonalization
      call dsyev('v', 'l', n, matrix, n, values, vectors, n*n, info)
      if (info.ne.0) 
     $   call errquit('cneo_diag: diagonalization failed', info, 0)

      ! Back-transformation
      call dgemm('n', 'n', n, n, n, 1d0, metric, n, matrix, n, 0d0,
     $            vectors, n)
      
      end

      subroutine cneo_fock_constrain(f, int1e_r, fock, n)
      implicit none
      integer n 
      double precision f(3), int1e_r(n**2,3), fock(n**2)

      fock(:) = f(1) * int1e_r(:,1) +
     $          f(2) * int1e_r(:,2) +
     $          f(3) * int1e_r(:,3)

      end

      subroutine cneo_blds12(g_evecs, svals, g_out, n, which)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
      integer n, g_evecs, g_out, which
      double precision svals(n)
      integer l_work1
      integer k_work1
      integer i, offset

      if(.not.ma_push_get(mt_dbl, n**2, 'work1', l_work1, k_work1))
     $  call errquit('cneo_blds12: allocation failed', 0, MA_ERR)

      call ga_get(g_evecs, 1, n, 1, n, dbl_mb(k_work1), n)
      offset = 0
      do i=1,n
        if(svals(i).lt.1.0d-9) then
          write(luout, *) 'Warning: small S eigenvalue: ', svals(i)
          call dscal(n, 0d0, dbl_mb(k_work1+offset), 1)
        else
          if(which.eq.1) then
            call dscal(n, 1d0/dsqrt(svals(i)),dbl_mb(k_work1+offset),1)
          elseif(which.eq.2) then
            call dscal(n, dsqrt(svals(i)), dbl_mb(k_work1+offset), 1)
          endif
        endif
        offset = offset + n
      enddo
      call ga_put(g_out, 1, n, 1, n, dbl_mb(k_work1), n)

      if(.not.ma_chop_stack(l_work1))
     $  call errquit('cneo_blds12: deallocation failed', 0, MA_ERR) 
      
      end

      subroutine cneo_symm(which, matrix, n)
      implicit none
      integer n
      character(len=1) which
      double precision matrix(n,n)

      integer i

      if(which.eq.'u') then
        do i=1,n-1
          matrix(i+1:n, i) = matrix(i, i+1:n)
        enddo
      else if (which.eq.'l') then
        do i=1,n-1 
          matrix(i, i+1:n) = matrix(i+1:n, i)
        enddo
      endif
      end




