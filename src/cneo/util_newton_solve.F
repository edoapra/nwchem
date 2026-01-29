      subroutine build_Jacobian(x, xlen, J, residual, delta, oreverse)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"

      integer xlen
      double precision x(xlen)
      external residual
      double precision delta
      logical oreverse

      double precision y(xlen), J(xlen,xlen)
      double precision x1(xlen),y1(xlen)
      integer i,ii
      double precision dx

      call residual(x, y, xlen)
      do i=1,xlen
         x1(:) = x(:)
         dx = sign(delta, y(i))
         if (oreverse) dx = -dx
         x1(i) = x1(i) + dx
         call residual(x1, y1, xlen)
         J(:,i) = (y1(:) - y(:))/dx
      enddo

      end

      subroutine inverse_3(matrix, inv_matrix, det, osvd)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"

      double precision matrix(3,3), inv_matrix(3,3)
      double precision det
      logical osvd
      integer i,j

      det = matrix(1, 1) * (matrix(2, 2) * matrix(3, 3) - 
     $                           matrix(2, 3) * matrix(3, 2))
     $    - matrix(1, 2) * (matrix(2, 1) * matrix(3, 3) - 
     $                           matrix(2, 3) * matrix(3, 1))
     $    + matrix(1, 3) * (matrix(2, 1) * matrix(3, 2) - 
     $                           matrix(2, 2) * matrix(3, 1))

      inv_matrix(1, 1) = matrix(2, 2) * matrix(3, 3) - 
     $                           matrix(2, 3) * matrix(3, 2)
      inv_matrix(1, 2) = matrix(1, 3) * matrix(3, 2) - 
     $                           matrix(1, 2) * matrix(3, 3)
      inv_matrix(1, 3) = matrix(1, 2) * matrix(2, 3) - 
     $                           matrix(1, 3) * matrix(2, 2)
      inv_matrix(2, 1) = matrix(2, 3) * matrix(3, 1) - 
     $                           matrix(2, 1) * matrix(3, 3)
      inv_matrix(2, 2) = matrix(1, 1) * matrix(3, 3) - 
     $                           matrix(1, 3) * matrix(3, 1)
      inv_matrix(2, 3) = matrix(1, 3) * matrix(2, 1) - 
     $                           matrix(1, 1) * matrix(2, 3)
      inv_matrix(3, 1) = matrix(2, 1) * matrix(3, 2) - 
     $                           matrix(2, 2) * matrix(3, 1)
      inv_matrix(3, 2) = matrix(1, 2) * matrix(3, 1) - 
     $                           matrix(1, 1) * matrix(3, 2)
      inv_matrix(3, 3) = matrix(1, 1) * matrix(2, 2) - 
     $                           matrix(1, 2) * matrix(2, 1)

      if (abs(det).ge.1.0d-8)  then
         inv_matrix(:,:) = inv_matrix(:,:)/det
      else if (osvd) then
         write(luout, *) 'Calling pseudoinverse...'
         call pseudoinverse(matrix, inv_matrix)
      else
         write(luout,*) 'small determinant:',det
      endif
      end


c     Modified from Zehua Chen's implementation on CP2K    
c     https://github.com/zc62/cp2k/blob/cneo/src/qs_cneo_methods.F  

      subroutine util_Newton_solve(x, xlen, residual, tol, maxiter, 
     $                             converged)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"

      integer xlen, maxiter
      logical converged
      double precision x(xlen)

      external residual
      double precision tol

      integer typex, xrow, xcol
      double precision dx_norm,step
      double precision g0,g1,g0p

      double precision rmax, det, rmax1
      integer iter,i,ii
      logical status

      double precision J(xlen, xlen), J1(xlen, xlen), Jinv(xlen, xlen) 
      double precision y(xlen), y1(xlen), x1(xlen), dx(xlen)
      double precision one,mone,zero
      parameter (one=1.0d0)
      parameter (mone=-1.0d0)
      parameter (zero=0.0d0)

      if (xlen.ne.3) call errquit('newton_solve:not implemented',0,0)

      ! Initialize
      J(:,:) = 0d0
      J1(:,:) = 0d0
      Jinv(:,:) = 0d0
      y(:) = 0d0
      y1(:) = 0d0
      x1(:) = 0d0
      dx(:) = 0d0
      
      do iter=1,maxiter
         call residual(x, y, xlen)
         rmax = maxval(abs(y(:)))
         if (rmax .le. tol) then
            converged = .true.
            goto 999
         endif

c         if (iter.eq.1) then
c            ! try if zero guess is better
c            call residual(x1, y1, xlen)
c            rmax1 = maxval(abs(y1(:)))
c            if (rmax1 .lt. rmax) then
c               x(:) = x1(:)
c               y(:) = y1(:)
c               rmax = rmax1
c            endif
c         endif

         call build_Jacobian(x, xlen, J, residual, 1.0d-3, .false.)
         call inverse_3(J, Jinv, det, .false.)
         if (abs(det).lt.1.0d-8) then
            call build_Jacobian(x, xlen, J1, residual, 1.0d-3, .true.)
            J(:,:) = 0.5d0*(J1(:,:) + J(:,:))
            call inverse_3(J, Jinv, det, .true.)
         endif

         dx(:) = 0d0
         do ii=1,xlen
            dx(:) = dx(:) - Jinv(:,ii) * y(ii)
         enddo
         dx_norm = dsqrt(dot_product(dx, dx))
         g0 = dsqrt(dot_product(y, y))

         x1(:) = x(:)
         y1(:) = y(:)
         x(:) = x(:) + dx(:)

         call residual(x, y, xlen)
         g1 = dsqrt(dot_product(y, y))

         ! line search
         step = one
 10      if (g1 .ge. g0) then
            if (step.lt.1.0d-2) then
               goto 99
            endif
            g0p = -g0/(step*dx_norm)
            step = step*max(-g0p/(2.0d0*(g1-g0-g0p)),1.0d-1)
            x(:) = x1(:) + step*dx(:)
            call residual(x, y, xlen)
            g1 = dsqrt(dot_product(y,y))
            goto 10
         endif

 99      continue

      enddo

 999  continue

      if (.not.converged) then
         write(luout,*) 'Warning:newton_solve failed to converge'
      endif

      end

      subroutine pseudoinverse(matrix, inv_matrix)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
      double precision matrix(3,3), inv_matrix(3,3)


      integer type,m,n,nsing
      integer g_u,g_vt,lval,kval
      integer g_test,g_v,g_inv
      integer g_tmp
      logical status
      double precision evals(3)
      double precision work1(3,3), work2(3,3)
      integer info, ival
      double precision threshold
      parameter threshold = 1d-8

c
c     Function not tested yet: currently didn't see pseudoinverse being called
c
      work1(:,:) = matrix(:,:)
      call dsyev('v', 'u', 3, work1, 3, evals, work2, 9, info)
      if(info.ne.0)
     $    call errquit('pseudoinverse: dsyev failed',info,0) 

      do ival=1, 3
        if(abs(evals(ival)).gt.1d-8) then
          work2(:,ival) = work1(:,ival)/evals(ival)
        else
          work2(:,ival) = 0d0
        endif 
      enddo

      call dgemm('n','t',3,3,3,1d0,work1,3,work2,3,0d0,inv_matrix,3)
      call dgemm('n','n',3,3,3,1d0,matrix,3,inv_matrix,3,0d0,work1,3)
      write(luout,*) 'debug:pseudoinverse:test:'
      do ival=1, 3
        write(luout, "(3(F12.8))") work1(ival, :)
      enddo

      end 