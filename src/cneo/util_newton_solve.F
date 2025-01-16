      subroutine build_Jacobian(g_x,xlen,g_J,residual,delta,oreverse)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"

      integer g_x, g_J, xlen
      external residual
      double precision delta
      logical oreverse

      double precision acc
      integer typex, xrow, xcol
      integer g_y, g_y1, g_x1
      double precision x(xlen), y(xlen), J(xlen,xlen)
      double precision x1(xlen),y1(xlen)
      integer i,ii
      logical status
      double precision dx

      call ga_get(g_x, 1,xlen,1,1, x,1)
      if (.not. ga_duplicate(g_x, g_y, 'tmp:y')) call errquit
     $      ('build:ga_duplicate?',0,0)
      if (.not. ga_duplicate(g_x, g_x1, 'tmp:x1')) call errquit
     $      ('build:ga_duplicate?',0,0)
      if (.not. ga_duplicate(g_x, g_y1, 'tmp:y1')) call errquit
     $      ('build:ga_duplicate?',0,0)

      call residual(acc,g_x,g_y)
      call ga_get(g_y,1,xlen,1,1,y,1)

      do i=1,xlen
         x1 = x
         dx = sign(delta,y(i))
         if (oreverse) dx = -dx
         x1(i) = x1(i) + dx
         call ga_zero(g_x1)
         call ga_put(g_x1,1,xlen,1,1,x1,1)
         call ga_zero(g_y1)
         call residual(acc,g_x1,g_y1)
         call ga_get(g_y1,1,xlen,1,1,y1,1)
         do ii=1,xlen
            J(ii,i) = (y1(ii)-y(ii))/dx
         enddo
      enddo

      call ga_put(g_J,1,xlen,1,xlen,J,xlen)

      status = .true.
      status = status .and. ga_destroy(g_y)
      status = status .and. ga_destroy(g_y1)
      status = status .and. ga_destroy(g_x1)
      if (.not. status) call errquit('build_jacobian:ga_destroy?',0,0)


      end

      subroutine ga_inverse_3(g_a, g_s, det,osvd)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

      integer g_a, g_s
      double precision matrix(3,3), inv_matrix(3,3)
      double precision det
      logical osvd
      integer i,j

      call ga_get(g_a,1,3,1,3,matrix,3)
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
         do i=1,3
            do j=1,3
               inv_matrix(i,j) = inv_matrix(i,j)/det
            enddo
         enddo
         call ga_put(g_s, 1,3,1,3,inv_matrix,3)
      else if (osvd) then
         write(6,*) 'Calling pseudoinverse...'
         call pseudoinverse(g_a,g_s)
      else
         write(6,*) 'small determinant:',det
         call ga_put(g_s, 1,3,1,3,inv_matrix,3)
      endif
      end


c     Modified from Zehua Chen's implementation on CP2K    
c     https://github.com/zc62/cp2k/blob/cneo/src/qs_cneo_methods.F  

      subroutine util_Newton_solve(g_x,residual,
     $            tol,maxiter,converged)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

      integer g_x
      integer xlen
      external residual
      integer maxiter
      double precision tol

      integer typex,xrow,xcol
      integer g_y,g_y1,g_J,g_inv,g_J1
      integer g_dx, g_x1
      double precision dx_norm,step
      double precision g0,g1,g0p

      double precision acc, rmax, det, rmax1
      integer iter,i
      logical converged, status

      integer g_list(15),ilist

      double precision one,mone,zero
      parameter (one=1.0d0)
      parameter (mone=-1.0d0)
      parameter (zero=0.0d0)

      ilist=1

      call ga_inquire(g_x,typex,xrow,xcol)
      if ((xcol.ne.1)) call errquit
     $      ('Newton_solve:input g_x shape not supported',0,0)
      xlen = xrow

      if(.not. ga_create(mt_dbl,xlen,xlen,'Newton:J',32,32,g_J))
     $      call errquit('newton_solve:ga_create?',0,0)
      call ga_zero(g_J)
      g_list(ilist) = g_J
      ilist = ilist+1
      if(.not. ga_duplicate(g_J,g_J1,'Newton:J_tmp')) call errquit
     $      ('Newton_solve:ga_duplicate?',0,0)
      call ga_zero(g_J1)
      g_list(ilist) = g_J1
      ilist = ilist+1

      if(.not. ga_duplicate(g_J,g_inv,'Newton:J_inv')) call errquit
     $      ('Newton_solve:ga_duplicate?',0,0)
      call ga_zero(g_inv)
      g_list(ilist) = g_inv
      ilist = ilist+1

      if(.not. ga_duplicate(g_x,g_y,'Newton:y')) call errquit
     $      ('Newton_solve:ga_duplicate?',0,0)
      call ga_zero(g_y)
      g_list(ilist) = g_y
      ilist = ilist+1

      if(.not. ga_duplicate(g_x,g_x1,'Newton:y')) call errquit
     $      ('Newton_solve:ga_duplicate?',0,0)
      call ga_zero(g_x1)
      g_list(ilist) = g_x1
      ilist = ilist+1

      if(.not. ga_duplicate(g_x,g_y1,'Newton:y1')) call errquit
     $      ('Newton_solve:ga_duplicate?',0,0)
      call ga_zero(g_y1)
      g_list(ilist) = g_y1
      ilist = ilist+1

      if(.not. ga_duplicate(g_x,g_dx,'Newton:y1')) call errquit
     $      ('Newton_solve:ga_duplicate?',0,0)
      call ga_zero(g_dx)
      g_list(ilist) = g_dx

      do iter=1,maxiter
         call residual(acc,g_x,g_y)
         call ga_maxelt(g_y,rmax)
         if (rmax .le. tol) then
            converged = .true.
            goto 999
         endif
         call ga_zero(g_J)
         call ga_zero(g_inv)
         if (iter.eq.1) then
            ! try if zero guess is better
            call ga_zero(g_x1)
            call residual(acc,g_x1,g_y1)
            call ga_maxelt(g_y1,rmax1)
            if (rmax1 .lt. rmax) then
               call ga_copy(g_x1,g_x)
               call ga_copy(g_y1,g_y)
            endif
         endif

         call build_Jacobian(g_x,xlen,g_J,residual,1.0d-3,.false.)
         if (xlen.ne.3) call errquit
     $               ('newton_solve:not implemented',0,0)
         call ga_inverse_3(g_J,g_inv,det,.false.)
         if (abs(det).lt.1.0d-8) then
            call build_Jacobian(g_x,xlen,g_J1,residual,1.0d-3,.true.)
            call ga_dadd(one,g_J,one,g_J1,g_J)
            call ga_dscal(g_J,5.0d-1)
            call ga_inverse_3(g_J,g_inv,det,.true.)
         endif

         call ga_zero(g_dx)

         call ga_dgemm('n','n',xlen,1,xlen,mone,g_inv,g_y,zero,g_dx)
         dx_norm = sqrt(ga_ddot(g_dx,g_dx))
         g0 = sqrt(ga_ddot(g_y,g_y))
         call ga_copy(g_x,g_x1)
         call ga_copy(g_y,g_y1)
         call ga_dadd(one,g_x,one,g_dx,g_x)
         call residual(acc,g_x,g_y)
         g1 = sqrt(ga_ddot(g_y,g_y))

         ! line search
         step = one

 10      if (g1 .ge. g0) then
            if (step.lt.1.0d-2) then
               goto 99
            endif
            g0p = -g0/(step*dx_norm)
            step = step*max(-g0p/(2.0d0*(g1-g0-g0p)),1.0d-1)
            call ga_dadd(one,g_x1,step,g_dx,g_x)
            call residual(acc,g_x,g_y)
            g1 = sqrt(ga_ddot(g_y,g_y))
            goto 10
         endif

 99      continue

      enddo

 999  continue

      if (.not.converged) then
         write(6,*) 'Warning:newton_solve failed to converge'
      endif

      status = .true.
      do i=ilist,1,-1
         status = status .and. ga_destroy(g_list(i))
      enddo
      if (.not.status) call errquit('newton_solve:ga_destroy?',0,0)

      end

      subroutine pseudoinverse(g_a,g_s)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer g_a    ![input]
      integer g_s    ![output]

      integer type,m,n,nsing
      integer g_u,g_vt,lval,kval
      integer g_test,g_v,g_inv
      integer g_tmp
      logical status

c
c     Function not tested yet: currently didn't see pseudoinverse being called
c

      call ga_inquire(g_a,type,m,n)
      nsing = min(n,m)
      if (.not. ga_create(mt_dbl,m,nsing,'pseudoinverse:u',32,32,g_u))
     $      call errquit('ga_create?',0,0)
      if (.not. ga_create(mt_dbl,m,m,'pseudoinverse:test',32,32,g_test))
     $      call errquit('ga_create?',0,0)
      if (.not. ga_create(mt_dbl,nsing,n,'pseudoinverse:vt',32,32,g_vt))
     $      call errquit('ga_create?',0,0)
      if (.not. ga_create(mt_dbl,nsing,nsing,'pseudoinverse:val',
     $      32,32,g_inv)) call errquit('ga_create?',0,0)
      if (.not.ma_push_get(mt_dbl,nsing,'pseudoinverse:vals',
     $      lval,kval)) call errquit('ma_push_get?',0,0)

      call ga_svd_seq(g_a,g_u,g_vt,dbl_mb(kval))

      call gen_reciprocal(g_inv,dbl_mb(kval),nsing)

      if (.not.ga_duplicate(g_v,g_tmp,'pseudoinverse:tmp')) call 
     $      errquit('ga_duplicate?',0,0)
      call ga_dgemm('t','n',n,nsing,nsing,1.0d0,g_vt,g_inv,0d0,g_tmp)
      call ga_dgemm('n','t',n,m,nsing,1.0d0,g_tmp,g_u,0d0,g_s)

      call ga_dgemm('n','n',m,m,n,1.0d0,g_a,g_s,0d0,g_test)
      print *, 'debug:pseudoinverse:test:'
      call ga_print(g_test)

      status = .true.
      status = status .and. ga_destroy(g_u)
      status = status .and. ga_destroy(g_test)
      status = status .and. ga_destroy(g_vt)
      status = status .and. ga_destroy(g_inv)
      status = status .and. ma_pop_stack(lval)
      if (.not. status) call errquit
     $      ('pseudoinverse:deallocate?',0,0)

      end 

      subroutine gen_reciprocal(g_a,vals,nelem)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

      integer g_a
      double precision vals(*)
      integer nelem
      double precision reciprocal(nelem,nelem)

      integer i
      reciprocal = 0d0
      do i=1,nelem
         reciprocal(i,i) = 1.0d0/vals(i)
      enddo

      call ga_put(g_a,1,nelem,1,nelem,reciprocal,nelem)

      end


      
      