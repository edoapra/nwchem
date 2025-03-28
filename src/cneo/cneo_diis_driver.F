      Subroutine cneo_diis_driver(toll_s, derr, ipol, icall, MXhist, 
     &                 nbf_e, nbf_n, geom, AO_bas_han, g_focks, g_nfock, 
     &            g_dens,g_ndens, g_tmp_n, g_tmp, g_svecs, svals, 
     &            g_svecs_n,svals_n, diis, nodiis,
     .     g_diis,keeps,g_s12m,g_s12p,reset,nconstr, mult, mult_diis)
c
C$Id$
c
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "ccneo.fh"
#include "ccneo_constrain.fh"
c      
      integer ispin
      integer ipol
      integer MXhist
      integer nbf_n, nbf_e, nbf
      integer g_focks(2), g_nfock(max_quantum_nuc)
      integer g_dens(2), g_ndens(max_quantum_nuc)
      integer g_tmp, g_tmp_n, g_tmp_local
      integer g_svecs, g_svecs_n(max_quantum_nuc)
      integer nconstr
      double precision mult(nconstr),mult_diis(nconstr,mxhist)
      integer g_diis(mxhist)
      double precision svals(*), svals_n(*)
      integer nhist(2), icall(*)
      integer geom, AO_bas_han
      ! double precision derr(2+max_quantum_nuc) ! Diis error norm [output]
      double precision derr(*)
      double precision toll_s
      integer lpulay, npulay, ipulay
      logical reset
      integer g_focks_en(2+max_quantum_nuc)
      integer g_dens_en(2+max_quantum_nuc)
      integer g_tmps(2+max_quantum_nuc)
      integer nbfs(2+max_quantum_nuc)
      integer svals_start

      Save lpulay, npulay, ipulay, nhist
c
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "dftps.fh"
c
c
      integer g_s12m, g_s12p, me, nproc, g_tmp2,g_tmp1
      integer lfwght, ifwght
c     
      logical diis, nodiis,keeps
      double precision avail, dw_agg_avail
      integer n_n2_objs,i, inuc
      double precision one, zero
      parameter (one = 1.d0, zero = 0.d0)

c     Preliminaries
      call ga_zero(g_tmp)
      call ga_zero(g_tmp_n)
c
c
c     Use the DIIS scheme to accelerate convergence.

      do i=1,ipol
         nbfs(i) = nbf_e
         g_tmps(i) = g_tmp
         g_focks_en(i) = g_focks(i)
         g_dens_en(i) = g_dens(i)
      enddo
      do i=ipol+1,ipol+nuc_num
         nbfs(i) = nbf_n
         g_tmps(i) = g_tmp_n
         g_focks_en(i) = g_nfock(i-ipol)
         g_dens_en(i) = g_ndens(i-ipol)
      enddo
c     
      if(reset) then
         if( icall(1).eq.0 ) return
         icall(1) = -1
         icall(2) = -1
         do i=1,nuc_num
            icall(ipol+i) = -1
         enddo
      endif
      if (icall(1).lt.0 )then
         if (.not.MA_Free_Heap(lpulay))
     &      call errquit('diis_driver: cannot free heap',0, MA_ERR)
         do i=1,2*mxhist*(ipol+nuc_num)
            if (.not. ga_destroy(g_diis(i))) call errquit
     &           ('diis_driver: could not destroy g_diis', i, GA_ERR)
         enddo
         if(reset) then
            icall(1) = 0
            icall(2) = 0
            do i=1,nuc_num
               icall(ipol+i) = 0
            enddo
         endif
         return
      endif
      if (odftps) call pstat_on(ps_diis)
      me = ga_nodeid()
      nproc = ga_nnodes()
c
c     Allocate scratch space.
c
      if ( icall(1).eq.0 )then
c
c       Create necessary DIIS arrays.
c
        npulay  = (MXhist+1)*(MXhist+2)/2
        if (.not.MA_Alloc_Get(MT_Dbl,npulay,'Pulay',lpulay,ipulay))
     &     call errquit('diis_driver: cannot allocate Pulay',0, MA_ERR)
        call dfill(npulay, 0.0d0, dbl_mb(ipulay), 1)
c
c       Inquire about amount of global memory left.  DIIS uses
c       several n^2 arrays to store fock and error matrices and
c       needs probably at least 3 of each to be effective.
c
        avail = ga_memory_avail()
        call ga_igop(msg_min_ga_avail, avail, 1, 'min')
c
c        - subtract off a couple of n^2 objects for later use,
c          if not(keeps) we use 2 more matrices for S^-1/2 and S^+1/2
c
        dw_agg_avail = dble(nproc)*dble(avail)/8.0d0
     $                                    -2.0d0*dble(nbf_e)**2
        if (.not.keeps) dw_agg_avail = dw_agg_avail
     $                     -2.0d0*dble(nbf_e)**2-2.0d0*dble(nbf_n)**2
c
        n_n2_objs = nint(dw_agg_avail/dble(nbf_e)**2)
        if (2*mxhist*ipol.gt.n_n2_objs)then
c
c          Reset mxhist.
c          
           mxhist = n_n2_objs/(2*ipol)
           if (mxhist.lt.3)then
              diis = .false.
              nodiis = .true.
              if (me.eq.0)write(LuOut,2118)
              if (odftps) call pstat_off(ps_diis)
              return
           else
              if (me.eq.0)write(LuOut,2119)mxhist
           endif
         endif
*ga:1:0
         do i=1,mxhist*ipol
            if (.not. ga_duplicate(g_focks(1), g_diis(i), 'g_diis'))
     &           call errquit('diis_driver: error creating ga',0,GA_ERR)
         enddo
         do i=mxhist*ipol+1,mxhist*(ipol+nuc_num)
            if (.not. ga_duplicate(g_nfock(1), g_diis(i), 'g_diis'))
     &           call errquit('diis_driver: error creating ga',0,GA_ERR)
         enddo
         do i=mxhist*(ipol+nuc_num)+1,mxhist*(ipol*2+nuc_num)
            if (.not. ga_duplicate(g_focks(1), g_diis(i), 'g_diis'))
     &           call errquit('diis_driver: error creating ga',0,GA_ERR)
         enddo
         ! return
         do i=mxhist*(ipol*2+nuc_num)+1,2*mxhist*(ipol+nuc_num)
            if (.not. ga_duplicate(g_nfock(1), g_diis(i), 'g_diis'))
     &           call errquit('diis_driver: error creating ga',0,GA_ERR)
         enddo
      endif
 2118 format(/,10x,'Not enough global memory currently available for ',
     &       /,10x,'DIIS.  DIIS option disabled. ')
 2119 format(/,10x,'The amount of global memory currently available',
     &       /,10x,'for DIIS only allows: ',i2,' Fock matrices to be',
     &       /,10x,'stored and used in the extrapolation.')


      if (.not.MA_Push_Get(MT_Dbl,MXhist+1,'Fock wghts',lfwght,ifwght))
     &   call errquit('diis_driver: cannot allocate Fock wghts',0,
     &       MA_ERR)
      do ispin=1,ipol+nuc_num
      call dfill(MXhist+1, 0.0d0, dbl_mb(ifwght), 1)
c
c     Increment call counter.
c
      icall(ispin) = icall(ispin) + 1
      i = mxhist*(ispin-1)+1
c
c     Construct a temporary global array.
c       
      if (.not. ga_duplicate(g_diis(i), g_tmp2, 'gtmp2'))
     &     call errquit('diis_driver: error creating ga',0,GA_ERR)
      call ga_zero(g_tmp2)

      if (.not. ga_duplicate(g_diis(i), g_tmp1, 'gtmp1'))
     &     call errquit('diis_driver: error creating ga',0,GA_ERR)
      call ga_zero(g_tmp1)

      g_tmp_local = g_tmps(ispin)

      if (odftps) call pstat_on(ps_bld12)
      if(.not.keeps) then
c
c     get S^(-1/2)
c
         if (.not. ga_duplicate(g_diis(i), g_s12m, 'gs12m'))
     &        call errquit('diis_driver: error creating ga',0,GA_ERR)
c
c     get S^(+1/2)
c
         if (.not. ga_duplicate(g_diis(i), g_s12p, 'gs12p'))
     &        call errquit('diis_driver: error creating ga',0,GA_ERR)

         if (ispin .gt. ipol) then
            inuc = ispin-ipol
            call cneo_set_cdft('n')
            svals_start = ((ispin-ipol)-1)*nbf_n+1
            call ga_zero(g_s12m)
            call ga_zero(g_s12p)
            call diis_bld12(1.0d-9,svals_n(svals_start),
     $                       g_svecs_n(inuc),g_s12m,g_tmp_local,2)
            call diis_bld12(1.0d-9,svals_n(svals_start),
     $                       g_svecs_n(inuc), g_s12p,g_tmp2,3)
            call cneo_set_cdft('e')
         else
            call diis_bld12(toll_s, svals, g_svecs, g_s12m, 
     $                                                g_tmp_local,2)
            call diis_bld12(toll_s, svals, g_svecs, g_s12p,g_tmp2,3)
         endif
      endif
      if (odftps) call pstat_off(ps_bld12)
c
c     Mat mult: F * S^(-1/2)
c
      nbf = nbfs(ispin)
      if (odftps) call pstat_on(ps_dgemm)
      call ga_dgemm('N','N',nbf,nbf,nbf,one,g_focks_en(ispin),
     .     g_s12m,zero,g_tmp1)
      if (odftps) call pstat_off(ps_dgemm)
      if(.not.keeps) then
         if (.not. ga_destroy(g_s12m)) call errquit
     &        ('diis_driver: could not destroy g_s12m', 0, GA_ERR)
      endif

c
c     Mat mult: S^(+1/2) * D
c
      if (odftps) call pstat_on(ps_dgemm)
      call ga_dgemm('N','N',nbf,nbf,nbf,one,g_s12p,g_dens_en(ispin),
     ,     zero,g_tmp2)
      if (odftps) call pstat_off(ps_dgemm)
      if(.not.keeps) then
         if (.not. ga_destroy(g_s12p)) call errquit
     &        ('diis_driver: could not destroy g_s12p', 0, GA_ERR)
      endif
c
C     Mat mult:  (S^(+1/2) * D) * (F * S^(-1/2))
c
      if (odftps) call pstat_on(ps_dgemm)
      call ga_zero(g_tmp_local)
      call ga_dgemm('N','N',nbf,nbf,nbf,one,g_tmp2,g_tmp1,zero,
     $                                             g_tmp_local)
      if (odftps) call pstat_off(ps_dgemm)
      if (.not. ga_destroy(g_tmp1)) call errquit
     &     ('diis_driver: could not destroy g_tmp1', 0, GA_ERR)

      
c
c     Assemble: S^(-1/2)*F*D*S^(1/2) - S^(1/2)*D*F*S^(-1/2)
c
      call ga_transpose(g_tmp_local,g_tmp2)
      call ga_add(-1d0,g_tmp_local,1d0,g_tmp2,g_tmp_local)

      if (.not. ga_destroy(g_tmp2)) call errquit
     &     ('diis_driver: could not destroy g_tmp2', 0, GA_ERR)
c
c
c     Compute Pulays B matrix.
      derr(ispin) = ga_ddot(g_tmp_local,g_tmp_local)
c
      call diis_core(derr(ispin), Dbl_MB(ipulay), MXhist, nhist,
     &     icall(ispin), ispin, ipol+nuc_num, g_diis, g_tmp_local)

c
      enddo

c
c     Solve DIIS equation for Hamiltonian weights.
c
      call diis_bsolve(Dbl_MB(ipulay), Dbl_MB(ifwght),
     &                 MXhist, nhist(1))
c
c     Compute an extrapolated Hamiltonian using the weights.
c
      call diis_hamwgt(Dbl_MB(ifwght), MXhist, nhist, icall, 
     &            ipol+nuc_num,g_focks_en, g_diis, nconstr,
     $            mult,mult_diis)


      do i=1,nuc_num
         call ga_dadd(1.0d0,g_nfock(i),-1.0d0,g_fock_add(i),g_fock0(i))
         call cneo_solve_constrain(i)
         call ga_dadd(1.0d0,g_fock0(i),1.0d0,g_fock_add(i),g_nfock(i))
      enddo
c
      if (.not.ma_pop_stack(lfwght))
     &   call errquit('diis_driver: cannot pop stack',0, MA_ERR)
c
      if (odftps) call pstat_off(ps_diis)

c     Clean up and leave
      call ga_zero(g_tmp)
      call ga_zero(g_tmp_n)

      return
      end



      subroutine cneo_set_cdft(particle)
      implicit none 
#include "errquit.fh"
#include "cdft.fh"
#include "ccneo.fh"
#include "ccneo_dft.fh"
      character particle

      if (particle .eq. 'e') then
         nbf_ao = nbf_e 
      else if (particle .eq. 'n') then
         nbf_ao = nbf_n
      else 
         call errquit('cneo_set_cdft:invalid particle type',0,0)
      endif

      end



