      subroutine dft_alcv(e,deltava,deltavb,
     n     basis,nbf,nmo,nocc,oskel,
     i     cetobfr,
     c     vmo,vmo_a,
     b     g_movecs,eval)
      implicit none
#include "errquit.fh"
      double precision e ! [out]
      integer deltava,deltavb ! [in]
      integer nbf,nmo
      double precision vmo(nmo,*)
      double precision vmo_a(nmo,*)
      integer basis
      integer g_movecs
      double precision eval(*)
      integer cetobfr(2,*)
      integer nocc
      logical oskel
c
      logical dova,dovb
      integer atom
c      

      dova=.false.
      dovb=.false.
      write(6,*) ' in : deltava ',deltava
      write(6,*) ' in : deltavb ',deltavb
      write(6,*) ' in : nocc ',nocc
      write(6,*) ' ev 1    ',eval(1)
      write(6,*) ' ev 2    ',eval(2)
      write(6,*) ' ev 3    ',eval(3)
      write(6,*) ' ev homo ',eval(nocc)
      write(6,*) ' ev lumo ',eval(nocc+1)
      write(6,*) ' ev nmo  ',eval(nmo)
      if(deltava.ne.0.and.deltavb.eq.0) then
         dova=.true.
      elseif(deltava.ne.0.and.deltavb.ne.0) then
         if(deltava.ne.deltavb) call errquit(
     c        ' dft_alcv: deltava.ne.deltavb ',0,0)
         dovb=.true.
      else
         write(6,*) ' va vb logic wrong'
         write(6,*) ' deltava ',deltava
         write(6,*) ' deltavb ',deltavb
         call errquit(' wrong ',0,0)
      endif
      if(dova) atom=deltava
      if(dovb) atom=deltavb
c
      write(6,*) ' atom ',atom
      call dft_vatom(atom, basis, oskel,
     c     cetobfr, nocc, nbf, nmo,
     G     g_movecs,vmo)
      write(6,*) ' dova ',dova
      write(6,*) ' dovb ',dovb
      if(dova) then
c     write va on disk
         call dft_readwriteva(nmo,nbf,'w',vmo)
         return
      endif
      if(dovb) then
c     read va from disk
         call dft_readwriteva(nmo,nbf,'r',vmo_a)
         call dft_altcdeltae(e,
     n        nocc, nbf, nmo,
     v        vmo, vmo_a, eval)
      endif
      return
      end
      subroutine dft_vatom(atom, basis, oskel,
     c     cetobfr, nocc, nbf, nmo,
     G     g_movecs,vmo)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer atom
      integer basis
      integer g_movecs
      integer cetobfr(2,*)
      integer nocc
      integer nbf
      integer nmo
      logical oskel
      double precision vmo(nmo,*)
c
      integer g_atv
      integer nbf_atom
      integer k_c,l_c
      integer k_temp,l_temp
      integer k_v,l_v
c      
      if(.not.ga_duplicate(g_movecs,g_atv,'gatv'))
     &   call errquit('dft_vatom: ga_duplicate failed',1, GA_ERR)
      call ga_zero(g_atv)
      call int_1e_ga(basis, basis, g_atv, 'kinetic', oskel)
      call int_1e_ga(basis, basis, g_atv, 'potential', oskel)

      nbf_atom=cetobfr(2,atom)-cetobfr(1,atom)+1
      write(6,*) ' lo hi ',cetobfr(1,atom),cetobfr(2,atom)
      if (.not. MA_Push_Get(MT_Dbl,nbf*nbf_atom,'eigenvec',l_c,k_c))
     &   call errquit('dft_vatom: failed to alloc',999, MA_ERR)
      if (.not. MA_Push_Get(MT_Dbl,nbf_atom*nbf_atom,'vat',l_v,k_v))
     &   call errquit('dft_vatom: failed to alloc',999, MA_ERR)
      if (.not. MA_Push_Get(MT_Dbl,nbf*nbf_atom,'t',l_temp,k_temp))
     &   call errquit('dft_vatom: failed to alloc',999, MA_ERR)
      if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2e'
      call dft_vatom_compute(g_atv,g_movecs,
     c     nbf,nmo,cetobfr(1,atom),cetobfr(2,atom),nocc,
     d     dbl_mb(k_c),dbl_mb(k_v),dbl_mb(k_temp),vmo)
      if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2e'

      if (.not.ma_chop_stack(l_c))
     &   call errquit('dft_vatom: cannot chop stack',98, MA_ERR)
      if(.not. ga_destroy(g_atv)) call errquit
     $     ('dft_vatom: destroy g_atv failed', 2700, GA_ERR)

      
      return
      end
      subroutine dft_vatom_compute(g_atv,g_movecs,
     c     nbf,nmo,nbf_lo,nbf_hi,nocc,
     c     c,v,temp,v_mo)
      implicit none
#include "mafdecls.fh"      
      integer g_atv,g_movecs
      integer nmo,nbf
      integer nbf_lo,nbf_hi
      integer nocc
      double precision c(nbf_hi-nbf_lo+1,nmo)
      double precision v(nbf_hi-nbf_lo+1,*)
      double precision v_mo(nmo,nmo)
c
      double precision temp(nmo,nbf_hi-nbf_lo+1)
      integer lgt
c
      lgt=nbf_hi-nbf_lo+1
      if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2e'
      write(6,*) ' zero c hi,nmo ',nbf_hi-nbf_lo+1,nmo
      c(nbf_hi-nbf_lo+1,nmo)=0d0
      if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2e'
      write(6,*) ' zeroed c hi,nmo '
      
      call ga_get(g_movecs,nbf_lo,nbf_hi,1,nbf,c,nbf_hi-nbf_lo+1)
c      call ga_get(g_movecs,nbf_lo,nbf_hi,1,nbf,c,1)
      if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2e'
      write(6,*) ' got c vec'
      
      call ga_get(g_atv,nbf_lo,nbf_hi,nbf_lo,nbf_hi,v,nbf_hi-nbf_lo+1)
      if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2e'
      
c     v_mo=c^T * v_at * c
c      
c     temp=c^t * v_at
      call dgemm('t','n',nmo,lgt,lgt,1d0,
     A     c,lgt,v,lgt,0d0,temp,nmo)
            if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2e'

c     v_mo=temp * c
      call dgemm('n','n',nmo,nmo,lgt,1d0,
     A           temp,nmo,c,lgt,0d0,v_mo,nmo)
            if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2e'

      return
      end
      subroutine dft_altcdeltae(e,
     n    nocc, nbf, nmo,
     v     vmo_b, vmo_a, eval)
      implicit none
      double precision e
      integer nocc, nbf, nmo
      double precision eval(*)
      double precision vmo_b(nmo,*)
      double precision vmo_a(nmo,*)
c
      integer i,j,a
c      
c     compute deltav = vb - va
      do j=1,nmo
         do i=1,nmo
            if(i.eq.1.and.j.eq.1)
     w           write(6,*) ' v ',i,j,vmo_b(i,j),vmo_a(i,j)
            vmo_b(i,j)= vmo_b(i,j)-vmo_a(i,j)
            if(i.eq.1.and.j.eq.1)
     w           write(6,*) ' deltav ',i,j,vmo_b(i,j)
         enddo
      enddo
c      
c     vmo_b aka deltav
c      
c     deltaV_ii
      e=0d0
      write(6,*) ' nocc ',nocc
      write(6,*) ' nmo ',nmo
      do i=1,nocc
         e=e+vmo_b(i,i)
      enddo
      write(6,*) ' e DeltaV  ',e
c     - (deltaV_ia)^2/(e_a - e_i)
      do i=1,nocc
         do a=nocc+1,nmo
            e = e - vmo_b(i,a)**2/(eval(a)-eval(i))
         enddo
      enddo
      write(6,*) ' e DeltaV2 ',e

      return
      end
#define VUNIT 64      
      subroutine dft_readwriteva(nmo,nbf,what,vmo)
      implicit none
      integer nmo,nbf
      character*1 what
      double precision vmo(nmo,nmo)
c
      integer i,j,iunit
      double precision mysum
c
      rewind(VUNIT)
      write(6,*) ' what is ',what
      if(what.eq.'w') then
         write(VUNIT) nbf,nmo
      elseif(what.eq.'r') then
         read(VUNIT) nbf,nmo
      else
         call errquit('rw wrong arg',0,0)
      endif
      mysum=0d0
      do i=1,nmo
         if(what.eq.'w') then
            call swrite(VUNIT, vmo(1,i), nmo)
         elseif(what.eq.'r') then
            call sread(VUNIT, vmo(1,i), nmo)
         endif
         do j=1,nmo
            mysum=mysum+vmo(j,i)**2
         enddo
      enddo
      write(6,*) 'mysum ',mysum
      return
      end
