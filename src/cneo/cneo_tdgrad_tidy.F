      subroutine cneo_tdgrad_finalize(rtdb,nsroots,ihdl_sroots,
     +          ntroots,ihdl_troots,nroots, ihdl_roots,
     +          ihdl_geom, ihdl_bfao, 
     +          nuc_num, ihdl_bfn)
c
      implicit none
c
#include "mafdecls.fh"
#include "geom.fh"
#include "bas.fh"
#include "errquit.fh"
#include "rtdb.fh"
c
      integer rtdb
      integer nsroots     ! the number of singlet roots
      integer ihdl_sroots ! the memory handle for the singlet roots
      integer ntroots     ! the number of triplet roots
      integer ihdl_troots ! the memory handle for the triplet roots
      integer nroots      ! the number of roots
      integer ihdl_roots  ! the memory handle for the roots
c
      integer ihdl_geom   ! the geometry handle
      integer ihdl_bfao   ! the basis set handle

      integer nuc_num 
      integer ihdl_bfn(nuc_num)  ! the nuclear basis set handles
      double precision tmp_charge
c
      character*32 pname
      integer i
      pname = "cneo_tdgrad_finalize: "
      call tddft_grad_finalize(nsroots,ihdl_sroots,
     +                           ntroots,ihdl_troots,
     +                           nroots, ihdl_roots,
     +                           ihdl_geom, ihdl_bfao)
      do i=1, nuc_num
         if(.not. bas_destroy(ihdl_bfn(i))) call errquit
     +      (pname//'failed to destroy nuc basis set',0,BASIS_ERR)
      enddo

      call cneo_int1er_tidy
c     Restore tot charge
      if (.not. rtdb_get(rtdb,'charge_real',MT_DBL,1,tmp_charge))
     $   call errquit(pname//'read tot charge',0,RTDB_ERR)
      if (.not. rtdb_put(rtdb,'charge',MT_DBL,1,tmp_charge))
     $   call errquit(pname//'restore tot charge',0,RTDB_ERR)

      end
      
      subroutine cneo_tdgrad_finalize_ints(rtdb,nuc_num,bas_e,bas_n)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
      integer rtdb
      integer nuc_num
      integer bas_e
      integer bas_n(nuc_num)

      integer bases(nuc_num+1)
      integer i
      character*32 pname
      integer calc_type
      pname = 'cneo_tdgrad_finailize_ints: '

      if (.not.rtdb_get(rtdb,'fock_xc:calc_type_save',mt_int,1,
     +  calc_type)) 
     + call errquit(pname//'failed to get calc_type_save',0,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'fock_xc:calc_type',mt_int,1,calc_type))
     + call errquit(pname//'failed to reset calc_type',0,RTDB_ERR)
      if (.not.rtdb_delete(rtdb,'fock_xc:calc_type_save')) 
     + call errquit(pname//'failed to delete calc_type_save',0,RTDB_ERR)
c
      call fock_2e_tidy(rtdb)
      call int_terminate()

      bases(1) = bas_e
      do i=1,nuc_num
         bases(i+1) = bas_n(i)
      enddo
      call cneo_schwarz_tidy(nuc_num+1,bases)

      end

      subroutine cneo_tdgrad_destory_mo_n(nuc_num,g_mo)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer nuc_num
      integer g_mo(nuc_num)

      integer i
      do i=1,nuc_num
         if(.not. ga_destroy(g_mo(i))) call errquit
     +         ('cneo_tdgrad_destroy_mo:ga_destroy?',0,GA_ERR)
      enddo

      end

c
c     Prints the resulting gradients for each of the excited states of interest.
c     
c     Copied from tddft_grad_print_results (change module name)
c
      subroutine cneo_tdgrad_print_result(nroots,iroot,nat,
     +                                    ihdl_geom,g_g)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "geom.fh"
c
c     Input:
c
      integer nroots        ! the number of roots
      integer iroot(nroots) ! the root numbers 
      integer nat           ! the number of atoms
      integer ihdl_geom     ! the geometry handle
      integer g_g           ! the global array containing the gradients
c
c     Local:
c
      double precision crd(3) ! coordinate of one atom
      integer ilo(3), ihi(3)  ! patch delimiters
      integer ld(3)           ! leading dimensions
      character*16 tag        ! atom label
      double precision q      ! atom charge
      integer j               ! counter
      integer ir              ! counter over roots
      integer iat             ! counter over atoms
      integer l_grd,k_grd     ! handle and index of gradient memory
c
c     Code:
c
      if (ga_nodeid().eq.0) then
        if (.not.ma_push_get(mt_dbl,3*nat,'grad',l_grd,k_grd)) call
     +    errquit('tddft_grad_print_results: could not allocate grad',
     +            0, MA_ERR)
        ld(1)  = 1
        ld(2)  = 3
        ilo(2) = 1
        ilo(3) = 1
        ihi(2) = 3
        ihi(3) = nat
        do ir=1,nroots
          write(luout,1000)iroot(ir),'x','y','z','x','y','z'
 1000     format(/,/,25x,'Root ',i3,/,/,25x,
     +            'CNEO-TDDFT ENERGY GRADIENTS',
     +           /,/,4X,'atom',15X,'coordinates',24X,'gradient',/,6X,
     +           2(1X,(3(10X,A1))))
          ilo(1) = ir
          ihi(1) = ir
          call nga_get(g_g,ilo,ihi,dbl_mb(k_grd),ld)
          do iat=1,nat
            if (.not. geom_cent_get(ihdl_geom,iat,tag,crd,q)) 
     +        call errquit
     +             ('tddft_grad_print_results: geometry corrupt?',0, 
     +              GEOM_ERR)
            write(luout,2000) iat, tag,(crd(j),j=1,3),
     +           (dbl_mb(k_grd+3*(iat-1)+j),j=0,2)
 2000       format(1X,I3,1X,A4,2(1X,3(1X,F10.6)))
          enddo
        enddo
        if (.not.ma_pop_stack(l_grd)) call
     +    errquit('tddft_grad_print_results: could not deallocate grad',
     +            0, MA_ERR)
      endif
      end
c $Id$
