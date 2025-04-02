      subroutine cneo_tdgrad_compute_w(rtdb,ihdl_geom,ihdl_bfao,tol2e,
     +           tda,ipol,nroots,nfc,naoc,nocc,nav,
     +           nfv,nao,g_mo,g_p,g_z,g_xpy,g_xmy,eps,omega,g_w,
     +           kfac,lhashf,otriplet,lhascd, nuc_num,
     +           ihdl_bfn, naoc_n,nav_n,nao_n,g_mo_n,g_z_n,
     +           tol2e_cneo)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer rtdb      ! runtime database handle
      integer ihdl_geom ! geometry handle
      integer ihdl_bfao ! basis set handle
      logical tda       ! .true. if Tamm-Dancoff approximation is used
      integer ipol      ! =1 (restricted), =2 (unrestricted)
      integer nfc(2)    ! the number of frozen core orbitals
      integer naoc(2)   ! the number of active occupied orbitals
      integer nocc(2)   ! the total number of occupied orbitals
      integer nav(2)    ! the number of active virtual orbitals
      integer nfv(2)    ! the number of frozen virtual orbitals
      integer nao       ! the number of basis functions
      integer nroots    ! the number of states to consider
      integer g_mo(2)   ! global array handle for MO coefficients
      integer g_p(2)    ! global array handle for Ppq
      integer g_z(2)    ! global array handle for Zia
      integer g_xpy(2)  ! global array handle for (X+Y)
      integer g_xmy(2)  ! global array handle for (X-Y)
c
      logical lhashf    ! =.true.  hybrid functionals
                        ! =.false. otherwise
      logical otriplet  ! =.true.  triplet excited states
                        ! =.false. singlet excited states
                        ! value does not matter for TDUDFT
      logical lhascd
c
      double precision kfac          ! the weight of the Hartree-Fock
                                     ! exchange contributions
      double precision eps(nao,2)    ! orbital eigenvalues
      double precision omega(nroots) ! the excitation energies
      double precision tol2e         ! integral tolerance
      double precision tol2e_cneo

c     === Nuclear variables ===
      integer nuc_num
      integer ihdl_bfn(nuc_num)
      integer naoc_n(nuc_num)
      integer nav_n(nuc_num)
      integer nao_n(nuc_num)
      integer g_mo_n(nuc_num)
      integer g_z_n(nuc_num)
c
c     Output:
c
      integer g_w(2)    ! global array handle for W
c
c     Local variables
c
      integer i,j
      integer idim(3)
      integer ichunk(3)
      integer g_tmp
      integer alo(3)
      character*32 pname
      logical oroot

      pname = 'cneo_tdgrad_compute_w: '

      call tddft_grad_compute_w(rtdb,ihdl_geom,ihdl_bfao,tol2e,tda,
     +     ipol,nroots,nfc,naoc,nocc,nav,nfv,nao,g_mo,g_p,g_z,
     +     g_xpy,g_xmy,eps,omega,g_w,
     +     kfac,lhashf,otriplet,lhascd)

      idim(1) = nroots
      idim(2) = naoc(1)
      idim(3) = naoc(1)
      ichunk(1) = nroots
      ichunk(2) = -1
      ichunk(3) = -1
      if(.not. nga_create(mt_dbl,3,idim,'w:g_tmp',ichunk,g_tmp))
     $      call errquit(pname//'create g_tmp',0,GA_ERR)

      do i=1,nuc_num
         call ga_zero(g_tmp)
         call cneo_tdgrad_response_ov(i,0,'enints',
     +            ihdl_bfao,ihdl_bfn(i),-2.0d0,g_z_n(i),g_mo(1),
     +            g_mo_n(i),nao,nao_n(i),naoc(1),naoc_n(i),
     +            nav(1),nav_n(i),nfc(1),0,
     +            nfv(1),0,nroots,g_tmp,ihdl_geom,
     +            nocc(1),naoc_n(i),'ij',tol2e_cneo)

         alo(1)=1
         alo(2)=1
         alo(3)=1
         call nga_add_patch(1.0d0,g_tmp,alo,idim,1.0d0,g_w(1),alo,idim,
     $            g_w(1),alo,idim)
      enddo

      if (util_print('cneo tdgrad w',print_debug)) then
            oroot = ga_nodeid().eq.0
            if (oroot) write(LuOut,*)'DEBUG: '//pname//'W'
            call tddft_grad_print_array(ipol,nroots,g_w,dble(ipol))
      endif

      if (.not. ga_destroy(g_tmp)) call errquit(pname//'destroy?',0,0)

      end

