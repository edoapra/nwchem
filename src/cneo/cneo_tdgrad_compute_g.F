      subroutine cneo_tdgrad_compute_g(rtdb,ihdl_geom,ihdl_bfao,tda,
     +  nat,nao,
     +  nfc,naoc,nocc,nav,nfv,ipol,nroots,oskel,tol2e,g_mo,
     +  g_p,g_wp,g_xpy,g_xmy,g_g,kfac,lhashf,otriplet,
     +  nuc_num,nuc_cent,ihdl_bfn,nao_n,naoc_n,nav_n,g_mo_n,g_z_n)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "stdio.fh"
#include "util.fh"
#include "rtdb.fh"
#include "msgids.fh"
      integer rtdb      ! the runtime database handle
      integer ihdl_geom ! the geometry handle
      integer ihdl_bfao ! the basis set handleA
      logical tda       ! True if Tamm-Dancoff approximation
      integer g_mo(2)   ! the molecular orbitals
      integer g_p(2)    ! one-particle difference density matrix in MO basis
      integer g_wp(2)   ! the energy weighted difference density matrix in MO basis
      integer g_xpy(2)  ! global arrays to hold X+Y
      integer g_xmy(2)  ! global arrays to hold X-Y
      integer nat       ! the number of atoms in the molecule
      integer nao       ! the number of atomic basis functions
      integer nfc(2)    ! the number of frozen cores
      integer naoc(2)   ! the number of active occupied orbitals
      integer nocc(2)   ! the number of occupied orbitals
      integer nav(2)    ! the number of active virtual orbitals
      integer nfv(2)    ! the number of frozen virtuals
      integer ipol      ! =1 (restricted), =2 (unrestricted)
      integer nroots    ! the number of roots
      logical oskel     ! Use of symmetry in Fock build
c
      logical lhashf    ! =.true.  hybrid functionals
                        ! =.false. otherwise
      logical otriplet  ! =.true.  triplet excited states
                        ! =.false. singlet excited states
                        ! set to .false. for TDUDFT
c
      double precision kfac  ! the weight of the Hartree-Fock exchange
                             ! contributions
      double precision tol2e ! 2-electron integral tolerance
c     Nulcear part:
      integer nuc_num
      integer nuc_cent(nuc_num)
      integer ihdl_bfn(nuc_num)
      integer nao_n(nuc_num), naoc_n(nuc_num), nav_n(nuc_num)
      integer g_mo_n(nuc_num)
      integer g_z_n(nuc_num)
c     Output:
      integer g_g
c     Local variables:
      integer i, ip, ir, j
      integer bases(nuc_num+1)
      integer max1e
      integer max2e
      integer mscratch_1e
      integer mscratch_2e
      integer lbuf
      integer lscratch
      integer nprocs

      integer g_dao(2)
      integer g_dao_n(nuc_num)
      integer g_p_ao
      integer g_z_n_ao(nuc_num)
      integer g_z_mo(nuc_num)
      integer idim(3), ichnk(3)
      integer ldim(2)   ! leading dimension for ga_get

      integer max_at_bf, max_at_bf_n, max_at_bf_tmp
      integer max_sh_bf, max_sh_bf_n, max_sh_bf_tmp
      integer blen, blen_n
      integer lsqatom, lsqatom_n
      integer maxblen, maxsh, maxq
      parameter (maxblen=36, maxsh=10)
      parameter (maxq=maxsh**4)

      integer l_buf, k_buf, l_scr, k_scr, l_labels, k_labels
      integer l_dens, k_dens, l_ndens, k_ndens, l_ndens2, k_ndens2
      integer l_pdm12, k_pdm12
      integer lforce
      integer l_force, k_force, l_frc_nuc, k_frc_nuc
      integer l_frc_en, k_frc_en, l_frc_nn, k_frc_nn

      integer alo(3),ahi(3), blo(3), bhi(3)
      integer g_tmp_n(nuc_num), g_tmp

      integer nsh, nsh_n(nuc_num)
      integer nshblocks, nshblocks_n(nuc_num)
      integer ld_e, lh_e, ld_n, lh_n, ld_n1, ld_n2, lh_n2
      integer ld_e2,lh_e2
      integer l_shmap, l_shglo, l_shghi, l_bfmap, l_rbfmap, 
     $     l_bfglo, l_bfghi, l_list, l_q4, l_bftoat
      integer k_shmap, k_shglo, k_shghi, k_bfmap, k_rbfmap, 
     $     k_bfglo, k_bfghi, k_list, k_q4, k_bftoat
      integer l_shbflo, l_shbfhi
      integer k_shbflo, k_shbfhi

      integer l_shmap_n(nuc_num), l_shglo_n(nuc_num)
      integer l_shghi_n(nuc_num), l_bfmap_n(nuc_num)
      integer l_rbfmap_n(nuc_num), l_bfglo_n(nuc_num)
      integer l_bfghi_n(nuc_num), l_labels_n(nuc_num)
      integer l_list_n(nuc_num), l_bftoat_n(nuc_num)
      integer k_shmap_n(nuc_num), k_shglo_n(nuc_num)
      integer k_shghi_n(nuc_num), k_bfmap_n(nuc_num)
      integer k_rbfmap_n(nuc_num), k_bfglo_n(nuc_num)
      integer k_bfghi_n(nuc_num), k_labels_n(nuc_num)
      integer k_list_n(nuc_num), k_bftoat_n(nuc_num)
      integer l_shbflo_n(nuc_num), l_shbfhi_n(nuc_num)
      integer k_shbflo_n(nuc_num), k_shbfhi_n(nuc_num)

c     Functions:
      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      character*32 pname
      logical status
      integer iat,lsqa

      pname = 'cneo_tdgrad_compute_g: '

      call tddft_grad_compute_g(rtdb,ihdl_geom,ihdl_bfao,tda,nat,
     +   nao,nfc,naoc,nocc,nav,nfv,ipol,nroots,oskel,tol2e,g_mo,
     +   g_p,g_wp,g_xpy,g_xmy,g_g,kfac,lhashf,otriplet)
      bases(1) = ihdl_bfao
      do i=1,nuc_num
         bases(i+1) = ihdl_bfn(i)
      enddo
      call int_init(rtdb,nuc_num+1,bases)
      call cneo_schwarz_init(ihdl_geom,nuc_num+1,bases)
      call int_terminate()
      call intd_init(rtdb,nuc_num+1,bases)
      call int_mem(max1e, max2e, mscratch_1e, mscratch_2e)
      call intb_mem_2e4c(max2e, mscratch_2e) ! blocking algorithm
      max2e = max(max2e,1296*100)          ! 100 D quartets 
      lbuf = max(max1e, max2e)
      lscratch = max(mscratch_1e, mscratch_2e)

      if (.not. bas_numcont(ihdl_bfao,nsh))
     $     call errquit(pname//'could not get nsh',0, BASIS_ERR)
      do i=1,nuc_num
         if (.not. bas_numcont(ihdl_bfn(i),nsh_n(i)))
     $     call errquit(pname//'could not get nsh',0, BASIS_ERR)
      enddo

      g_p_ao = ga_create_atom_blocked(ihdl_geom,ihdl_bfao,'g_p_ao')
      g_tmp = ga_create_atom_blocked(ihdl_geom,ihdl_bfao,'g_tmp')
      do ip=1,ipol
         g_dao(ip)=ga_create_atom_blocked(ihdl_geom,ihdl_bfao,'g_d_ao')
      enddo
      call tddft_grad_compute_dao(ipol,nao,nocc,g_mo,g_dao)
      if (ipol.eq.1) call ga_scale(g_dao(1),2.0d0)
      do i=1,nuc_num
         g_dao_n(i) = ga_create_atom_blocked(ihdl_geom,ihdl_bfn(i),
     +                                                'g_d_ao_n')
         g_z_n_ao(i) = ga_create_atom_blocked(ihdl_geom,ihdl_bfn(i),
     +                                                'g_z_ao_n')
         call tddft_grad_compute_dao(1,nao_n(i),naoc_n(i),
     +                              g_mo_n(i),g_dao_n(i))
         g_tmp_n(i) = ga_create_atom_blocked(ihdl_geom,ihdl_bfn(i),
     +                                                'g_tmp_n')
      enddo

      do i=1,nuc_num
         idim(1)=nroots
         idim(2)=nao_n(i)
         idim(3)=nao_n(i)
         ichnk(1)=nroots
         ichnk(2)=-1
         ichnk(3)=-1
         if (.not. nga_create
     +            (mt_dbl,3,idim,pname//'g_z_mo',ichnk,g_z_mo(i)))
     +            call errquit(pname//'allocate g_z_mo',0,GA_ERR)
         do ir=1,nroots
            alo(1)=ir
            ahi(1)=ir
            alo(2)=1
            ahi(2)=naoc_n(i)
            alo(3)=1
            ahi(3)=nav_n(i)
            blo(1)=ir
            bhi(1)=ir
            blo(2)=1
            bhi(2)=naoc_n(i)
            blo(3)=naoc_n(i)+1
            bhi(3)=nao_n(i)
            call nga_copy_patch('n', g_z_n(i), alo, ahi, 
     +                                 g_z_mo(i), blo, bhi)
            blo(2)=naoc_n(i)+1
            bhi(2)=nao_n(i)
            blo(3)=1
            bhi(3)=naoc_n(i)
            call nga_copy_patch('t', g_z_n(i), alo, ahi, 
     +                                 g_z_mo(i), blo, bhi)
         enddo ! ir
         call ga_scale(g_z_mo(i),5.0d-1)
      enddo ! i


      if (.not. bas_nbf_ce_max(ihdl_bfao,max_at_bf))
     $     call errquit(pname//' could not get max_at_bf',0, BASIS_ERR)
      if (.not. bas_nbf_cn_max(ihdl_bfao,max_sh_bf))
     $     call errquit(pname//'could not get max_sh_bf',0, BASIS_ERR)

      max_at_bf_n = 0
      do i=1,nuc_num
         if (.not. bas_nbf_ce_max(ihdl_bfn(i),max_at_bf_tmp))
     $     call errquit(pname//' could not get max_at_bf',0, BASIS_ERR)
         if (max_at_bf_tmp .ge. max_at_bf_n) max_at_bf_n=max_at_bf_tmp
      enddo
      max_sh_bf_n = 0
      do i=1,nuc_num
         if (.not. bas_nbf_cn_max(ihdl_bfn(i),max_sh_bf_tmp))
     $     call errquit(pname//' could not get max_at_bf',0, BASIS_ERR)
         if (max_sh_bf_tmp .ge. max_sh_bf_n) max_sh_bf_n=max_sh_bf_tmp
      enddo
      blen = min(nao, maxblen, 6*max_sh_bf)
      blen_n = min(maxblen, 6*max_sh_bf_n)
      do i=1,nuc_num
         if (nao_n(i) .le. blen_n) blen_n = nao_n(i)
      enddo
      lsqatom = max_at_bf * max_at_bf
      lsqatom_n = max_at_bf_n * max_at_bf_n

      if (.not. ma_push_get(mt_dbl,lbuf,'deriv buffer',l_buf,k_buf))
     $     call errquit(pname//'could not allocate buffer',lbuf,
     &       MA_ERR)
*      call ma_summarize_allocated_blocks()

      if (.not. ma_push_get(mt_dbl,lscratch,'deriv scratch', 
     $     l_scr, k_scr))call errquit(pname//' scratch alloc failed',
     $     lscratch, MA_ERR)
      if (.not. ma_push_get(mt_int,lbuf/3,'labels',l_labels,k_labels))
     $     call errquit(pname//' could not allocate labels',lbuf/3,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,lsqatom,'local_density',
     $     l_dens,k_dens))
     $     call errquit(pname//'could not allocate l_dens',lsqatom,
     &       MA_ERR)

      if (.not. ma_push_get(mt_dbl,lsqatom_n,'local_density:n',
     $     l_ndens,k_ndens))
     $     call errquit(pname//'could not allocate l_ndens',lsqatom_n,
     &       MA_ERR)

      if (nuc_num .gt. 1) then 
         if (.not. ma_push_get(mt_dbl,lsqatom_n,'local_density:n2',
     $     l_ndens2,k_ndens2))
     $     call errquit(pname//'could not allocate l_ndens',lsqatom_n,
     &       MA_ERR)
      endif

      lforce = nat*3
      if (.not. ma_push_get(mt_dbl, lforce, 'frc_nuc', l_frc_nuc, 
     $     k_frc_nuc)) call errquit
     &    (pname//'could not allocate l_force',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'frc_en', l_frc_en, 
     $     k_frc_en)) call errquit
     &    (pname//'could not allocate l_force',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'frc_nn', l_frc_nn, 
     $     k_frc_nn)) call errquit
     &    (pname//'could not allocate l_force',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl,lforce,'forces',l_force,k_force))
     $     call errquit(pname//'could not allocate l_force',lforce,
     &       MA_ERR)

      if (.not. ma_push_get(mt_dbl, blen**2*blen_n**2, 'pdm2d',
     $     l_pdm12, k_pdm12)) call errquit
     $     (pname//'failed allocating pdm2d',blen**4, MA_ERR)
      

      if (.not. ma_push_get(mt_dbl, blen*blen, 'ld_e', lh_e, ld_e))
     $     call errquit(pname//'could not allocate ld_e',0, MA_ERR)

      if (.not. ma_push_get(mt_dbl, blen_n*blen_n, 'ld_n', lh_n, ld_n))
     $     call errquit(pname//'could not allocate ld_n',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen_n*blen_n, 'ld_n2', lh_n2, 
     $             ld_n2))
     $     call errquit(pname//'could not allocate ld_n',0, MA_ERR)

      if (.not. ma_push_get(mt_int, nsh, 'shmap', l_shmap, k_shmap))
     $     call errquit(pname//'could not allocate shmap',nsh,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shglo', l_shglo, k_shglo))
     $     call errquit(pname//'could not allocate blo',nsh, MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shghi', l_shghi, k_shghi))
     $     call errquit(pname//'could not allocate bhi',nsh, MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shbflo', l_shbflo, k_shbflo))
     $     call errquit(pname//'could not allocate bflo',nsh, MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shbfhi', l_shbfhi, k_shbfhi))
     $     call errquit(pname//'could not allocate bfhi',nsh, MA_ERR)
      if (.not. ma_push_get(mt_int, nao, 'bfglo', l_bfglo, k_bfglo))
     $     call errquit(pname//'could not allocate blo',nao, MA_ERR)
      if (.not. ma_push_get(mt_int, nao, 'bfghi', l_bfghi, k_bfghi))
     $     call errquit(pname//'could not allocate bhi',nao, MA_ERR)
      if (.not. ma_push_get(mt_int, nao, 'bfmap', l_bfmap, k_bfmap))
     $     call errquit(pname//'could not allocate bfmap',nao, MA_ERR)
      if (.not. ma_push_get(mt_int, nao, 'rbfmap', l_rbfmap, k_rbfmap))
     $     call errquit(pname//'could not allocate rbfmap',nao, MA_ERR)
      if (.not. ma_push_get(mt_int, nao, 'bftoat', l_bftoat, k_bftoat))
     $     call errquit(pname//'could not allocate bftoat',nao, MA_ERR)

      call cneo_schwarz_set(ihdl_bfao)
      call grad_shorder(ihdl_bfao, nsh, nao, maxsh, blen,
     $     nshblocks, int_mb(k_shglo), int_mb(k_shghi),
     &     int_mb(k_shmap), 
     $     int_mb(k_bfmap), int_mb(k_rbfmap), int_mb(k_bfglo),
     $     int_mb(k_bfghi), int_mb(k_shbflo), int_mb(k_shbfhi))

      do i = 1, nao
         if (.not. bas_bf2ce(ihdl_bfao, i, iat)) call errquit('bf2ce',i,
     &       BASIS_ERR)
         int_mb(k_bftoat+i-1) = iat
      end do

      do i=1,nuc_num
      if (.not. ma_push_get(mt_int, nsh_n(i), 'shmap', l_shmap_n(i), 
     $      k_shmap_n(i)))
     $     call errquit(pname//'could not allocate shmap',nsh_n(i),
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nsh_n(i), 'shglo', l_shglo_n(i), 
     $      k_shglo_n(i)))
     $    call errquit(pname//'could not allocate blo',nsh_n(i),MA_ERR)
      if (.not. ma_push_get(mt_int, nsh_n(i), 'shghi', l_shghi_n(i), 
     $      k_shghi_n(i)))
     $    call errquit(pname//'could not allocate bhi',nsh_n(i),MA_ERR)
      if (.not. ma_push_get(mt_int, nsh_n(i), 'shbflo', l_shbflo_n(i), 
     $      k_shbflo_n(i)))
     $     call errquit(pname//'could not allocate bflo',nsh_n(i), 
     $         MA_ERR)
      if (.not. ma_push_get(mt_int, nsh_n(i), 'shbfhi', l_shbfhi_n(i), 
     $      k_shbfhi_n(i)))
     $     call errquit(pname//'could not allocate bfhi',nsh_n(i), 
     $         MA_ERR)
      if (.not. ma_push_get(mt_int, nao_n(i), 'bfglo', l_bfglo_n(i), 
     $      k_bfglo_n(i)))
     $    call errquit(pname//'could not allocate blo',nao_n(i),MA_ERR)
      if (.not. ma_push_get(mt_int, nao_n(i), 'bfghi', l_bfghi_n(i), 
     $      k_bfghi_n(i)))
     $    call errquit(pname//'could not allocate bhi',nao_n(i),MA_ERR)
      if (.not. ma_push_get(mt_int, nao_n(i), 'bfmap', l_bfmap_n(i), 
     $      k_bfmap_n(i)))
     $     call errquit(pname//'could not allocate bfmap',nao_n(i), 
     $         MA_ERR)
      if (.not. ma_push_get(mt_int, nao_n(i), 'rbfmap', l_rbfmap_n(i), 
     $      k_rbfmap_n(i)))
     $     call errquit(pname//'could not allocate rbfmap',nao_n(i), 
     $         MA_ERR)
      if (.not. ma_push_get(mt_int, nao_n(i), 'bftoat', l_bftoat_n(i), 
     $      k_bftoat_n(i)))
     $     call errquit(pname//'could not allocate bftoat',nao_n(i), 
     $         MA_ERR)

      call cneo_schwarz_set(ihdl_bfn(i))
      call grad_shorder(ihdl_bfn(i), nsh_n(i), nao_n(i), maxsh,blen_n,
     $     nshblocks_n(i), int_mb(k_shglo_n(i)), int_mb(k_shghi_n(i)),
     &     int_mb(k_shmap_n(i)), int_mb(k_bfmap_n(i)), 
     $     int_mb(k_rbfmap_n(i)), int_mb(k_bfglo_n(i)),
     $     int_mb(k_bfghi_n(i)), int_mb(k_shbflo_n(i)), 
     $     int_mb(k_shbfhi_n(i)))

      do j = 1, nao_n(i)
         if (.not. bas_bf2ce(ihdl_bfn(i), j, iat)) call 
     &       errquit('bf2ce',i,BASIS_ERR)
         int_mb(k_bftoat_n(i)+j-1) = iat
      end do   ! j
      enddo ! nuc_num

      if (.not. ma_push_get(mt_int, 4*maxq, 'list', l_list, k_list))
     $     call errquit(pname//'could not allocate list',4*maxq, MA_ERR)
      if (.not. ma_push_get(mt_dbl, maxq, 'q4', l_q4, k_q4))
     $     call errquit(pname//'could not allocate q4',maxq, MA_ERR)

      nprocs=ga_nnodes()
      do ir=1,nroots
         call dfill(lforce, 0.0D0, dbl_mb(k_force), 1)
         call dfill(lforce, 0.0D0, dbl_mb(k_frc_nuc), 1)
         call dfill(lforce, 0.0D0, dbl_mb(k_frc_en), 1)
         call dfill(lforce, 0.0d0, dbl_mb(k_frc_nn), 1)
c
c     1. load elec grad
c
         alo(1) = ir
         ahi(1) = ir
         alo(2) = 1
         ahi(2) = 3
         alo(3) = 1
         ahi(3) = nat
         ldim(1) = 1
         ldim(2) = 3
         call nga_get(g_g,alo,ahi,dbl_mb(k_force),ldim)
c
c     2. calculate g_p_ao
c
         call ga_zero(g_p_ao)
         do ip = 1, ipol
            call tddft_grad_trans_mo2ao(1,nao,nfc(ip),naoc(ip),nocc(ip),
     +          nav(ip),nfv(ip),ir,1.0d0,1.0d0,"pq",g_mo(ip),g_p(ip),
     +          "pq",g_p_ao)
         enddo
c
c     3. calculate g_z_ao(nuc_num)
c
         do i=1,nuc_num
            call ga_zero(g_z_n_ao(i))
            call tddft_grad_trans_mo2ao(1,nao_n(i),0,naoc_n(i),
     +               naoc_n(i),nav_n(i),0,ir,1.0d0,1.0d0,'pq',g_mo_n(i),
     +               g_z_mo(i),'pq',g_z_n_ao(i))
         enddo
c
c     4. Calculate nuclear hcore grad \sum hcore (Dn+Zn)
c
         do i=1,nuc_num
            call ga_zero(g_tmp_n(i))
            call ga_add(1.0d0, g_z_n_ao(i),1.0d0,g_dao_n(i),g_tmp_n(i))
            call cneo_schwarz_set(ihdl_bfn(i))
            call grad1n(i,nuc_num,nuc_cent,dbl_mb(k_buf), lbuf, 
     $            dbl_mb(k_scr), lscratch,dbl_mb(k_ndens),
     $            dbl_mb(k_frc_nuc),g_tmp_n(i),ihdl_bfn(i),
     $            ihdl_geom,nprocs,nat,max_at_bf_n,rtdb)
         enddo
c
c     5. calculate en grad
c
c     5.1 \sum (Pe+De) grad(e|n) Dn
c
         do i=1,nuc_num
            call ga_zero(g_tmp)
            call ga_add(1.0d0,g_p_ao,1.0d0,g_dao(1),g_tmp)
            if (ipol.eq.2) call ga_add(1d0,g_tmp,1d0,g_dao(2),g_tmp)
            call ga_reorder(g_tmp, .true., int_mb(k_rbfmap),
     $        .true., int_mb(k_rbfmap))
            call ga_zero(g_tmp_n(i))
            call ga_copy(g_dao_n(i),g_tmp_n(i))
            call ga_reorder(g_tmp_n(i), .true., int_mb(k_rbfmap_n(i)),
     $        .true., int_mb(k_rbfmap_n(i)))

            call grad12(
     $      dbl_mb(ld_e), dbl_mb(ld_n), dbl_mb(k_frc_en), blen,
     $      blen_n, g_tmp, g_dao_n(i), ihdl_geom,ihdl_bfao,ihdl_bfn(i),
     $      nprocs, nat, lscratch, dbl_mb(k_scr), lbuf/12,dbl_mb(k_buf),
     $      int_mb(k_labels), maxq, int_mb(k_list), dbl_mb(k_q4),
     $      1.0d-9, nsh, nsh_n(i), 0, nao, nao_n(i),
     $      -1.0d0, nshblocks,
     $      int_mb(k_shmap), int_mb(k_shglo), int_mb(k_shghi), 
     $      int_mb(k_bfglo), int_mb(k_bfghi), 
     $      int_mb(k_bfmap), int_mb(k_rbfmap), 
     $      int_mb(k_bftoat), int_mb(k_shbflo), int_mb(k_shbfhi),
     $      nshblocks_n(i), 
     $      int_mb(k_shmap_n(i)), int_mb(k_shglo_n(i)), 
     $      int_mb(k_shghi_n(i)), 
     $      int_mb(k_bfglo_n(i)), int_mb(k_bfghi_n(i)), 
     $      int_mb(k_bfmap_n(i)), int_mb(k_rbfmap_n(i)), 
     $      int_mb(k_bftoat_n(i)), int_mb(k_shbflo_n(i)), 
     $      int_mb(k_shbfhi_n(i)))
c
c     5.2 \sum Zn grad(n|e) De
c
            call ga_zero(g_tmp)
            call ga_copy(g_dao(1),g_tmp)
            if (ipol.eq.2) call ga_add(1.0d0,g_dao(2),
     +                                 1.0d0,g_tmp,g_tmp)
            call ga_reorder(g_tmp, .true., int_mb(k_rbfmap),
     $        .true., int_mb(k_rbfmap))
            call ga_zero(g_tmp_n(i))
            call ga_copy(g_z_n_ao(i),g_tmp_n(i))
            call ga_reorder(g_tmp_n(i), .true., int_mb(k_rbfmap_n(i)),
     $        .true., int_mb(k_rbfmap_n(i)))
            call grad12(
     $      dbl_mb(ld_e), dbl_mb(ld_n), dbl_mb(k_frc_en), blen,
     $      blen_n, g_tmp, g_z_n_ao(i),ihdl_geom,ihdl_bfao,ihdl_bfn(i),
     $      nprocs, nat, lscratch, dbl_mb(k_scr), lbuf/12,dbl_mb(k_buf),
     $      int_mb(k_labels), maxq, int_mb(k_list), dbl_mb(k_q4),
     $      1.0d-9, nsh, nsh_n(i), 0, nao, nao_n(i),
     $      -1.0d0, nshblocks,
     $      int_mb(k_shmap), int_mb(k_shglo), int_mb(k_shghi), 
     $      int_mb(k_bfglo), int_mb(k_bfghi), 
     $      int_mb(k_bfmap), int_mb(k_rbfmap), 
     $      int_mb(k_bftoat), int_mb(k_shbflo), int_mb(k_shbfhi),
     $      nshblocks_n(i), 
     $      int_mb(k_shmap_n(i)), int_mb(k_shglo_n(i)), 
     $      int_mb(k_shghi_n(i)), 
     $      int_mb(k_bfglo_n(i)), int_mb(k_bfghi_n(i)), 
     $      int_mb(k_bfmap_n(i)), int_mb(k_rbfmap_n(i)), 
     $      int_mb(k_bftoat_n(i)), int_mb(k_shbflo_n(i)), 
     $      int_mb(k_shbfhi_n(i)))
         enddo ! end nuc_num
c
c     6. calculate nn grad
c        \sum Zn1 grad(n1|n2) Dn2      n1 != n2
c        \sum Dn1 grad(n1|n2) Dn2      n1 < n2 (ground state)
c        therefore,
c           \sum (Zn1+Dn1) grad(n1|n2) Dn2   n1 < n2
c           \sum Zn1 grad(n1|n2) Dn2         n1 > n2
c
         do i=1,nuc_num
            do j=1,nuc_num
               call ga_zero(g_tmp_n(i))
               call ga_copy(g_z_n_ao(i),g_tmp_n(i))
               if (j.gt.i) then
                  call ga_add(1.0d0,g_z_n_ao(i),
     +                        1.0d0,g_dao_n(i),g_tmp_n(i))
               else if (j.eq.i) then
                  goto 100
               endif
               call ga_reorder(g_tmp_n(i),.true.,int_mb(k_rbfmap_n(i)),
     $        .true., int_mb(k_rbfmap_n(i)))
               call ga_zero(g_tmp_n(j))
               call ga_copy(g_dao_n(j),g_tmp_n(j))
               call ga_reorder(g_tmp_n(j),.true.,int_mb(k_rbfmap_n(j)),
     $        .true., int_mb(k_rbfmap_n(j)))
            
               call grad12(
     $         dbl_mb(ld_n), dbl_mb(ld_n2), dbl_mb(k_frc_nn),blen_n,
     $         blen_n, g_tmp_n(i), g_tmp_n(j), ihdl_geom, 
     $         ihdl_bfn(i), ihdl_bfn(j), nprocs, nat, 
     $         lscratch, dbl_mb(k_scr), lbuf/12,dbl_mb(k_buf),
     $         int_mb(k_labels), maxq, int_mb(k_list), dbl_mb(k_q4),
     $         1.0d-9, nsh_n(i), nsh_n(j), 0, nao_n(i), nao_n(j),
     $         1.0d0, nshblocks_n(i),
     $         int_mb(k_shmap_n(i)), int_mb(k_shglo_n(i)), 
     $         int_mb(k_shghi_n(i)), 
     $         int_mb(k_bfglo_n(i)), int_mb(k_bfghi_n(i)), 
     $         int_mb(k_bfmap_n(i)), int_mb(k_rbfmap_n(i)), 
     $         int_mb(k_bftoat_n(i)), int_mb(k_shbflo_n(i)), 
     $         int_mb(k_shbfhi_n(i)),
     $         nshblocks_n(j), 
     $         int_mb(k_shmap_n(j)), int_mb(k_shglo_n(j)), 
     $         int_mb(k_shghi_n(j)), 
     $         int_mb(k_bfglo_n(j)), int_mb(k_bfghi_n(j)), 
     $         int_mb(k_bfmap_n(j)), int_mb(k_rbfmap_n(j)), 
     $         int_mb(k_bftoat_n(j)), int_mb(k_shbflo_n(j)), 
     $         int_mb(k_shbfhi_n(j)))

 100           continue

            enddo ! end j: nuc_num
         enddo ! end i: nuc_num

         call ga_sync()
         call ga_dgop(msg_grad_nuc, dbl_mb(k_frc_en), 3*nat, '+')
         call ga_dgop(msg_grad_nuc, dbl_mb(k_frc_nn), 3*nat, '+')
         call ga_sync()

         if (ga_nodeid() .eq. 0) then
         status = rtdb_parallel(.false.)
         write(luout,2200) 'loaded elec gradient',' ',' ',
     $           ((dbl_mb(k_force+i-1+3*(j-1)),i=1,3),j=1,nat)
         do i=0, 3*nat-1
            dbl_mb(k_force+i) = dbl_mb(k_force+i) +
     $           dbl_mb(k_frc_nuc+i) + dbl_mb(k_frc_en+i) +
     $           dbl_mb(k_frc_nn+i)
         end do
         write(luout,2200) 'quantum elec-nuc gradient',' ',' ',
     $           ((dbl_mb(k_frc_en+i-1+3*(j-1)),i=1,3),j=1,nat)
         write(luout,2200) 'total nuclear core gradient',' ',' ',
     $           ((dbl_mb(k_frc_nuc+i-1+3*(j-1)),i=1,3),j=1,nat)
         write(luout,2200) 'quantum nuc-nuc gradient',' ',' ',
     $           ((dbl_mb(k_frc_nn+i-1+3*(j-1)),i=1,3),j=1,nat)
         write(luout,2200) 'total ',
     $           'cneo tddft', ' gradient',
     $           ((dbl_mb(k_force+i-1+3*(j-1)),i=1,3),j=1,nat)

         if (.not. rtdb_put(rtdb, 'cneo_tddft:gradient', mt_dbl,3*nat, 
     $        dbl_mb(k_force)))call errquit
     $        ('gradients: could not store gradients',1, RTDB_ERR)

         status = rtdb_parallel(.true.)
         endif

         alo(1)=ir
         ahi(1)=ir
         alo(2)=1
         ahi(2)=3
         alo(2)=1
         ahi(2)=nat
         ldim(1)=1
         ldim(2)=3
         call nga_put(g_g,alo,ahi,dbl_mb(k_force),ldim)
         
      enddo ! end nroot loop
 2200       format(A,A,A/,1000(3(1x,F12.6),/))

      if (.not. ma_verify_allocator_stuff())
     $     call errquit('grad: ma corrupt',0, MA_ERR)

      if (.not. ma_chop_stack(l_buf)) 
     $     call errquit('grad: failed chopping MA stack',0, MA_ERR)

      call int_terminate()
      call cneo_schwarz_tidy(nuc_num+1, bases)
c
c     destroy global arrays
c
      status = .true.
      status = status .and. ga_destroy(g_p_ao)
      status = status .and. ga_destroy(g_tmp)
      do ip=1,ipol
         status = status .and. ga_destroy(g_dao(ip))
      enddo
      do i=1,nuc_num
         status = status .and. ga_destroy(g_dao_n(i))
         status = status .and. ga_destroy(g_z_n_ao(i))
         status = status .and. ga_destroy(g_tmp_n(i))
         status = status .and. ga_destroy(g_z_mo(i))
      enddo
      if (.not. status) call errquit(pname//'ga_destoy?',0,GA_ERR)


      end