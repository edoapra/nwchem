      subroutine cneo_dft_energy_gradient(rtdb)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "nwc_const.fh"
#include "ccneo.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"
      integer rtdb
      logical cneo_dft_energy
      external cneo_dft_energy
      logical int_normalize
      external int_normalize
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      integer geom, geom_cneo,basis
      integer bases(max_quantum_nuc+1)
      integer bases_nuc(max_quantum_nuc)
      integer i,j

      integer l_buf,k_buf,l_scr,k_scr,l_dens,k_dens
      integer l_frc_nuc,k_frc_nuc,nproc,nat
      integer lbuf,lscratch,lsqatom,lforce
      integer max1e, max2e, mscratch_1e, mscratch_2e
      character*255 movecs_in
      integer max_at_bf, max_at_bf_n, max_at_bf_tmp
      integer max_sh_bf, max_sh_bf_n, max_sh_bf_tmp
      logical status

      ! integer bas1, bas2, nbf1, nbf2
      integer lsqatom_n, l_ndens, k_ndens, l_ndens2, k_ndens2
      integer g_ndens(max_quantum_nuc), g_dens
      integer nbf, nbf_n(max_quantum_nuc)
      integer nsh, nsh_n(max_quantum_nuc)
      character*5 str_i
      integer nopen, nclosed
      integer iat, lsqa
      character*8 scftype

      integer l_frc_en, k_frc_en, l_pdm12, k_pdm12
      integer l_frc_nn, k_frc_nn
      integer nshblocks, nshblocks_n(max_quantum_nuc)
      integer blen, blen_n
      integer maxblen, maxsh, maxq
      parameter (maxblen=36, maxsh=10)
      parameter (maxq=maxsh**4)

      integer ld_e, lh_e, ld_n, lh_n, ld_n1, ld_n2, lh_n2
      integer ld_e2,lh_e2
      integer l_shmap, l_shglo, l_shghi, l_bfmap, l_rbfmap, 
     $     l_bfglo, l_bfghi, l_labels, l_list, l_q4, l_bftoat
      integer k_shmap, k_shglo, k_shghi, k_bfmap, k_rbfmap, 
     $     k_bfglo, k_bfghi, k_labels, k_list, k_q4, k_bftoat
      integer l_shbflo, l_shbfhi
      integer k_shbflo, k_shbfhi

      integer l_shmap_n(max_quantum_nuc), l_shglo_n(max_quantum_nuc)
      integer l_shghi_n(max_quantum_nuc), l_bfmap_n(max_quantum_nuc)
      integer l_rbfmap_n(max_quantum_nuc), l_bfglo_n(max_quantum_nuc)
      integer l_bfghi_n(max_quantum_nuc), l_labels_n(max_quantum_nuc)
      integer l_list_n(max_quantum_nuc), l_bftoat_n(max_quantum_nuc)
      integer k_shmap_n(max_quantum_nuc), k_shglo_n(max_quantum_nuc)
      integer k_shghi_n(max_quantum_nuc), k_bfmap_n(max_quantum_nuc)
      integer k_rbfmap_n(max_quantum_nuc), k_bfglo_n(max_quantum_nuc)
      integer k_bfghi_n(max_quantum_nuc), k_labels_n(max_quantum_nuc)
      integer k_list_n(max_quantum_nuc), k_bftoat_n(max_quantum_nuc)
      integer l_shbflo_n(max_quantum_nuc), l_shbfhi_n(max_quantum_nuc)
      integer k_shbflo_n(max_quantum_nuc), k_shbfhi_n(max_quantum_nuc)

      integer l_force, k_force
      logical grad_cneo_elec
      external grad_cneo_elec
      integer ii, jj
      character*16 tag
      double precision crd(3), q
      double precision tmp_charge


      status = cneo_dft_energy(rtdb)
      status = rtdb_cput(rtdb,'task:theory', 1, 'dft')

      status = geom_create(geom,'geometry')
      status = geom_rtdb_load(rtdb,geom,'geometry')
      
      call cneo_nuc_init(rtdb,geom)
      call cneo_geom_init(rtdb,geom_cneo)

      status = bas_create(basis,'ao basis')
      status = bas_rtdb_load(rtdb,geom_cneo,basis,'ao basis')
      status = int_normalize(rtdb,basis)
      status = grad_cneo_elec(rtdb,geom_cneo,basis)
      bases(1) = basis
      call cneo_bas_init(rtdb,geom_cneo,bases_nuc)
      do i=1,nuc_num
         status = int_normalize(rtdb,bases_nuc(i))
         bases(i+1) = bases_nuc(i)
      enddo

      if (.not. bas_numbf(basis,nbf))
     $     call errquit('grad: could not get nbf',0, BASIS_ERR)
      if (.not. bas_numcont(basis,nsh))
     $     call errquit('grad: could not get nsh',0, BASIS_ERR)
      do i=1,nuc_num
         if (.not. bas_numbf(bases_nuc(i),nbf_n(i)))
     $     call errquit('grad: could not get nbf',0, BASIS_ERR)
         if (.not. bas_numcont(bases_nuc(i),nsh_n(i)))
     $     call errquit('grad: could not get nsh',0, BASIS_ERR)
      enddo
      if (.not. geom_ncent(geom,nat))
     $     call errquit('grad: could not get natoms',0, GEOM_ERR)


      call int_init(rtdb,nuc_num+1,bases(1:nuc_num+1))
      call cneo_schwarz_init(geom_cneo,nuc_num+1,bases)
      call int_terminate()
      call intd_init(rtdb, nuc_num+1, bases)
      if (.not. geom_ncent(geom_cneo,nat))
     $     call errquit('grad: could not get natoms',0, GEOM_ERR)

      call int_mem(max1e, max2e, mscratch_1e, mscratch_2e)
      call intb_mem_2e4c(max2e, mscratch_2e) ! blocking algorithm
      max2e = max(max2e,1296*100)          ! 100 D quartets 
      lbuf = max(max1e, max2e)
      lscratch = max(mscratch_1e, mscratch_2e)
      if (.not. bas_nbf_ce_max(basis,max_at_bf))
     $     call errquit('grad: could not get max_at_bf',0, BASIS_ERR)

      max_at_bf_n = 0
      do i=1,nuc_num
         if (.not. bas_nbf_ce_max(bases_nuc(i),max_at_bf_tmp))
     $     call errquit('grad: could not get max_at_bf',0, BASIS_ERR)
         if (max_at_bf_tmp .ge. max_at_bf_n) max_at_bf_n=max_at_bf_tmp
      enddo

      max_sh_bf_n = 0
      do i=1,nuc_num
         if (.not. bas_nbf_cn_max(bases_nuc(i),max_sh_bf_tmp))
     $     call errquit('grad: could not get max_at_bf',0, BASIS_ERR)
         if (max_sh_bf_tmp .ge. max_sh_bf_n) max_sh_bf_n=max_sh_bf_tmp
      enddo

      blen = min(nbf, maxblen, 6*max_sh_bf)
      blen_n = min(maxblen, 6*max_sh_bf_n)
      do i=1,nuc_num
         if (nbf_n(i) .le. blen_n) blen_n = nbf_n(i)
      enddo

      lsqatom = max_at_bf * max_at_bf
      lsqatom_n = max_at_bf_n * max_at_bf_n

      if (.not. ma_push_get(mt_dbl,lbuf,'deriv buffer',l_buf,k_buf))
     $     call errquit('grad:could not allocate buffer',lbuf,
     &       MA_ERR)
*      call ma_summarize_allocated_blocks()

      if (.not. ma_push_get(mt_dbl,lscratch,'deriv scratch', 
     $     l_scr, k_scr))call errquit('grad: scratch alloc failed',
     $     lscratch, MA_ERR)
      if (.not. ma_push_get(mt_int,lbuf/3,'labels',l_labels,k_labels))
     $     call errquit('grad: could not allocate labels',lbuf/3,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,lsqatom,'local_density',
     $     l_dens,k_dens))
     $     call errquit('grad:could not allocate l_dens',lsqatom,
     &       MA_ERR)

      if (.not. ma_push_get(mt_dbl,lsqatom_n,'local_density:n',
     $     l_ndens,k_ndens))
     $     call errquit('grad:could not allocate l_ndens',lsqatom_n,
     &       MA_ERR)

      if (nuc_num .gt. 1) then 
         if (.not. ma_push_get(mt_dbl,lsqatom_n,'local_density:n2',
     $     l_ndens2,k_ndens2))
     $     call errquit('grad:could not allocate l_ndens',lsqatom_n,
     &       MA_ERR)
      endif

      lforce = nat*3
      if (.not. ma_push_get(mt_dbl, lforce, 'frc_nuc', l_frc_nuc, 
     $     k_frc_nuc)) call errquit
     &    ('grad:could not allocate l_force',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'frc_en', l_frc_en, 
     $     k_frc_en)) call errquit
     &    ('grad:could not allocate l_force',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lforce, 'frc_nn', l_frc_nn, 
     $     k_frc_nn)) call errquit
     &    ('grad:could not allocate l_force',lforce, MA_ERR)
      if (.not. ma_push_get(mt_dbl,lforce,'forces',l_force,k_force))
     $     call errquit('grad:could not allocate l_force',lforce,
     &       MA_ERR)

      if (.not. rtdb_get(rtdb, 'dft:gradient',mt_dbl,3*nat,
     $      dbl_mb(k_force))) call errquit
     $      ('cneo_grad: cannot load dft gradients',0,0)
      write(luout,2200) 'loaded ',
     $           'dft elec', ' gradient',
     $           ((dbl_mb(k_force+i-1+3*(j-1)),i=1,3),j=1,nat)

      do i=1,nuc_num
         movecs_in = ' '
         write(str_i, '(I0)') i
         g_ndens(i)=ga_create_atom_blocked(geom_cneo,bases_nuc(i),
     $            'grad1n:tmp ndens'//trim(str_i))
         if (.not. rtdb_cget(rtdb,'cneo:output vectors:n'//trim(str_i),
     $         1,movecs_in))
     $        call errquit('gradients: DFT MO vectors not defined',0,
     &       RTDB_ERR)
         call cneo_dens_from_file(movecs_in,geom_cneo,bases_nuc(i),
     $                                   g_ndens(i),nbf_n(i),1,0,'RHF')
         call ga_dscal(g_ndens(i), 5.0d-1)
      enddo
      movecs_in = ' '
      if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
     $        call errquit('gradients: DFT MO vectors not defined',0,
     &       RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'scf:nopen', mt_int, 1, nopen))
     $        nopen = 0
      if (.not. rtdb_get(rtdb, 'scf:nclosed', mt_int, 1, nclosed))
     $        call errquit('gradients: get of nclosed? ', 0, RTDB_ERR)
      if (.not. rtdb_cget(rtdb, 'dft:scftype', 1, scftype))
     &        call errquit('gradients: scftype not defined',0,
     &       RTDB_ERR)

      g_dens = ga_create_atom_blocked(geom_cneo, basis, 
     $        'density matrix')
      call cneo_dens_from_file(movecs_in,geom_cneo,basis,
     $                    g_dens,nbf,nclosed,nopen,scftype)

      if (.not. ma_push_get(mt_dbl, blen**2*blen_n**2, 'pdm2d',
     $     l_pdm12, k_pdm12)) call errquit
     $     ('gradient:failed allocating pdm2d',blen**4, MA_ERR)

      nproc = ga_nnodes()
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_nuc), 1)
      call dfill(lforce, 0.0D0, dbl_mb(k_frc_en), 1)
      call dfill(lforce, 0.0d0, dbl_mb(k_frc_nn), 1)

      do i=1,nuc_num
         call cneo_schwarz_set(bases_nuc(i))
         call grad1n(i,nuc_num,nuc_cent,dbl_mb(k_buf), lbuf, 
     $            dbl_mb(k_scr), lscratch,dbl_mb(k_ndens),
     $            dbl_mb(k_frc_nuc),g_ndens(i),bases_nuc(i),
     $            geom_cneo,nproc,nat,max_at_bf_n,rtdb)
      enddo


      if (.not. ma_push_get(mt_dbl, blen*blen, 'ld_e', lh_e, ld_e))
     $     call errquit('grad:could not allocate ld_e',lsqa, MA_ERR)

      if (.not. ma_push_get(mt_dbl, blen_n*blen_n, 'ld_n', lh_n, ld_n))
     $     call errquit('grad:could not allocate ld_n',lsqa, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen_n*blen_n, 'ld_n2', lh_n2, 
     $             ld_n2))
     $     call errquit('grad:could not allocate ld_n',lsqa, MA_ERR)

      if (.not. ma_push_get(mt_int, nsh, 'shmap', l_shmap, k_shmap))
     $     call errquit('grad:could not allocate shmap',nsh,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shglo', l_shglo, k_shglo))
     $     call errquit('grad:could not allocate blo',nsh, MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shghi', l_shghi, k_shghi))
     $     call errquit('grad:could not allocate bhi',nsh, MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shbflo', l_shbflo, k_shbflo))
     $     call errquit('grad:could not allocate bflo',nsh, MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shbfhi', l_shbfhi, k_shbfhi))
     $     call errquit('grad:could not allocate bfhi',nsh, MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfglo', l_bfglo, k_bfglo))
     $     call errquit('grad:could not allocate blo',nbf, MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfghi', l_bfghi, k_bfghi))
     $     call errquit('grad:could not allocate bhi',nbf, MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfmap', l_bfmap, k_bfmap))
     $     call errquit('grad:could not allocate bfmap',nbf, MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'rbfmap', l_rbfmap, k_rbfmap))
     $     call errquit('grad:could not allocate rbfmap',nbf, MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bftoat', l_bftoat, k_bftoat))
     $     call errquit('grad:could not allocate bftoat',nbf, MA_ERR)

      call cneo_schwarz_set(basis)
      call grad_shorder(basis, nsh, nbf, maxsh, blen,
     $     nshblocks, int_mb(k_shglo), int_mb(k_shghi),
     &     int_mb(k_shmap), 
     $     int_mb(k_bfmap), int_mb(k_rbfmap), int_mb(k_bfglo),
     $     int_mb(k_bfghi), int_mb(k_shbflo), int_mb(k_shbfhi))

      do i = 1, nbf
         if (.not. bas_bf2ce(basis, i, iat)) call errquit('bf2ce',i,
     &       BASIS_ERR)
         int_mb(k_bftoat+i-1) = iat
      end do
      call ga_reorder(g_dens, .true., int_mb(k_rbfmap),
     $        .true., int_mb(k_rbfmap))

      do i=1,nuc_num
      if (.not. ma_push_get(mt_int, nsh_n(i), 'shmap', l_shmap_n(i), 
     $      k_shmap_n(i)))
     $     call errquit('grad:could not allocate shmap',nsh_n(i),
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nsh_n(i), 'shglo', l_shglo_n(i), 
     $      k_shglo_n(i)))
     $     call errquit('grad:could not allocate blo',nsh_n(i), MA_ERR)
      if (.not. ma_push_get(mt_int, nsh_n(i), 'shghi', l_shghi_n(i), 
     $      k_shghi_n(i)))
     $     call errquit('grad:could not allocate bhi',nsh_n(i), MA_ERR)
      if (.not. ma_push_get(mt_int, nsh_n(i), 'shbflo', l_shbflo_n(i), 
     $      k_shbflo_n(i)))
     $     call errquit('grad:could not allocate bflo',nsh_n(i), 
     $         MA_ERR)
      if (.not. ma_push_get(mt_int, nsh_n(i), 'shbfhi', l_shbfhi_n(i), 
     $      k_shbfhi_n(i)))
     $     call errquit('grad:could not allocate bfhi',nsh_n(i), 
     $         MA_ERR)
      if (.not. ma_push_get(mt_int, nbf_n(i), 'bfglo', l_bfglo_n(i), 
     $      k_bfglo_n(i)))
     $     call errquit('grad:could not allocate blo',nbf_n(i), MA_ERR)
      if (.not. ma_push_get(mt_int, nbf_n(i), 'bfghi', l_bfghi_n(i), 
     $      k_bfghi_n(i)))
     $     call errquit('grad:could not allocate bhi',nbf_n(i), MA_ERR)
      if (.not. ma_push_get(mt_int, nbf_n(i), 'bfmap', l_bfmap_n(i), 
     $      k_bfmap_n(i)))
     $     call errquit('grad:could not allocate bfmap',nbf_n(i), 
     $         MA_ERR)
      if (.not. ma_push_get(mt_int, nbf_n(i), 'rbfmap', l_rbfmap_n(i), 
     $      k_rbfmap_n(i)))
     $     call errquit('grad:could not allocate rbfmap',nbf_n(i), 
     $         MA_ERR)
      if (.not. ma_push_get(mt_int, nbf_n(i), 'bftoat', l_bftoat_n(i), 
     $      k_bftoat_n(i)))
     $     call errquit('grad:could not allocate bftoat',nbf_n(i), 
     $         MA_ERR)

      call cneo_schwarz_set(bases_nuc(i))
      call grad_shorder(bases_nuc(i), nsh_n(i), nbf_n(i), maxsh,blen_n,
     $     nshblocks_n(i), int_mb(k_shglo_n(i)), int_mb(k_shghi_n(i)),
     &     int_mb(k_shmap_n(i)), int_mb(k_bfmap_n(i)), 
     $     int_mb(k_rbfmap_n(i)), int_mb(k_bfglo_n(i)),
     $     int_mb(k_bfghi_n(i)), int_mb(k_shbflo_n(i)), 
     $     int_mb(k_shbfhi_n(i)))

      do j = 1, nbf_n(i)
         if (.not. bas_bf2ce(bases_nuc(i), j, iat)) call 
     &       errquit('bf2ce',i,BASIS_ERR)
         int_mb(k_bftoat_n(i)+j-1) = iat
      end do   ! j
      call ga_reorder(g_ndens(i), .true., int_mb(k_rbfmap_n(i)),
     $        .true., int_mb(k_rbfmap_n(i)))
      enddo ! nuc_num

      if (.not. ma_push_get(mt_int, 4*maxq, 'list', l_list, k_list))
     $     call errquit('grad:could not allocate list',4*maxq, MA_ERR)
      if (.not. ma_push_get(mt_dbl, maxq, 'q4', l_q4, k_q4))
     $     call errquit('grad:could not allocate q4',maxq, MA_ERR)

      do i=1,nuc_num
         call grad12(
     $      dbl_mb(ld_e), dbl_mb(ld_n), dbl_mb(k_frc_en), blen,
     $      blen_n, g_dens, g_ndens(i), geom_cneo, basis, bases_nuc(i),
     $      nproc, nat, lscratch, dbl_mb(k_scr), lbuf/12,dbl_mb(k_buf),
     $      int_mb(k_labels), maxq, int_mb(k_list), dbl_mb(k_q4),
     $      1.0d-9, nsh, nsh_n(i), nopen, nbf, nbf_n(i),
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

         do j=i+1,nuc_num
            call grad12(
     $         dbl_mb(ld_n), dbl_mb(ld_n2), dbl_mb(k_frc_nn),blen_n,
     $         blen_n, g_ndens(i), g_ndens(j), geom_cneo, 
     $         bases_nuc(i), bases_nuc(j), nproc, nat, 
     $         lscratch, dbl_mb(k_scr), lbuf/12,dbl_mb(k_buf),
     $         int_mb(k_labels), maxq, int_mb(k_list), dbl_mb(k_q4),
     $         1.0d-9, nsh_n(i), nsh_n(j), 0, nbf_n(i), nbf_n(j),
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
         enddo ! j loop for nn energy grad
      enddo ! i loop for en energy grad


      if (ga_nodeid() .eq. 0) then
         status = rtdb_parallel(.false.)
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
     $           'cneo dft', ' gradient',
     $           ((dbl_mb(k_force+i-1+3*(j-1)),i=1,3),j=1,nat)

         if (.not. rtdb_put(rtdb, 'cneo:gradient', mt_dbl, 3*nat, 
     $        dbl_mb(k_force)))call errquit
     $        ('gradients: could not store gradients',1, RTDB_ERR)

         status = rtdb_parallel(.true.)
      endif

      call ga_sync
      if (ga_nodeid() .eq. 0 
     $     .and. util_print('forces', print_low)) then
         write(luout,1000) 'CNEO',
     $        'x','y','z','x','y','z'
         do 30, i=1, nat
            if (.not. geom_cent_get(geom, i, tag, crd, q)) call errquit
     $           ('gradients: geometry corrupt?',0, GEOM_ERR)
            write(luout,2000) i, tag,(crd(j),j=1,3),
     $           (dbl_mb(k_force+3*(i-1)+j),j=0,2)
 30      continue
         write(luout,*)
 1000    format(/,/,25X,A,' ENERGY GRADIENTS',/,/,4X,'atom',15X,
     $        'coordinates',
     $        24X,'gradient',/,6X,2(1X,(3(10X,A1))))
 2000    format(1X,I3,1X,A4,2(1X,3(1X,F10.6)))
         call util_flush(luout)
      end if

 2200       format(A,A,A/,1000(3(1x,F12.6),/))

      if (.not. ma_verify_allocator_stuff())
     $     call errquit('grad: ma corrupt',0, MA_ERR)

      if (.not. ma_chop_stack(l_buf)) 
     $     call errquit('grad: failed chopping MA stack',0, MA_ERR)
c     
      do i=1, nuc_num
         if (.not. ga_destroy(g_ndens(i))) call errquit 
     $        ('error destroying density', 1, GA_ERR)
      end do
      if (.not. ga_destroy(g_dens)) call errquit 
     $        ('error destroying density', 1, GA_ERR)

      call cneo_schwarz_tidy(nuc_num+1,bases)
      call intd_terminate()
      if (.not. geom_destroy(geom)) call errquit
     $      ('cneo_grad:geom_destroy?',0,0)
      if (.not. geom_destroy(geom_cneo)) call errquit
     $      ('cneo_grad:geom_destroy?',0,0)
      do i=1,nuc_num+1
         if (.not. bas_destroy(bases(i))) call errquit
     $      ('cneo_grad:bas_destroy?',0,0)
      enddo
      if (.not. rtdb_get(rtdb,'charge_real',MT_DBL,1,tmp_charge))
     $   call errquit('cneo_dftgrad:read tot charge',0,0)
      if (.not. rtdb_put(rtdb,'charge',MT_DBL,1,tmp_charge))
     $   call errquit('cneo_dftgrad:restore tot charge',0,0)


      end


c      Modified from grad1.F

      subroutine grad1n(idx, nuc_num, nuc_cent,
     $            H, lbuf, scr, lscr, ndens,
     $            frc_nuc,
     $            g_ndens, basis, geom, nproc,
     $            nat, max_at_bf, rtdb)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"

      integer idx, nuc_num, nuc_cent(nuc_num)
      integer lbuf, lscr, g_ndens, g_force
      integer basis, geom, nat, nproc, rtdb
      integer max_at_bf
      double precision H(lbuf)
      double precision frc_kin(3,nat), frc_nuc(3,nat)
      double precision ndens(max_at_bf, max_at_bf)
      double precision scr(lscr)

      integer qmnuc, iat, i
      logical foundit
      integer ijatom, next, iat1, iat2, iat3, ish1, ish2,
     $     iab1f, iab1l, iab2f, iab2l, iac1f, iac1l, iac2f, iac2l,
     $     if1, il1, if2, il2,
     $     icart, ic, nint, ip1, ip2

      double precision crd1, crd2 ! atomic coordinates
      dimension crd1(3), crd2(3)

      integer idatom
      dimension idatom(2)

      double precision dE, dx, dy, dz, qfac, fact, q1, q2

      logical status, pointforce,dobq

      character*16 name

      integer bq_ncent
      integer i_qbq,i_cbq
      double precision r12
c
      integer nxtask, task_size
      external nxtask
      integer ecp_channels

      ecp_channels = 1
      qmnuc = nuc_cent(idx)

      iat1 = qmnuc
      status = bas_ce2cnr(basis,iat1,iac1f,iac1l)
      status = bas_ce2bfr(basis,iat1,iab1f,iab1l)
      call ga_get(g_ndens,iab1f,iab1l,iab1f,iab1l,
     $         ndens,max_at_bf)
      do ish1 = iac1f,iac1l
         do ish2 = iac1f,ish1
            status = bas_cn2bfr(basis,ish1,if1,il1)
            if1 = if1 - iab1f + 1
            il1 = il1 - iab1f + 1
            status = bas_cn2bfr(basis,ish2,if2,il2)
            if2 = if2 - iab1f + 1
            il2 = il2 - iab1f + 1
            call intd_1eh1(basis,ish1,basis,ish2,lscr,scr,
     $               lbuf,H)
            ! TODO:  scale H with charge
            ic=1
            do iat2=1,nat
               do icart=1,3
                  dE=0d0
                  do ip1=if1,il1
                     do ip2=if2,il2
                        ! minus sigh for positive nuc charge?
                        dE=dE-ndens(ip1,ip2) * H(ic)
                        ic = ic+1
                     enddo !ip2
                  enddo! ip1
                  if (ish1.ne.ish2) dE = dE+dE
                  frc_nuc(icart,iat2) = frc_nuc(icart,iat2)+dE
               enddo !icart
            enddo !iat2

         enddo !ish2
      enddo !ish1

      end

c      Modified from grad2.F
      subroutine grad12(d1_ij, d2_kl, frc, blen1, blen2,
     $         g_dens1, g_dens2, geom, bas1, bas2,
     $         nproc, nat, lscr, scr, leri, eri, labels,
     $         maxq, list, q4, tol2e, nsh1, nsh2,
     $         nopen, nbf1, nbf2, jfac,
     $         nshblocks1, shmap1, shglo1, shghi1, bfglo1, 
     $         bfghi1, bfmap1, rbfmap1, bftoat1, shbflo1, shbfhi1,
     $         nshblocks2, shmap2, shglo2, shghi2, bfglo2, 
     $         bfghi2, bfmap2, rbfmap2, bftoat2, shbflo2, shbfhi2)

      implicit none 
#include "errquit.fh"
#include "sym.fh"
#include "global.fh"
#include "schwarz.fh"
#include "util.fh"
#include "bas.fh"

      integer blen1, blen2
      double precision d1_ij(blen1,blen1),d2_kl(blen2,blen2)
      double precision frc(3,nat)
      integer g_dens1, g_dens2, g_force
      integer geom, bas1, bas2
      integer nproc, nat, lscr, leri
      double precision scr(lscr), eri(3,4,leri)
      integer labels(leri,4)
      integer maxq
      integer list(maxq,4)
      double precision q4(maxq)
      double precision tol2e
      integer nsh1, nsh2
      logical oactive(2)
      integer nopen, nbf1, nbf2
      double precision jfac

      ! copied from grad2.F
      double precision pdm2d(blen2,blen2,blen1,blen1) ! [scratch] AO 2-pdm
      integer nshblocks1         ! [input] No. of shell blocks
      integer shmap1(nsh1)        ! [input] Map from new to old shell order
      integer shglo1(nshblocks1)  ! [input] First new shell in block
      integer shghi1(nshblocks1)  ! [input] Last new shell in block
      integer bfglo1(nshblocks1)  ! [input] First new basis func in block
      integer bfghi1(nshblocks1)  ! [input] Last new basis func in block
      integer bfmap1(nbf1)        ! [input] Map from new to old basis func order
      integer rbfmap1(nbf1)       ! [input] Map from old to new basis func order
      integer bftoat1(nbf1)       ! [input] Map from old basis func to atom
      integer shbflo1(nsh1)       ! [input] First new bf in new shell
      integer shbfhi1(nsh1)       ! [input] Last new bf in new shell
      integer nshblocks2         ! [input] No. of shell blocks
      integer shmap2(nsh2)        ! [input] Map from new to old shell order
      integer shglo2(nshblocks2)  ! [input] First new shell in block
      integer shghi2(nshblocks2)  ! [input] Last new shell in block
      integer bfglo2(nshblocks2)  ! [input] First new basis func in block
      integer bfghi2(nshblocks2)  ! [input] Last new basis func in block
      integer bfmap2(nbf2)        ! [input] Map from new to old basis func order
      integer rbfmap2(nbf2)       ! [input] Map from old to new basis func order
      integer bftoat2(nbf2)       ! [input] Map from old basis func to atom
      integer shbflo2(nsh2)       ! [input] First new bf in new shell
      integer shbfhi2(nsh2)       ! [input] Last new bf in new shell
C-------------------------local variables-----------------------------
      integer  next, nint, ijklblock
      double precision scale, q4max

      integer nxtask, task_size
      external nxtask

      double precision psum
      integer i, j, k, l
      integer iilo, jjlo, kklo, lllo
      integer iihi, jjhi, kkhi, llhi
      integer ish, jsh, ksh, lsh, idim, jdim, kdim, ldim
      integer ibflo, ibfhi, jbflo, jbfhi, kbflo, kbfhi, lbflo, lbfhi
      integer ishlo, ishhi, jshlo, jshhi, kshlo, kshhi, lshlo, lshhi
      integer lblockhi, jshtop, kshtop, lshtop
      integer ishblock, jshblock, kshblock, lshblock
      integer iish, jjsh, kksh, llsh, iiat, iat
      integer nq, integ, atoms(4)
*      integer numq
      double precision smax, sij, sijkl, block_eff, p
*     double precision integ_acc
      
      integer ibflo_prev, jbflo_prev, kbflo_prev, lbflo_prev ! for caching
      
      logical oij, okl, oikjl, ouhf, orohf, omcscf, omore, exso
      logical intbd_init4c, test_intbd_2e4c, intbd_2e4c
      external intbd_init4c, test_intbd_2e4c, intbd_2e4c

      logical oskel
      double precision dtmp(blen1,blen1)
      double precision pdmtmp(blen1)
      double precision coeff
      oskel = .true.
      oactive(1) = .true.
      oactive(2) = .true.

      smax = schwarz_max()
      q4max = 8.0d0*(sym_number_ops(geom) + 1)*10000.0d0
      
      call int_acc_std()
      ibflo_prev = -1
      jbflo_prev = -1
      kbflo_prev = -1
      lbflo_prev = -1

      task_size = nshblocks1*nshblocks2
      task_size = task_size*(task_size+1)/2
      task_size = max(1,task_size/(20*nproc))

      ijklblock = 0
      next = nxtask(nproc,task_size)
      do ishblock = nshblocks1, 1, -1
         ishlo = shglo1(ishblock)
         ishhi = shghi1(ishblock)
         ibflo = bfglo1(ishblock)
         ibfhi = bfghi1(ishblock)
         idim  = ibfhi - ibflo + 1
         do jshblock = 1, ishblock
            jshlo = shglo1(jshblock)
            jshhi = shghi1(jshblock)
            jbflo = bfglo1(jshblock)
            jbfhi = bfghi1(jshblock)
            jdim  = jbfhi - jbflo + 1
            do kshblock = nshblocks2, 1, -1
               kshlo = shglo2(kshblock)
               kshhi = shghi2(kshblock)
               kbflo = bfglo2(kshblock)
               kbfhi = bfghi2(kshblock)
               kdim  = kbfhi - kbflo + 1
               lblockhi = kshblock
               do lshblock = 1, lblockhi
                  lshlo = shglo2(lshblock)
                  lshhi = shghi2(lshblock)
                  lbflo = bfglo2(lshblock)
                  lbfhi = bfghi2(lshblock)
                  ldim  = lbfhi - lbflo + 1
                  if (next .eq. ijklblock) then
c     
c     Get blocks of the one-particle densities
c     
                     call grad_get_dens_12(
     $                    d1_ij,  d2_kl, blen1, blen2, 
     $                    ibflo, ibfhi, jbflo, jbfhi, 
     $                    kbflo, kbfhi, lbflo, lbfhi, 
     $                    ibflo_prev,jbflo_prev,kbflo_prev,lbflo_prev,
     $                    g_dens1, g_dens2)
c     
c     Build the list of integral shell quartets in NWChem shell labelling
c     
                     oij = ishblock.eq.jshblock
                     okl = kshblock.eq.lshblock
                     oikjl = .false.
c     
                     nq = 0
c     
                     do iish = ishlo, ishhi
                        atoms(1)=bftoat1(bfmap1(shbflo1(iish)))
                        iilo = shbflo1(iish)-ibflo+1
                        iihi = shbfhi1(iish)-ibflo+1
                        ish = shmap1(iish)
                        jshtop = jshhi
                        if (oij) jshtop = iish
                        do jjsh = jshlo, jshtop
                           atoms(2)=bftoat1(bfmap1(shbflo1(jjsh)))
                           jjlo = shbflo1(jjsh)-jbflo+1
                           jjhi = shbfhi1(jjsh)-jbflo+1
                           jsh = shmap1(jjsh)
                           call cneo_schwarz_set(bas1)
                           sij = schwarz_shell(ish,jsh)
                           if (sij*smax*q4max .lt. tol2e) goto 100 ! Next jjsh
                           kshtop = kshhi
                           do kksh = kshlo, kshtop
                              atoms(3)=bftoat2(bfmap2(shbflo2(kksh)))
                              kklo = shbflo2(kksh)-kbflo+1
                              kkhi = shbfhi2(kksh)-kbflo+1
                              ksh = shmap2(kksh)
                              lshtop = lshhi
                              if (okl) lshtop = kksh
                              do llsh = lshlo, lshtop
                                 atoms(4)=bftoat2(bfmap2(shbflo2(llsh)))
                                 lllo = shbflo2(llsh)-lbflo+1
                                 llhi = shbfhi2(llsh)-lbflo+1
                                 lsh = shmap2(llsh)
c     
c     Check on sparsity, symmetry, active list, single center
c     
                                 call cneo_schwarz_set(bas2)
                                 sijkl = sij*schwarz_shell(ksh,lsh)
c     
                                 if (sijkl*q4max .lt. tol2e) goto 200 ! next llsh
                                 scale = 1.0d0
c     
c     Scale according to permutation symmetry of the shell labels
c     
                                 if (.not.oij .or. ish.ne.jsh) 
     $                                scale = scale*2.0d0
                                 if (.not.okl .or. ksh.ne.lsh) 
     $                                scale = scale*2.0d0
c     
c     Make the twopdm just for the current shell block.  Note that
c     pdm2d will have holes with junk in it for quartets we did not
c     compute ... cannot zero this since it will be too expensive.
c     
                                 call grad_make_twopdm_12(
     $                                d1_ij,  d2_kl,  
     $                                blen1, blen2, 
     $                                iilo, jjlo, kklo, lllo,
     $                                iihi, jjhi, kkhi, llhi,
     $                                pdm2d, psum, jfac)
c     
c     Final screening including magnitude of the density
c     
                                 sijkl = sijkl*psum*scale
                                 if (sijkl .gt. tol2e) then
                                    nq = nq + 1
                                    q4(nq) = scale
                                    list(nq,1) = ish
                                    list(nq,2) = jsh
                                    list(nq,3) = ksh
                                    list(nq,4) = lsh
                                 endif
c     
 200                          end do ! next llsh
                           end do
 100                    end do  ! next jjsh
                     end do

*                     numq = numq + nq

c     
                     if (.not. intbd_init4c(
     $                    bas1, list(1,1), list(1,2),
     $                    bas2, list(1,3), list(1,4),
     $                    nq, q4, .true., lscr, scr, leri, 
     $                    block_eff)) call errquit('grad2:txs init?',nq,
     &       INT_ERR)
c     
 10                  omore = intbd_2e4c(
     $                    bas1, list(1,1), list(1,2),
     $                    bas2, list(1,3), list(1,4),
     $                    nq, q4, .true., tol2e, .false.,
     $                    labels(1,1),labels(1,2), 
     $                    labels(1,3), labels(1,4), 
     $                    eri, leri, nint, lscr, scr)
c     
*     write(6,*) ' nint ', nint, omore
c     
                     if (nint .gt. 0) then
                        do integ = 1, nint
                           i = labels(integ,1)
                           j = labels(integ,2)
                           k = labels(integ,3)
                           l = labels(integ,4)
*     write(6,*) ' b ' ,i,j,k,l
                           atoms(1) = bftoat1(i)
                           atoms(2) = bftoat1(j)
                           atoms(3) = bftoat2(k)
                           atoms(4) = bftoat2(l)
*     write(6,*) ' a ' ,atoms
                           i = rbfmap1(i) - ibflo + 1
                           j = rbfmap1(j) - jbflo + 1
                           k = rbfmap2(k) - kbflo + 1
                           l = rbfmap2(l) - lbflo + 1
*     write(6,*) ' c ', i,j,k,l
                           p = pdm2d(l,k,j,i)
*     write(6,*) ' p ', p
                           do iiat = 1, 4
                              iat = atoms(iiat)
                              frc(1,iat) =  frc(1,iat) + 
     $                             p*eri(1,iiat,integ)
                              frc(2,iat) =  frc(2,iat) + 
     $                             p*eri(2,iiat,integ)
                              frc(3,iat) =  frc(3,iat) + 
     $                             p*eri(3,iiat,integ)
                           end do
                        end do
                     end if
                     if (omore) goto 10
c     
                     next = nxtask(nproc,task_size)
                  end if
c     
                  ijklblock = ijklblock + 1
c     
               end do
            end do
         end do
      end do
c
      call int_acc_std()
c
      next = nxtask(-nproc,task_size)
      call ga_sync()
c     
*      write(6,*) ' numq ', numq
c     
      end

c     Modified from grad_dft.F
c     To calculate elec gradient without the consideration of 
c     quantum nuc.
c     Redirect geom/ao_bas_han to geom_cneo/bas_e
c     Remove geom_destroy/bas_destroy for use in cneo dft gradients

      logical function grad_cneo_elec(rtdb,geom_cneo,bas_e)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "schwarz.fh"
#include "stdio.fh"
#include "cdft.fh"
c
c     dftgrad module.
c
c     Context is '...:dftgrad' --> changed to '...:dft'   JAN
c
c     Assumes DFT has been completed, MO vectors stored
c     and all information is still in the RTDB
c
      integer geom_cneo,bas_e
      logical int_normalize
      external int_normalize
      logical int_norm_2c
      external int_norm_2c
c
      integer rtdb              ! [input] database handle
      integer nbases
      logical converged, status
      integer bases(3), plevel
      character*80 theory
c
c----------------------------------------------------------------------- 
c
c     Push context down to DFT
c
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
      call util_print_push()
      call util_print_rtdb_load(rtdb,'dft')
      call ecce_print_module_entry ('dft')
c
      if (.not. rtdb_cget(rtdb, 'dft:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)
      if(theory .eq. 'dft')then 
         status = rtdb_get(rtdb, 'dft:converged', MT_LOG, 1, converged) 
      else if(theory .eq. 'sodft')then 
         status = rtdb_get(rtdb, 'sodft:converged', MT_LOG,1,converged)
      endif
      if (.not.( status .and. converged ))then
         if ((ga_nodeid() .eq. 0) 
     &      .and. util_print('information', print_none) )then
            write (luout,*)'status: ', status, '   converged: ',
     &                     converged
        endif
        call errquit(
     &       'dft grad_dft: no converged DFT wavefunction available',
     &     0, UNKNOWN_ERR)
      endif 
c
c     Extract high level info from the data-base setting defaults
c
c     load DFT parameters into common; turn off DFT printing
c
      call dft_pstat_init(rtdb)
c
      call util_print_get_level(plevel)
      call util_print_set_level(print_none)
      call dft_rdinput(rtdb)
      
      if (.not. bas_destroy(ao_bas_han)) call errquit
     $      ('grad_cneo_elec:destroy ao_bas_han?',0,0)
      ao_bas_han = bas_e
      
      call util_print_set_level(plevel)
c
      if (XCFIT.and.CDFIT.and.(.not.ADFT))then
         nbases = 3
         bases(1) = AO_bas_han
         bases(2) = CD_bas_han
         bases(3) = XC_bas_han
      elseif (XCFIT.and.CDFIT.and.ADFT)then
         nbases = 2
         bases(1) = AO_bas_han
         bases(2) = CD_bas_han
      elseif((.not.XCFIT).and.CDFIT)then
         nbases = 2
         bases(1) = AO_bas_han
         bases(2) = CD_bas_han
      elseif((.not.CDFIT).and.XCFIT)then
         nbases = 2
         bases(1) = AO_bas_han
         bases(2) = XC_bas_han
      else
         nbases = 1
         bases(1) = AO_bas_han
      endif
      if (.not. geom_destroy(geom)) call errquit
     $      ('grad_cneo_elec:geom_destroy?',0,0)
      geom = geom_cneo
      call dft_inpana(rtdb)
c
c     initialize for schwarz screening
c
      if (nbases .gt. 1) call int_app_set_no_texas(rtdb)
      call int_init(rtdb, 1, ao_bas_han)
      call schwarz_init (geom, ao_bas_han)
      call int_terminate()
c
c     initialize for derivative integrals
c
      call intd_init(rtdb,nbases,bases)
c
      if (ga_nodeid() .eq. 0)then
         if (util_print('information',print_default) )then
            write(luout,*)
            write(luout,*)
            call util_print_centered(luout,
     &          'NWChem DFT Gradient Module',40,.true.)
            write(luout,*)
            write(luout,*)
            if (title .ne. ' ')then
               call util_print_centered(luout, title, 40, .false.)
               write(luout,*)
               write(luout,*)
            endif

            if (ipol .eq. 1)then
               write(luout,1) rcharge, 'closed shell'
            else 
               write(luout,1) rcharge, 'open shell'
            endif 
            call util_flush(luout)

         endif

 1       format(/
     &        '  charge          = ', f6.2/
     &        '  wavefunction    = ', a/)
         if (util_print('debug',print_debug) )then
            if (.not. geom_print(geom)) 
     &         call errquit('grad_dft: geom_print ?',0, GEOM_ERR)
            if (.not. bas_print(ao_bas_han)) 
     &         call errquit('grad_dft: bas_print ?',0, BASIS_ERR)
         endif
      endif
c
c     Compute the gradients arising from CD fit and XC
c
      if(theory .eq. 'dft')call dft_gradients(rtdb)
      if(theory .eq. 'sodft')call dft_gradients_so(rtdb)
c
      call schwarz_tidy ()
c
      if(CDFIT)then
        if (.not. bas_destroy(cd_bas_han))
     &     call errquit('grad_dft:not able to destroy CD_bas:',86,
     &       BASIS_ERR)
      endif
      if(XCFIT.and.(.not.ADFT))then
        if (.not.bas_destroy(XC_bas_han))
     &     call errquit('grad_dft:not able to destroy XC_bas:',86,
     &       BASIS_ERR)
      endif

      call intd_terminate()
      if (nbases .gt. 1) call int_app_unset_no_texas(rtdb)
c
c     Compute the rest of the gradients (1-e and 2-e)
c
      if(theory .eq. 'dft')call grad_force(rtdb, ao_bas_han, geom)
      if(theory .eq. 'sodft')call grad_force_so(rtdb, ao_bas_han, geom)
c
      call dft_pstat_print
c
c     terminate integral scope
c
      call ecce_print_module_exit ('dft','ok')
      call util_print_pop
      grad_cneo_elec = .true.
c
      return
      end

      subroutine grad_get_dens_12(d1_ij, d2_kl, blen1, blen2,
     $         ilo, ihi, jlo, jhi, klo, khi, llo, lhi, 
     $         ilo_prev, jlo_prev, klo_prev, llo_prev,
     $         g_dens1, g_dens2)
      implicit none 
#include "global.fh"

      double precision d1_ij(*),d2_kl(*)
      integer blen1, blen2
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      integer ilo_prev, jlo_prev, klo_prev, llo_prev
      integer g_dens1, g_dens2

      logical ogetij
      logical ogetkl

      ogetij = ilo.ne.ilo_prev .or. jlo.ne.jlo_prev
      ogetkl = klo.ne.klo_prev .or. llo.ne.llo_prev

      ilo_prev = ilo
      jlo_prev = jlo
      klo_prev = klo
      llo_prev = llo

      if (ogetij) then
         call ga_get(g_dens1, ilo, ihi, jlo, jhi,d1_ij,blen1)
      endif
      if (ogetkl) then
         call ga_get(g_dens2, klo, khi, llo, lhi,d2_kl,blen2)
      endif

      end

      subroutine grad_make_twopdm_12(d1_ij, d2_kl, blen1, blen2,
     $         ilo, jlo, klo, llo, ihi, jhi, khi, lhi, 
     $         pdm2d, psum, jfac)
      implicit none 
#include "errquit.fh"
      double precision d1_ij(blen1,blen1), d2_kl(blen2,blen2)
      integer blen1, blen2
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      double precision pdm2d(blen2,blen2,blen1,blen1)
      double precision psum, jfac

      double precision den2
      integer i, j, k, l

      psum = 0.0d0
      do i = ilo, ihi
         do j = jlo, jhi
            do k = klo, khi
               do l = llo, lhi
                  den2 = jfac*(d1_ij(i,j)*d2_kl(k,l))
                  pdm2d(l,k,j,i) = den2
                  psum = psum+den2**2
               enddo
            enddo
         enddo
      enddo
      psum = sqrt(psum)

      end

      subroutine cneo_dens_from_file(movecs,geom,basis,g_dens,
     $            nbf,nclosed,nopen,scftype)
      implicit none 
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
      character*(*) movecs
      integer geom,basis,g_dens, nclosed, nbf,nopen
      character*(*) scftype

      integer g_vecs, g_vecs2, g_dens2
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical movecs_read
      external movecs_read

      double precision occ(nbf)
      double precision eval(nbf)
      logical status
      integer g_odens
      double precision docca, doccb
      integer nocc, nocca, noccb, i

      g_vecs = ga_create_atom_blocked (geom, basis, 'mo vectors' )
      if (.not. movecs_read(movecs,1,occ,eval,g_vecs)) call errquit
     $      ('cneo_grad_dens:movecs_read?',0,0)

      docca=0.0d0
      doccb=0.0d0
      do i=1, nbf
         docca = docca+occ(i)
      enddo
      status = .true.
      if (scftype .eq. 'RHF') then
         nocc = nclosed
         call ga_matmul_patch('n', 't', 2.0d0, 0.0d0,
     $        g_vecs,    1, nbf, 1, nocc,
     $        g_vecs,    1, nocc, 1, nbf,
     $        g_dens, 1, nbf, 1, nbf)
         call ga_symmetrize(g_dens)

      else if (scftype .eq. 'UHF') then
         if (.not. ga_duplicate(g_dens,g_dens2,'tmp uhf dens'))
     $      call errquit('cneo_grad_dens:ga_duplicate?',0,0)
         g_vecs2 = ga_create_atom_blocked (geom, basis, 
     $                                    'uhf mo vectors' )

         occ = 0.0d0
         eval = 0.0d0
         if (.not. movecs_read(movecs, 2, occ, eval,
     $        g_vecs2)) call errquit
     $        ('grad_dens: could not read beta UHF MOs',0, DISK_ERR)
         do i = 1, nbf
           doccb = doccb + occ(i)
         enddo

         if (docca.ge.doccb) then
            nocca = nclosed + nopen ! == nalpha
            noccb = nclosed         ! nbeta
         else
            nocca = nclosed         ! == nalpha
            noccb = nclosed + nopen ! nbeta
         endif
         call ga_matmul_patch('n', 't', 1.0d0, 0.0d0,
     $       g_vecs,    1, nbf,           1, nocca,
     $       g_vecs,    1, nocca, 1, nbf,
     $       g_dens, 1, nbf,           1, nbf)
         call ga_matmul_patch('n', 't', 1.0d0, 0.0d0,
     $       g_vecs2,      1, nbf,     1, noccb,
     $       g_vecs2,      1, noccb, 1, nbf,
     $       g_dens2,       1, nbf,     1, nbf)

         call ga_dadd(1.0d0, g_dens, 1.0d0, g_dens2, g_dens)

         status = status .and. ga_destroy(g_dens2)
         status = status .and. ga_destroy(g_vecs2)
      else
         call errquit('cneo_grad_dens:invalid scftype',0,0)
      endif

      status = status .and. ga_destroy(g_vecs)
      if (.not. status) call errquit
     $         ('cneo_grad_dens:ga_destroy?',0,0)

      end