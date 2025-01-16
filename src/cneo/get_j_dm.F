      subroutine get_j_dm(geom,bas1,bas2,jfac,
     $     tol2e,oskel,g_dens,g_fock,inttype,
     $     idx1, idx2, ofock)
c     This subroutine implicitly assumed that the input
c     density matrix is symmetric, so as the output fock
c     matrix.
c     As a result, this subroutine cannot be utilized 
c     in cneo TDDFT grad for the calculation of en response.
c     Similar subroutine is implemented for TDDFT gradient
c     calculation
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "geom.fh"
#include "cfock.fh"
#include "util.fh"
#include "stdio.fh"

      integer geom, bas1, bas2
      double precision jfac, tol2e
      logical oskel
      integer g_dens, g_fock

      logical int2e_file_open, int2e_file_rewind, int2e_file_close
      external int2e_file_open, int2e_file_rewind, int2e_file_close
      logical int2e_test_mem,int2e_buf_write
      external int2e_test_mem,int2e_buf_write

      integer nbf1, nbf2, nsh1, nsh2
      integer natom
      logical status

      integer l_q4,k_q4,l_ijkl,k_ijkl
      integer l_plist1,k_plist1,l_ps1,k_ps1,l_pinfo1,k_pinfo1
      integer l_plist2,k_plist2,l_ps2,k_ps2,l_pinfo2,k_pinfo2
      integer tdim
      parameter (tdim = 100000)
      integer memlist(15),lmemlist

      character*(*) inttype
      integer l_fock, k_fock, l_dens, k_dens, l_rdens, k_rdens
      integer l_dens_tmp,k_dens_tmp
      double precision tol2e_local
      integer idx1, idx2, iatom1, iatom2
      logical ofock
      character*5 tmp_idx1,tmp_idx2 
      character*20 tmp_filename

      write(tmp_idx1, '(I0)') idx1
      write(tmp_idx2, '(I0)') idx2
      if (inttype .eq. 'aoints') then
         write(6,*) 'warning: 2e4c int under lower symmetry'
         call util_file_name('aoints',.true.,.false.,int2efilename)
      else if (inttype .eq. 'neints') then
         if (ofock) then
            if (.not. oreadfile)
     $         call errquit('cneo:get_j_dm: init enints first',0,0)
            tmp_filename = 'enints'//trim(tmp_idx1)
         else 
            tmp_filename = 'neints'//trim(tmp_idx1)
         endif
      else if (inttype .eq. 'enints') then
         tmp_filename = 'enints'//trim(tmp_idx1)
      else if (inttype .eq. 'nnints') then
         if (idx1 .lt. idx2) then
            tmp_filename = 'nnints'//trim(tmp_idx1)//trim(tmp_idx2)
         else if (idx1 .gt. idx2) then
            if (.not. oreadfile)
     $            call errquit('get_j_dm: init nnints',0,0)
            tmp_filename = 'nnints'//trim(tmp_idx2)//trim(tmp_idx1)
         else 
            call errquit('get_j_dm: same idx1, idx2?',0,0)
         endif
      else
         call errquit('cneo:get_j_dm: int type not supported',0,0)
      endif

      call util_file_name(tmp_filename,.true.,.false.,int2efilename)

      ! print *, 'cneo:get_j_dm:int2efilename:', int2efilename, inttype
      maxquartet = 10000
      maxeri = 1296*100

      tol2e_local = min(tol2e,1d-8)

      if (oreadfile) then
         if (.not. int2e_file_open(int2efilename,0,0,tol2e,.true.))
     $         call errquit('cneo:eri_dm:open int2e file',0,0)
      else
         if (.not. int2e_file_open(int2efilename,0,0,tol2e,.false.))
     $         call errquit('cneo:eri_dm:open int2e file',0,0)
      end if

      status = int2e_file_rewind()

      status = status .and. bas_numbf(bas1,nbf1)
      status = status .and. bas_numbf(bas2,nbf2)
      status = status .and. bas_numcont(bas1,nsh1)
      status = status .and. bas_numcont(bas2,nsh2)
      status = status .and. geom_ncent(geom,natom)

      if (.not. status) call errquit('get_j_dm: geom/bas',0,0)

      lmemlist = 0
      status = .true.
      status = status .and. ma_push_get(MT_DBL, nbf2*nbf2, 
     $     'en_rep: density', l_dens, k_dens)
      if (status) then
        lmemlist = lmemlist + 1
        memlist(lmemlist) = l_dens
      endif
      status = status .and. ma_push_get(MT_DBL, nbf1*nbf1, 
     $     'en_rep: fock',    l_fock, k_fock)
      if (status) then
        lmemlist = lmemlist + 1
        memlist(lmemlist) = l_fock
      endif
      status = status .and. ma_push_get(MT_DBL, nsh2*nsh2, 
     $     'en_rep: fock',    l_rdens, k_rdens)
      if (status) then
        lmemlist = lmemlist + 1
        memlist(lmemlist) = l_rdens
      endif
      status = status .and. ma_push_get(MT_DBL, maxquartet,
     $     'en_rep: q4',    l_q4, k_q4)
      if (status) then
        lmemlist = lmemlist + 1
        memlist(lmemlist) = l_q4
      endif
      status = status .and. ma_push_get(MT_INT, 4*maxquartet,
     $     'en_rep: ijkl',    l_ijkl, k_ijkl)
      if (status) then
        lmemlist = lmemlist + 1
        memlist(lmemlist) = l_ijkl
      endif
      status = status .and. ma_push_get(MT_INT, nbf1*(nbf1+1),
     $     'en_rep: plist1',    l_plist1, k_plist1)
      if (status) then
        lmemlist = lmemlist + 1
        memlist(lmemlist) = l_plist1
      endif
      status = status .and. ma_push_get(MT_DBL, nbf1*(nbf1+1)/2,
     $     'en_rep: plist1',    l_ps1, k_ps1)
      if (status) then
        lmemlist = lmemlist + 1
        memlist(lmemlist) = l_ps1
      endif
      status = status .and. ma_push_get(MT_INT, 2*tdim,
     $     'en_rep: pinfo1',    l_pinfo1, k_pinfo1)
      if (status) then
         lmemlist = lmemlist + 1
         memlist(lmemlist) = l_pinfo1
      endif
c     For bas2
      status = status .and. ma_push_get(MT_INT, nbf2*(nbf2+1),
     $     'en_rep: plist2',    l_plist2, k_plist2)
      if (status) then
        lmemlist = lmemlist + 1
        memlist(lmemlist) = l_plist2
      endif
      status = status .and. ma_push_get(MT_DBL, nbf2*(nbf2+1)/2,
     $     'en_rep: plist2',    l_ps2, k_ps2)
      if (status) then
        lmemlist = lmemlist + 1
        memlist(lmemlist) = l_ps2
      endif
      status = status .and. ma_push_get(MT_INT, 2*tdim,
     $     'en_rep: pinfo2',    l_pinfo2, k_pinfo2)
      if (status) then
         lmemlist = lmemlist + 1
         memlist(lmemlist) = l_pinfo2
      endif
      ! add error info later (line 143 in ao_replicated)

      call ga_sync()
      if (ga_nodeid() .eq. 0) then
         status = ma_push_get(MT_DBL, nbf2*nbf2, 
     $     'en_rep: density tmp', l_dens_tmp, k_dens_tmp)
         if (.not. status) call errquit('cneo:fock:tmp dens',0,0)
         call ga_get(g_dens, 1, nbf2, 1, nbf2, 
     $           dbl_mb(k_dens_tmp), nbf2)
         call scf_dens_reorder(1, nbf2,
     $        dbl_mb(k_dens_tmp), dbl_mb(k_dens))
         if(.not. ma_pop_stack(l_dens_tmp))
     $         call errquit('cneo:fock:free tmp dens',0,0)
      end if

      call ga_sync()
      call ga_brdcst(373, dbl_mb(k_dens), mdtob(nbf2*nbf2), 0)
      call ga_sync()
      call dfill(nbf1*nbf1, 0.0d0, dbl_mb(k_fock), 1)

      if (util_print('semi_direct',print_debug) .and.
     $     ga_nodeid().eq.0) then
         write(6,*) ' oread, owrite, incmax, incmin, incomplete ',
     $        oreadfile, owritefile, incmax, incmin, incomplete
         call util_flush(6)
      endif
c
      if (oreadfile) then
         if (util_print('semi_direct',print_debug) .and.
     $        ga_nodeid().eq.0) then
            write(6,*) ' Building from file ', tol2e_local
            call util_flush(6)
         endif
         call fock_en_rep_from_file(geom, bas1, bas2, nbf1, nbf2,
     $             jfac, tol2e_local, dbl_mb(k_dens), 
     $             dbl_mb(k_fock), inttype, idx1, idx2)
      else

         call fock_en_rep(geom, bas1,bas2, nbf1,nbf2, nsh1,nsh2,
     $     natom, jfac, tol2e_local, oskel,
     $     dbl_mb(k_dens), dbl_mb(k_fock), dbl_mb(k_rdens),
     $     int_mb(k_plist1), int_mb(k_plist2),
     $     dbl_mb(k_ps1), dbl_mb(k_ps2), maxquartet,
     $     dbl_mb(k_q4), int_mb(k_ijkl), tdim,
     $     int_mb(k_pinfo1), int_mb(k_pinfo2))
      endif

      call ga_sync()
      call ga_dgop(374, dbl_mb(k_fock), nbf1*nbf1, '+')
      call ga_sync()
      
      call ga_put(g_fock, 1, nbf1, 1, nbf1,
     $         dbl_mb(k_fock),nbf1)
      
      status = ma_pop_stack(l_pinfo2)
      status = ma_pop_stack(l_ps2) .and. status
      status = ma_pop_stack(l_plist2) .and. status
      status = ma_pop_stack(l_pinfo1) .and. status
      status = ma_pop_stack(l_ps1) .and. status
      status = ma_pop_stack(l_plist1) .and. status
      status = ma_pop_stack(l_ijkl)  .and. status
      status = ma_pop_stack(l_q4)    .and. status
      status = ma_pop_stack(l_rdens) .and. status
      status = ma_pop_stack(l_fock)  .and. status
      status = ma_pop_stack(l_dens)  .and. status
      if (.not. status) call errquit('en_eri: ma_pop?', 0, MA_ERR)

      call ga_sync()
      call ga_dscal(g_fock, 4.0d0)
      if (inttype .ne. 'aoints') then
         call ga_symmetrize(g_fock)
      endif

      if (owritefile) then
         if (.not. int2e_buf_write(.true.))
     $         call errquit('cneo:get_j_dm: write file',0,0)
         owritefile = .false.
         oreadfile = .true.
      end if

      if (ofock) then
         if (.not. int2e_file_close(.true.))
     $         call errquit('cneo:get_j_dm: close file',0,0)
      else
         if (.not. int2e_file_close(.false.))
     $         call errquit('cneo:get_j_dm: close file',0,0)
      endif

      ! call ga_print(g_fock)

      end

      subroutine fock_en_rep(geom, bas1, bas2, nbf1,nbf2, nsh1, nsh2,
     $      natom, jfac, tol2e, oskel, dens, fock, rdens,
     $      plist1, plist2, ps1, ps2, qdim, q4, ijkl, tdim, 
     $      pinfo1, pinfo2)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "geom.fh"
#include "bas.fh"
#include "global.fh"
#include "cfock.fh"
#include "util.fh"
#include "mafdecls.fh"

      integer geom, bas1, bas2
      ! integer nfock
      integer nbf1, nbf2
      integer nsh1, nsh2
      integer natom
      double precision jfac
      double precision tol2e
      logical oskel
      double precision dens(nbf2*nbf2)
      double precision fock(nbf1*nbf1)
      double precision rdens(nsh2,nsh2)
      integer plist1(2,*), plist2(2,*)        ! Pairs sorted by type
      double precision ps1(*), ps2(*)    ! Schwarz screening for pair
      integer qdim
      double precision q4(qdim) ! Q4 factors for each quartet
      integer ijkl(qdim, 4)     ! Indices for quartet
      integer tdim              ! Max no. of pair types (guess)
      integer pinfo1(2,tdim), pinfo2(2,tdim)     ! Start and end of each pair type
      integer ntp               ! No. of tasks per processor
c     
      integer natom_unique
      integer nshpair1, nshpair2
      integer ntype1, ntype2             ! No. of pair types
      integer max_p             ! Max pairs per block
      integer chunk             ! Chunking factor for tasks
c     
      integer ijklt, ijt, klt, ij, kl, ish, jsh, ksh, lsh
      integer  iat, jat, kat, lat, imat
      double precision q2, qq4, qjunk, sijkl, denmax, dijkl
      double precision dentol, integ_acc
      character*80 buff
c
      double precision jfacmax
      double precision nqraw, nqsch, nqsym
c
      integer nq                ! No. of quartets generated
c
      integer ntreq
      double precision neri_direct, neri_cached
      double precision nqdone, neri_done
      common/fred/nqdone, neri_done, neri_direct, neri_cached, ntreq
c
      integer cnmax
      parameter (cnmax=nw_max_shells) ! Max no. of shells 
      integer cntoce1(cnmax), cntoce2(cnmax)          ! Map from contractions to centers
      integer cntobfr1(2,cnmax), cntobfr2(2,cnmax)       ! Map from contractions to bf range
c
      integer nops              ! No. of group operators = max value of q4
      double precision dnops, dentol_pairs
c
      integer next, current, nproc, incomplete_top
c
      logical oprint, osavewritefile, odotask
      logical odonea
c     
      logical int2e_set_bf_range
      external int2e_set_bf_range
      integer nxtask
      external nxtask

      next = 1
      incomplete_top = 0

      ! nops = sym_number_ops(geom)+1
      ! dnops = nops
      dnops = 1
      nops = 1

      if (.not. int2e_set_bf_range(
     $      1,nbf1,1,nbf1,1,nbf2,1,nbf2)) call errquit
     $      ('fock_en_rep: set range failed',0,0)
      incomplete = -1

      if (.not. geom_ncent_unique(geom, natom_unique)) call errquit
     $     ('fock_2e_rep: bad geom handle?',0, GEOM_ERR)
      if (nsh1 .gt. cnmax)   call errquit('fock_en_rep: cnmax',nsh1,
     &       UNKNOWN_ERR)
      if (nsh2 .gt. cnmax)   call errquit('fock_en_rep: cnmax',nsh2,
     &       UNKNOWN_ERR)

      do ish = 1, nsh1
         if (.not. bas_cn2ce(bas1, ish, cntoce1(ish)))
     $        call errquit('ao_rep: bad basis', 0, BASIS_ERR)
      end do
      do ish = 1, nsh2
         if (.not. bas_cn2ce(bas2, ish, cntoce2(ish)))
     $        call errquit('ao_rep: bad basis', 0, BASIS_ERR)
      end do
      do ish = 1, nsh1
         if (.not. bas_cn2bfr(bas1,ish,cntobfr1(1,ish),cntobfr1(2,ish)))
     $        call errquit('ao_rep: bad basis', 0, BASIS_ERR)
      end do
      do ish = 1, nsh2
         if (.not. bas_cn2bfr(bas2,ish,cntobfr2(1,ish),cntobfr2(2,ish)))
     $        call errquit('ao_rep: bad basis', 0, BASIS_ERR)
      end do

      call dfill(nsh2*nsh2,0.0d0,rdens,1)
      jfacmax = max(0.0d0, abs(jfac))
      call util_mat_reduce(nbf2, nsh2, cntobfr2, dens(1), 1, 
     $         rdens,'rms')

      denmax = 0.0d0
      do ish=1,nsh2
         do jsh=1,nsh2
            denmax = max(denmax,rdens(jsh,ish))
         end do
      end do

      if(abs(denmax).gt.0) then
         dentol = min(dentolmax,tol2e/denmax,tol2e/denmax**2) ! Threshold to screen integs only
      else
         dentol = tol2e
      endif

      if (intacc .eq. 0.0d0) then
cedo         integ_acc = 1d-15 ! Not used
         integ_acc = 1d-11 !   This seems fine with txs EA 1/2003
         continue
*         integ_acc = min(1d-15,max(0.01d0*dentol,1d-30)) ! Variable
*         call int_acc_set(integ_acc)
      else
         integ_acc = intacc     ! User controlled
         call int_acc_set(integ_acc)
      end if
      if (util_print('en_rep_info', print_debug) .and.
     $       ga_nodeid().eq.0) then
         write(6,111) tol2e, denmax, dentol, integ_acc, odensityscreen,
     $        oreadfile, owritefile
 111     format(1x,1p,' tol2e=',d8.1,' denmax=',d8.1,
     $        ' dentol=',d8.1,' integ_acc=',d8.1,
     $        ' odenscr=',l1,' oread=',l1,' owrite=',l1)
      endif

      ntp = 10
      nshpair1 = nsh1*(nsh1+1)/2
      nshpair2 = nsh2*(nsh2+1)/2
      if (ga_nnodes() .gt. 1) then
         max_p = nshpair1*nshpair2
         max_p = max_p/sqrt(dble(ntp)*ga_nnodes())
         max_p = min(nint(sqrt(dble(nops*qdim))),max_p)
         max_p = max(10,max_p)
      else
         max_p = int(sqrt(dble(nops*qdim)))
      end if

      dentol_pairs = dentol/dnops

 10   call cneo_schwarz_set(bas1)
      call fock_pairs(bas1, max_p, .true., 1, natom, 1, natom, 
     $     dentol_pairs, plist1, ps1, pinfo1, tdim, ntype1)
      call cneo_schwarz_set(bas2)
      call fock_pairs(bas2, max_p, .true., 1, natom, 1, natom, 
     $     dentol_pairs, plist2, ps2, pinfo2, tdim, ntype2)
      if (ga_nodeid() .eq. 0) call util_flush(6)

      chunk = max(1, (ntype1*ntype2)/(10*ga_nnodes()))
      if (max_p .lt. nint(sqrt(dble(qdim))) .and.
     $     chunk .gt. 4) then
         max_p = min(nint(sqrt(dble(qdim))),2*max_p)
         ! if (oprint .and. ga_nodeid().eq.0) write(6,*) ' A ', max_p
         odonea = .true.
         goto 10
      else if ((ntype1*ntype2)/ga_nnodes().lt.ntp  .and.
     $        max_p.gt.10 .and. (.not. odonea)) then
         max_p = max(10,max_p/2)
         ! if (oprint .and. ga_nodeid().eq.0) write(6,*) ' B ', max_p
         goto 10
      end if
      if ((util_print('semi_direct',print_debug))
     $     .and. ga_nodeid().eq.0) then
         write(6,*) ' max_p ', max_p, ' ntype1 ', ntype1,
     $     ' ntype2 ',ntype2,' chunk ',chunk
         call util_flush(6)
      end if

      
      nproc = ga_nnodes()
      ! next = nxtask(nproc, chunk)

      next = 0

      current = next
      nq = 0
c
      nqraw = 0
      nqsch = 0
      nqsym = 0
c
      nqdone = 0
      neri_done = 0
      ntreq  = 1
      neri_direct = 0
      neri_cached = 0
c
      do ijklt = ntype1+ntype2, 2, -1
         do ijt = min(ntype1,ijklt-1), max(1,ijklt-ntype2), -1
            klt = ijklt - ijt
            odotask = next.eq.current
            if (odotask) then
c
               call int2e_file_record_position()
               do ij = pinfo1(1,ijt), pinfo1(2,ijt)
                  ish = plist1(1,ij)
                  jsh = plist1(2,ij)
                  iat = cntoce1(ish)
                  jat = cntoce1(jsh)
                  q2  = 1.0d0
                  if (ish .eq. jsh) q2 = 0.5d0
                  do kl = pinfo2(1,klt),pinfo2(2,klt)
                     qq4 = q2
                     ksh = plist2(1,kl)
                     lsh = plist2(2,kl)
                     kat = cntoce2(ksh)
                     lat = cntoce2(lsh)

                     if (ksh.eq.lsh) qq4 = qq4*0.5d0
c
                     nqraw = nqraw + 1.0d0
c     schwarz screening
                     sijkl = ps1(ij)*ps2(kl)
                     if (sijkl*1*qq4 .lt. 0) goto 200 ! Next kl
                     nqsch = nqsch + 1.0d0
                     nqsym = nqsym + 1.0d0
c     save it for later
                     nq = nq + 1
                     ijkl(nq,1) = ish
                     ijkl(nq,2) = jsh
                     ijkl(nq,3) = ksh
                     ijkl(nq,4) = lsh
                     q4(nq)     = qq4
                     if (nq .eq. qdim) then
                        nqdone = nqdone + nq
                        call fock_rep_en_txs(bas1,bas2,nbf1,nbf2,
     $                       jfac, dentol, dens, fock,
     $                       nq, ijkl(1,1), ijkl(1,2), ijkl(1,3),
     $                       ijkl(1,4), q4, current)
                        nq = 0
                     end if
 200              end do        ! Next KL pair
 100           end do           ! Next IJ pair
c     
c     Empty out interactions for this type
c     
               if (nq .gt. 0) then
                  nqdone = nqdone + nq
                  call fock_rep_en_txs(bas1,bas2,nbf1,nbf2,
     $                       jfac, dentol, dens, fock, 
     $                       nq, ijkl(1,1), ijkl(1,2), ijkl(1,3),
     $                       ijkl(1,4), q4, current)
                  nq = 0
               end if
c
c     Get next task ... messy because of semidirect
c
               if (oreadfile) then ! Semi-direct
                  continue
               else
                  ! next = nxtask(nproc,chunk)
                  next = next + 1
                  ntreq = ntreq + 1
               endif
            end if
            current = current + 1
         end do                 ! Next IJ type
      end do                    ! Next IJKL type

      call ga_sync()
      end

      subroutine fock_rep_en_txs(bas1, bas2, nbf1, nbf2, jfac, 
     $      tol2e, dens, fock, nq, iq, jq, kq, lq, q4, current)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cfock.fh"
#include "util.fh"
      integer bas1, bas2
      integer nbf1, nbf2
      double precision jfac
      double precision tol2e
      double precision dens(nbf2*nbf2)
      double precision fock(nbf1*nbf1)
      integer nq
      integer iq(nq), jq(nq), kq(nq), lq(nq)
      double precision q4(nq)
      integer current           ! Current task number for semi-direct
c     
c     Given a list of integral quartets, notionally all of the
c     same type, compute the integrals and add them into the
c     fock matrix, optionally caching them on disk/memory.
c     
      integer lenscr, leneri, lenlab
      integer l_eri, k_eri, l_lab, k_lab, l_scr, k_scr
      integer k_i, k_j, k_k, k_l 
*      integer kkk
c     
      integer neri, ntries
      double precision block_eff
      logical more
c
      integer maxg
c
      integer ntreq
      double precision neri_direct, neri_cached
      double precision nqdone, neri_done
      common/fred/nqdone, neri_done, neri_direct, neri_cached,ntreq
c     
      logical int2e_file_write, intb_2e4c, modified_intb_init4c
      external int2e_file_write, intb_2e4c, modified_intb_init4c
c     
      if (nq .le. 0) return
c
c     Figure out dimensions for scratch and eri buffer space.
c     The eri buffer space is defaulted in scf_get_info to some
c     suitable non-zero number.  The scratch space defaults to zero
c     in which case int_mem_2e4c is used, otherwise the user input
c     is taken.
c
      call intb_mem_2e4c(maxg, lenscr) ! blocking algorithm
      if (maxscr .ne. 0) then
         lenscr = maxscr
      end if
      leneri = max(maxeri,maxg)
      lenlab = 4*leneri
c
      if (util_print('fock_rep_en_txs', print_debug)) then
         write(6,11) nq, lenscr, leneri
 11      format(' nq=',i6,', lenscr=',i8,', leneri=',i8)
         call util_flush(6)
      end if
c
c     Allow two goes at reducing memory before fail completely
c
c     Ooops ... cannot reduced the texas scratch memory.  Only
c     try reducing the ERI and LAB memory
c
      ntries = 1                ! No. of tries at allocating memory
 12   if (ntries.eq.1) then
         continue
      else if (ntries.le.3) then
****         lenscr = lenscr/2
         leneri = leneri/2
         lenlab = 4*leneri
      else
         call errquit('fock_rep_txs: insufficient memory ',
     $        lenscr+leneri+lenlab, MEM_ERR)
      endif
      ntries = ntries + 1
c 
      if (.not. ma_push_get(mt_dbl, leneri, 'eri', l_eri, k_eri)) then
         goto 12
      endif
      if (.not. ma_push_get(mt_dbl, lenscr, 'scr', l_scr, k_scr)) then
         if (.not. ma_pop_stack(l_eri))
     $        call errquit('fock_rep_txs: ma corrupt',0, MA_ERR)
         goto 12
      endif
      if (.not. ma_push_get(mt_int, lenlab, 'lab', l_lab, k_lab)) then
         if (.not. ma_pop_stack(l_scr))
     $        call errquit('fock_rep_txs: ma corrupt',0, MA_ERR)
         if (.not. ma_pop_stack(l_eri))
     $        call errquit('fock_rep_txs: ma corrupt',0, MA_ERR)
         goto 12
      endif
c     
      k_i = k_lab
      k_j = k_i   + leneri
      k_k = k_j   + leneri
      k_l = k_k   + leneri
c     
      if (.not. modified_intb_init4c(bas1, iq, jq, bas2, kq, lq,
     $     nq, q4, .true., lenscr, dbl_mb(k_scr), leneri, block_eff))
     $     call errquit('intb_init?',nq, UNKNOWN_ERR)
c     
 50   continue
c
      more = intb_2e4c(bas1, iq, jq, bas2, kq, lq,
     $     nq, q4, .true., tol2e, .false.,
     $     int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $     dbl_mb(k_eri), leneri, neri, lenscr, dbl_mb(k_scr))
*
*
c     
      if (neri .gt. 0) then
c
         neri_done = neri_done + neri
c     
c     Cache integrals if desired
c     
         if (.not. int2e_file_write(neri,
     $           int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $           dbl_mb(k_eri))) then
            call int2e_file_reposition_truncate()
            owritefile = .false.
            incomplete = current
            call ga_sync()
         else
            neri_cached = neri_done
         endif
c     
c     
         call fock_en_rep_1_label(nbf1, nbf2, jfac, tol2e, neri,
     $        int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $        dbl_mb(k_eri), dens, fock)
      endif
c     
      if (more) then
         if (util_print('fock_rep_en_txs', print_never)) then
            write(6,*) ' Texas has split the request '
         end if
         goto 50
      end if
c     
 100  if (.not. ma_pop_stack(l_lab))
     $     call errquit('fock_rep_txs: cannot free lab',lenlab, MA_ERR)
      if (.not. ma_pop_stack(l_scr))
     $     call errquit('fock_rep_txs: cannot free scr',lenlab, MA_ERR)
      if (.not. ma_pop_stack(l_eri))
     $     call errquit('fock_rep_txs: cannot free eri',lenlab, MA_ERR)
c     
      end

      subroutine fock_en_rep_1_label(nbf1, nbf2, jfac, tol2e,
     $      neri, ilab, jlab, klab, llab, eri, dens, fock)
      implicit none
#include "nwc_const.fh"

      integer nbf1, nbf2
      double precision jfac, tol2e
      integer neri
      integer ilab(neri), jlab(neri),klab(neri),llab(neri)
      double precision eri(neri)
      double precision fock(*), dens(*)

      integer i, j, k, l, ind
      integer ij_fock, lk_dens, ik, li, lj, jk
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
c
            i = ilab(ind)
            j = jlab(ind)
            k = klab(ind)
            l = llab(ind)
cc
            ij_fock = i + im1(j)
            lk_dens = l + im2(k)
c
            fock(ij_fock) = fock(ij_fock) + gj*dens(lk_dens)
         end if
      end do
c     
      end