      subroutine cneo_nuc_init(rtdb, geom)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "case.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "ccneo.fh"


      integer rtdb  ![input]
      integer geom  ![input]
      character*255 tmp
      character tmp_smb
      character quantum_nuc_smb
      integer nat
      logical status
      integer i,ii,k

      status = geom_ncent(geom,nat)
      tmp = 'cneo:quantum_nuc:'
      k = inp_strlen(tmp)+1

      tmp(k:) = ' '
      tmp(k:) = 'nuc_num'
      if (rtdb_get(rtdb,tmp,mt_int,1,nuc_num)) then
        tmp(k:) = ' '
        tmp(k:) = 'idx'
        status = rtdb_get(rtdb,tmp,mt_int,nuc_num,nuc_cent)
      else
        tmp(k:) = ' '
        tmp(k:) = 'symbol'
        if (rtdb_cget(rtdb,tmp,1,quantum_nuc_smb)) then
 110      ii = 0
          nuc_cent = 0
          do 100 i=1,nat
            tmp_smb = tags(i,geom)
            if (tmp_smb .eq. quantum_nuc_smb) then
              ii = ii+1
              nuc_cent(ii) = i
            end if
 100      continue
          tmp(k:) = ' '
          tmp(k:) = 'nuc_num'
          nuc_num = ii
          status = rtdb_put(rtdb,tmp,mt_int,1,ii)
          tmp(k:) = ' '
          tmp(k:) = 'idx'
          status = rtdb_put(rtdb,tmp,mt_int,nuc_num,nuc_cent)
          go to 120
        else
          print *, 'no idx/symbol defined. By default: H'
          quantum_nuc_smb = 'H'
          go to 110
        end if
      end if

 120  continue

      end


      subroutine cneo_geom_init(rtdb,geom_cneo)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "nwc_const.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "geomP.fh"
#include "util.fh"
#include "stdio.fh"
#include "inp.fh"
#include "ccneo.fh"

      integer rtdb   ![input]
      integer geom_cneo   ![output]
      logical status
      character*255 tmp
      integer nat
      double precision charge_cneo, tmp_charge
      double precision masselec, avogadro, mass_au
      integer i, iatom
      parameter(masselec=9.1093826d-31)
      parameter(avogadro=6.022140857d23)
      mass_au = 1.0d-3/avogadro

      status = geom_create(geom_cneo, 'geometry_cneo')
      status = geom_rtdb_load(rtdb, geom_cneo, 'geometry')
      write(6,*) 'cneo_geom_init: new geom handle:', geom_cneo

c  read quantum_nuc
      status = geom_ncent(geom_cneo,nat)
c  update mass and charge

      tmp_charge = 0.0d0

      do i=1,nuc_num
         iatom = nuc_cent(i)
         if ((iatom.le.0) .or. (iatom .gt. nat))
     $         call errquit('cneo_geom_init:invalid nuc index',0,0)
         tmp_charge = tmp_charge-charge(iatom,geom_cneo)
         geom_mass(iatom,geom_cneo) = geom_mass(iatom,geom_cneo)
     $             -charge(iatom,geom_cneo)*masselec/mass_au
         charge(iatom,geom_cneo) = 0.0d0
      enddo

      if (.not. rtdb_get(rtdb,'charge',MT_DBL,1,charge_cneo))
     $      charge_cneo = 0.0d0

      if (.not. rtdb_put(rtdb,'charge_real',mt_dbl,1,charge_cneo))
     $      call errquit('cneo_geom_init:cannot save real charge',0,0)

      charge_cneo = charge_cneo+tmp_charge

      if (.not. rtdb_put(rtdb,'charge',MT_DBL,1,charge_cneo))
     $     call errquit('cneo_geom_init:charge update failed',0,0)

	   call geom_compute_values(geom_cneo)
      
      end

      subroutine cneo_bas_init(rtdb,geom_cneo,bases_nuc)
      implicit none

#include "errquit.fh"
#include "inp.fh"
#include "bas.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "ccneo.fh"

      integer rtdb, geom_cneo
      integer bases_nuc(max_quantum_nuc)
      logical status
      integer i, iatom

      bases_nuc = -999
      do i=1,nuc_num
         iatom = nuc_cent(i)
         write(6,*), 'cneo_bas_init: building basis for ',iatom
         if (.not. bas_create(bases_nuc(i), 'nuc basis'))
     $      call errquit('cneo init: bas create',0,0)
         call cneo_bas_load(rtdb,geom_cneo,bases_nuc(i),'nuc basis',
     $                                       iatom)
      enddo

      end
  
      subroutine cneo_bas_load(rtdb,geom_cneo,basisin,name,iatom)
      implicit none
      
#include "errquit.fh"

#include "rtdb.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "geom.fh"
#include "context.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "inp.fh"
#include "global.fh"
#include "stdio.fh"
#include "bas.fh"

#include "bas_exndcf.fh"


      integer rtdb,geom_cneo,basisin,iatom
      integer bas_cneo, lentmp, nexcf
      character*(*) name
      character*256 tmp
      logical rtdb_status
      integer h_tmp, k_tmp, h_new, k_new
      integer iunique, uniquecent, istar , ntag_read
      integer istart, iexcpt, ilen, junique
      integer i_array,l_array
      logical oIsexcept
      character*16 tag_to_add, tag_in_lib
      integer star_nr_tags
      integer ecXtra      ! amount of extra space required for zero shell info
      parameter (ecXtra = 2)
      character*25 nameex


      bas_cneo = basisin + Basis_Handle_Offset

      ibs_geom(bas_cneo) = geom_cneo
      bs_name(bas_cneo) = name
      len_bs_name(bas_cneo) = inp_strlen(name)

      if (.not.context_rtdb_match(rtdb,name,bs_trans(bas_cneo)))
     &       bs_trans(bas_cneo) = name
      len_bs_trans(bas_cneo) = inp_strlen(bs_trans(bas_cneo))

      rtdb_status = .true.

      tmp = 'basis:'//bs_trans(bas_cneo)(1:len_bs_trans(bas_cneo))
      lentmp = inp_strlen(tmp) + 1
c
      ntag_read=0
      tmp(lentmp:) = ' '
      tmp(lentmp:) = ':bs_nr_tags'
      rtdb_status = rtdb_status .and. 
     &              rtdb_get(rtdb,tmp,mt_int,1,ntag_read)
      if (ntag_read .gt. 0) then
          tmp(lentmp:) = ' ' 
          tmp(lentmp:) = ':bs_tags'
          rtdb_status = rtdb_status .and.
     &       rtdb_cget(rtdb, tmp, ntags_bsmx, bs_tags(1,bas_cneo))
c
          tmp(lentmp:) = ' ' 
          tmp(lentmp:) = ':bs_stdname'
          rtdb_status = rtdb_status .and.
     &       rtdb_cget(rtdb, tmp, ntags_bsmx, bs_stdname(1,bas_cneo))
      endif
c
      tmp(lentmp:) = ' ' 
      tmp(lentmp:) = ':assoc ecp name'
      rtdb_status = rtdb_status .and.
     &       rtdb_cget(rtdb, tmp, 2, name_assoc(1,bas_cneo))
c
      nexcf = 0
      tmp(lentmp:) = ' ' 
      tmp(lentmp:) = ':number of exps and coeffs'
      rtdb_status = rtdb_status .and.
     &    rtdb_get(rtdb,tmp,mt_int,1,nexcf)
c
      write(nameex,'(a23,i2)')' bas_cneo exps and coeffs ',bas_cneo
c
      if (exndcf(H_exndcf,bas_cneo) .ne. -1) then
          if (.not. ma_free_heap(exndcf(H_exndcf,bas_cneo)))
     $        call errquit('bas_rtdb_do_load: ma is corrupted',
     $        exndcf(H_exndcf,bas_cneo), RTDB_ERR)
      endif

      if (.not.ma_alloc_get(mt_dbl,(nexcf+ecXtra),nameex,
     &    h_tmp, k_tmp)) then
        write(LuOut,*)' not enough memory'
        call errquit
     &      (' bas_rtdb_do_load: error allocating space'//
     &      ' for exndcf',911, RTDB_ERR)
      else
        call dfill((nexcf+ecXtra),0.0d00,dbl_mb(k_tmp),1)
        exndcf(H_exndcf,bas_cneo) = h_tmp
        exndcf(K_exndcf,bas_cneo) = k_tmp
        exndcf(SZ_exndcf,bas_cneo) = (nexcf+ecXtra)
      endif
c
      if (nexcf .gt. 0) then
          tmp(lentmp:) = ' ' 
          tmp(lentmp:) = ':exps and coeffs'
          rtdb_status = rtdb_status .and.
     &       rtdb_get(
     &       rtdb, tmp, mt_dbl, nexcf, dbl_mb(mb_exndcf(1,bas_cneo)))
      endif
c
      tmp(lentmp:) = ' ' 
      tmp(lentmp:) = ':header'
      rtdb_status = rtdb_status .and.
     &       rtdb_get(
     &       rtdb, tmp, mt_int, ndbs_head, infbs_head(1,bas_cneo))
c
      tmp(lentmp:) = ' ' 
      tmp(lentmp:) = ':tags info'
      rtdb_status = rtdb_status .and.
     &       rtdb_get(
     &       rtdb, tmp, mt_int,
     &       ndbs_tags*ntags_bsmx, infbs_tags(1,1,bas_cneo))
c
      tmp(lentmp:) = ' ' 
      tmp(lentmp:) = ':contraction info'
      rtdb_status = rtdb_status .and.
     &       rtdb_get(
     &       rtdb, tmp, mt_int,
     &       ndbs_ucont*nucont_bsmx, infbs_cont(1,1,bas_cneo))
c
c now deal with the star tags from the input. We have the geometry, so
c we extract the tags from the geometry, load the star tag info, and
c add the bas_cneo sets
c
c first check if we have any star tags to deal with at all
c account for old rtdbs without star tag info
c
      tmp(lentmp:) = ' '
      tmp(lentmp:) = ':star nr tags'
      star_nr_tags = 0
      if (rtdb_status) rtdb_status = rtdb_status .and. rtdb_get(rtdb,
     &       tmp,mt_int,1,star_nr_tags)
c
      if (star_nr_tags .ne. 0)
     $   call errquit('cneo:star input currently not supported',0,0)
c
      star_nr_tags = 0
c
c read the basis now get check status of read operations
c
      if (.not.rtdb_status) then
          if (exndcf(H_exndcf,bas_cneo) .ne. -1) then
            if (.not. ma_free_heap(exndcf(H_exndcf,bas_cneo)))
     $           call errquit('bas_rtdb_load: ma free failed?',0,
     &       MA_ERR)
            exndcf(H_exndcf,bas_cneo) = -1
            exndcf(K_exndcf,bas_cneo) = 0
            exndcf(SZ_exndcf,bas_cneo)= 0
          endif

          return
      endif

      call cneo_bas_geobas_build(basisin, iatom)

      
      end

      subroutine cneo_bas_geobas_build(basisin, iatom)
      implicit none
#include "errquit.fh"
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "geom.fh"
#include "context.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "inp.fh"
#include "global.fh"
#include "bas_ibs_dec.fh"
#include "stdio.fh"
c::passed
      integer basisin, iatom
c::local
      integer basis
      integer geom
      integer nat
      integer i, idum_cont, idum_at
      integer j, jstart, jend, jsize
      integer kstart, kend, ksize, lsize, icount
      integer nbf, iu_cont, myang
      integer my_gen, my_type
      integer atn
      character*16 element
      character*2  symbol
      logical status
      logical foundit
      logical found_any
      logical is_bq
      logical ecpORso, relbas
      integer int_dummy, num_elec
      integer h_tmp, k_tmp
      character*2 tag12
      character*16 name_tmp
      double precision erep_save
      integer uce, idbstag
*debug:mem      integer inode
c::functions
      logical bas_high_angular
      external bas_high_angular
      integer nbf_from_ucont
      external nbf_from_ucont
      logical ecp_get_num_elec
      external ecp_get_num_elec
      logical bas_match_tags
      external bas_match_tags
      logical basis_is_rel
      external basis_is_rel

      integer nuc_num
      character*256 tmp
      integer quantum_nuc(2)
c
#include "bas_ibs_sfn.fh"
c
      basis = basisin + Basis_Handle_Offset
      geom  = ibs_geom(basis)
      ecpORso = .false.
      relbas = basis_is_rel(basisin)
c
*debug:mem      do inode = 0,(ga_nnodes() - 1)
*debug:mem        if (inode.eq.ga_nodeid()) call MA_summarize_allocated_blocks()
*debug:mem        call ga_sync()
*debug:mem      enddo
c

      status = geom_ncent(geom, nat)
      
      if (nat.eq.0.or..not.status) then
        write(LuOut,*)' bas_geobas_build: ERROR '
        write(LuOut,*)' number of centers is zero or weird'
        write(LuOut,*)' nat = ',nat
c..... add diagnostics later
        return
      endif
c.... set spherical flag
      if (infbs_head(HEAD_SPH,basis).eq.1) then
        bas_spherical(basis) = .true.
      else
        bas_spherical(basis) = .false.
      endif
c.... set flag if any general contractions are present
      bas_any_gc(basis)       = .false.
      if(.not. ecpORso) then
        do i = 1,(infbs_head(Head_Ncont,basis))
          if (.not.bas_any_gc(basis)) then
            my_gen  = infbs_cont(Cont_Ngen,i,basis)
            my_type = infbs_cont(Cont_Type,i,basis)
            if (my_gen.gt.1.and.my_type.ge.0)
     &          bas_any_gc(basis) = .true.
          endif
        enddo
      endif
c.... set flag if any sp (spd,spdf) shells are present
      bas_any_sp_shell(basis) = .false.
      if(.not. ecpORso) then
        do i = 1,(infbs_head(Head_Ncont,basis))
          if (.not.bas_any_sp_shell(basis)) then
            my_gen  = infbs_cont(Cont_Ngen,i,basis)
            my_type = infbs_cont(Cont_Type,i,basis)
            if (my_type.lt.0) then
              if (my_gen.ne.2) then
                write(luout,*)' sp shell with n_gen = ',my_gen
                call errquit ('bas_geobas_build: fatal error',911,
     &       BASIS_ERR)
              endif
              bas_any_sp_shell(basis) = .true.
            endif
          endif
        enddo
      endif
c
c... clear old ibs_ce2uce if it exists
      if (ibs_ce2uce(SZ_ibs,basis).gt.0) then
*debug:mem        write(LuOut,*)' clearing old ce2uce data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        h_tmp = ibs_ce2uce(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('bas_geobas_build: error freeing ibs_ce2uce',911,
     &       BASIS_ERR)
        ibs_ce2uce(H_ibs,basis)  = 0
        ibs_ce2uce(K_ibs,basis)  = 0
        ibs_ce2uce(SZ_ibs,basis) = 0 
      endif
c                                 123456789012
      write(name_tmp,'(a12,i4)') ' ibs_ce2uce ',basis
      if (.not.ma_alloc_get(mt_int,nat,name_tmp,h_tmp,k_tmp)) then
        call errquit
     &      ('bas_geobas_build: error ma_alloc ibs_ce2uce',911,
     &       MA_ERR)
      else
*debug:mem        write(LuOut,*)' generating ce2uce data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        call ifill(nat,0,int_mb(k_tmp),1)
        ibs_ce2uce(H_ibs,basis)  = h_tmp
        ibs_ce2uce(K_ibs,basis)  = k_tmp
        ibs_ce2uce(SZ_ibs,basis) = nat
*debug:mem        do inode = 0,(ga_nnodes() - 1)
*debug:mem          if (inode.eq.ga_nodeid()) call MA_summarize_allocated_blocks()
*debug:mem          call ga_sync()
*debug:mem          call util_flush(LuOut)
*debug:mem        enddo
      endif
c
c build center to unique center map
c
      found_any = .false.
      do 00100 i=1,nat
        foundit = .false.
*before match_tags:        do 00101 j = 1,infbs_head(HEAD_NTAGS,basis)
*before match_tags:          if(bas_match_tags(tags(i,geom),bs_tags(j,basis))) then
*before match_tags:            int_mb(mb_ibs_ce2uce(i,basis)) = j
*before match_tags:            foundit = .true.
*before match_tags:            goto 00102
*before match_tags:          endif
*before match_tags:00101   continue
          if (i .eq. iatom) then

            if (bas_match_tags(tags(i,geom),basisin,j)) then
                int_mb(mb_ibs_ce2uce(i,basis)) = j
                foundit = .true.
                found_any = .true.
                goto 00102
            endif
            if (.not. foundit .and. .not. (ecpORso .or. relbas)) then
                if (geom_tag_to_element(tags(i,geom), symbol, element,
     $            atn)) then
                  if (ga_nodeid().eq.0)
     &              write(LuOut,10) 
     &              tags(i,geom)(1:inp_strlen(tags(i,geom))),
     &              element(1:inp_strlen(element)),
     &              bs_name(basis)(1:len_bs_name(basis))
 10              format(/' ERROR: geometry tag ',a,' (',a,
     &             ') is a quantum nuclei ',
     $              'but has no functions in basis "',a,'"'/
     $              ' ERROR: only bq* centers can have no functions')
                  if (ga_nodeid().eq.0) call util_flush(LuOut)
                  call errquit
     &              ('bas_geobas_build: basis/geometry mismatch', 0,
     &           BASIS_ERR)
                else
                  tag12 = tags(i,geom)(1:2)
                  is_bq = inp_compare(.false.,'bq',tag12) .or.
     $               (inp_compare(.false.,'X',tags(i,geom)(1:1)) .and. 
     $               (.not. inp_compare(.false.,'e',tags(i,geom)(2:2))))
                  if (ga_nodeid().eq.0 .and.(.not.is_bq))
     &              write(LuOut,11) i,
     &              tags(i,geom)(1:inp_strlen(tags(i,geom))),
     $              bs_name(basis)(1:len_bs_name(basis))
 11              format(/'WARNING: geometry tag ',i4, ' ', a,
     $              ' not found in basis "',a,'"'/)
                  int_mb(mb_ibs_ce2uce(i,basis)) = 0
                  if (ga_nodeid().eq.0) call util_flush(LuOut)
                endif
            endif
00102       continue
          else
            int_mb(mb_ibs_ce2uce(i,basis)) = 0
          endif
00100 continue
*
      if (.not.found_any) then
        if (ecpORso) then
          return
        else
          if (ga_nodeid().eq.0) then
            write(luout,*)' none of the geometry tags matched any ',
     &          'basis set tag in the basis "',
     &          bs_name(basis)(1:len_bs_name(basis)),
     &          '"'
          endif
          call errquit('bas_geobas_build: fatal error',911,
     &       BASIS_ERR)
        endif
      endif

c
c build total # of contractions 
c      
*debug:mem        do inode = 0,(ga_nnodes() - 1)
*debug:mem          if (inode.eq.ga_nodeid()) then
*debug:mem            call MA_summarize_allocated_blocks()
*debug:mem            if (MA_verify_allocator_stuff()) then
*debug:mem              write(LuOut,*)' no errors'
*debug:mem            else
*debug:mem              write(LuOut,*)' errors'
*debug:mem            endif
*debug:mem          endif
*debug:mem          call ga_sync()
*debug:mem          call util_flush(LuOut)
*debug:mem        enddo
*debug:mem      call bas_print_allocated_info('bas_geobas_build 1')
      ncont_tot_gb(basis)  = 0
      do 10200 i=1,nat
          uce = sf_ibs_ce2uce(i,basis)
*debug:mem         write(LuOut,*)' myuce = ',uce,i,ga_nodeid()
*debug:mem         call util_flush(LuOut)
          if (uce.gt.0) then
            idum_cont = infbs_tags(TAG_NCONT,uce,basis)
            ncont_tot_gb(basis)   = idum_cont + ncont_tot_gb(basis)
          endif
10200 continue
c
c allocate space for center -> contraction range map
c
c... clear old ibs_ce2cnr if it exists
      if (ibs_ce2cnr(SZ_ibs,basis).gt.0) then
*debug:mem        write(LuOut,*)' clearing old ce2cnr data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        h_tmp = ibs_ce2cnr(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('bas_geobas_build: error freeing ibs_ce2cnr',911,
     &       MEM_ERR)
        ibs_ce2cnr(H_ibs,basis)  = 0
        ibs_ce2cnr(K_ibs,basis)  = 0
        ibs_ce2cnr(SZ_ibs,basis) = 0 
      endif
c                                 123456789012
      write(name_tmp,'(a12,i4)') ' ibs_ce2cnr ',basis
      if (.not.ma_alloc_get(mt_int,(2*nat),
     &      name_tmp,h_tmp,k_tmp)) then
        call errquit
     &      ('bas_geobas_build: error ma_alloc ibs_ce2cnr',911, MA_ERR)
      else
*debug:mem        write(LuOut,*)' generating ce2cnr data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        call ifill((2*nat),0,int_mb(k_tmp),1)
        ibs_ce2cnr(H_ibs,basis)  = h_tmp
        ibs_ce2cnr(K_ibs,basis)  = k_tmp
        ibs_ce2cnr(SZ_ibs,basis) = 2*nat
*debug:mem        do inode = 0,(ga_nnodes() - 1)
*debug:mem          if (inode.eq.ga_nodeid()) call MA_summarize_allocated_blocks()
*debug:mem          call ga_sync()
*debug:mem          call util_flush(LuOut)
*debug:mem        enddo
      endif
c
c build center -> contraction range map 
c
      int_dummy = 0
      do 00200 i=1,nat
          if (sf_ibs_ce2uce(i,basis).gt.0) then
            idum_cont =
     &           infbs_tags(TAG_NCONT,sf_ibs_ce2uce(i,basis),basis)
            int_mb(mb_ibs_ce2cnr(1,i,basis)) = int_dummy + 1
            int_mb(mb_ibs_ce2cnr(2,i,basis)) = int_dummy + idum_cont
            int_dummy = idum_cont + int_dummy
          else
*. . . . . . . . . . . . . . . . . . . . ! No functions on this center
            int_mb(mb_ibs_ce2cnr(1,i,basis)) = 0 
            int_mb(mb_ibs_ce2cnr(2,i,basis)) = -1
          endif
00200 continue
c
      if (ncont_tot_gb(basis) .eq. 0) call errquit
     $    ('bas_geobas_build: no functions in basis set', 0,
     &       BASIS_ERR)
      if (ncont_tot_gb(basis) .gt. ncont_mx) then
        write(LuOut,*)' number of contractions     = ',
     &      ncont_tot_gb(basis)
        write(LuOut,*)' number of contractions max = ',ncont_mx
        call errquit ('bas_geobas_build: ncont.gt.ncont_mx ',911,
     &       BASIS_ERR)
      endif
c
c allocate space for contraction -> center map
c
c... clear old ibs_cn2ce if it exists
      if (ibs_cn2ce(SZ_ibs,basis).gt.0) then
*debug:mem        write(LuOut,*)' clearing old cn2ce data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        h_tmp = ibs_cn2ce(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('bas_geobas_build: error freeing ibs_cn2ce',911,
     &       MEM_ERR)
        ibs_cn2ce(H_ibs,basis)  = 0
        ibs_cn2ce(K_ibs,basis)  = 0
        ibs_cn2ce(SZ_ibs,basis) = 0 
      endif
c                                 123456789012
      write(name_tmp,'(a12,i4)') ' ibs_cn2ce  ',basis
      if (.not.ma_alloc_get(mt_int,(1+ncont_tot_gb(basis)),
     &      name_tmp,h_tmp,k_tmp)) then
        call errquit('bas_geobas_build: error ma_alloc ibs_cn2ce',911,
     &       MA_ERR)
      else
*debug:mem        write(LuOut,*)' generating cn2ce data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        call ifill((1+ncont_tot_gb(basis)),0,int_mb(k_tmp),1)
        ibs_cn2ce(H_ibs,basis)  = h_tmp
        ibs_cn2ce(K_ibs,basis)  = k_tmp
        ibs_cn2ce(SZ_ibs,basis) = 1+ncont_tot_gb(basis)
      endif
c      
c build contraction -> center map
c
      do 00300 i=1,nat
          if (sf_ibs_ce2uce(i,basis).gt.0) then
            jstart = sf_ibs_ce2cnr(1,i,basis)
            jend   = sf_ibs_ce2cnr(2,i,basis)
            do 00400 j=jstart,jend
                int_mb(mb_ibs_cn2ce(j,basis)) = i
00400       continue
          endif
00300 continue
c
c set zero element of cn2ce to something useless
c
      int_mb(mb_ibs_cn2ce(0,basis)) = -1
c
c allocate space for ibs_cn2ucn map
c
c... clear old ibs_cn2ucn if it exists
      if (ibs_cn2ucn(SZ_ibs,basis).gt.0) then
*debug:mem        write(LuOut,*)' clearing old cn2ucn data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        h_tmp = ibs_cn2ucn(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('bas_geobas_build: error freeing ibs_cn2ucn',911,
     &       BASIS_ERR)
        ibs_cn2ucn(H_ibs,basis)  = 0
        ibs_cn2ucn(K_ibs,basis)  = 0
        ibs_cn2ucn(SZ_ibs,basis) = 0 
      endif
c                                 123456789012
      write(name_tmp,'(a12,i4)') ' ibs_cn2ucn ',basis
      if (.not.ma_alloc_get(mt_int,(1+ncont_tot_gb(basis)),
     &      name_tmp,h_tmp,k_tmp)) then
        call errquit
     &      ('bas_geobas_build: error ma_alloc ibs_cn2ucn',911,
     &       MA_ERR)
      else
*debug:mem        write(LuOut,*)' generating cn2ucn data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        call ifill((1+ncont_tot_gb(basis)),0,int_mb(k_tmp),1)
        ibs_cn2ucn(H_ibs,basis)  = h_tmp
        ibs_cn2ucn(K_ibs,basis)  = k_tmp
        ibs_cn2ucn(SZ_ibs,basis) = 1+ncont_tot_gb(basis)
      endif
c
c build contraction -> unique contraction map
c
      do 00500 i=1,nat
          jstart = sf_ibs_ce2cnr(1,i,basis)
          jend   = sf_ibs_ce2cnr(2,i,basis)
          jsize  = jend - jstart + 1
          if (jsize .gt. 0) then
            idum_at = sf_ibs_ce2uce(i,basis)
            kstart = infbs_tags(TAG_FCONT,idum_at,basis)
            kend   = infbs_tags(TAG_LCONT,idum_at,basis)
            ksize  = kend - kstart + 1
            lsize  = infbs_tags(TAG_NCONT,idum_at,basis)
            if (jsize.eq.ksize.and.ksize.eq.lsize) then
                icount = 0
                do 00600 j=jstart,jend
                  int_mb(mb_ibs_cn2ucn(j,basis)) = kstart + icount
                  icount = icount + 1
00600          continue
            else
               write(LuOut,*)' bas_geobas_build: ERROR '
               write(LuOut,*)' contraction range size mismatch'
               write(LuOut,*)'        cont. range (',jstart,':',jend,')'
               write(LuOut,*)' unique cont. range (',kstart,':',kend,')'
               write(LuOut,*)'        cont. size: ',jsize
               write(LuOut,*)' calculated unique cont. size: ',ksize
               write(LuOut,*)'     lookup unique cont. size: ',lsize
               return
            endif
          endif
00500 continue
c
c set zero element
c
      int_mb(mb_ibs_cn2ucn(0,basis)) = 0
c
c allocate space for ibs_cn2bfr map
c
c... clear old ibs_cn2bfr if it exists
      if (ibs_cn2bfr(SZ_ibs,basis).gt.0) then
*debug:mem        write(LuOut,*)' clearing old cn2bfr data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        h_tmp = ibs_cn2bfr(H_ibs,basis)
        if (.not.ma_free_heap(h_tmp)) call errquit
     &        ('bas_geobas_build: error freeing ibs_cn2bfr',911,
     &       BASIS_ERR)
        ibs_cn2bfr(H_ibs,basis)  = 0
        ibs_cn2bfr(K_ibs,basis)  = 0
        ibs_cn2bfr(SZ_ibs,basis) = 0 
      endif
c                                 123456789012
      write(name_tmp,'(a12,i4)') ' ibs_cn2bfr ',basis
      if (.not.ma_alloc_get(mt_int,(2*(1+ncont_tot_gb(basis))),
     &      name_tmp,h_tmp,k_tmp)) then
        call errquit
     &      ('bas_geobas_build: error ma_alloc ibs_cn2bfr',911,
     &       BASIS_ERR)
      else
*debug:mem        write(LuOut,*)' generating cn2bfr data ',ga_nodeid()
*debug:mem        call util_flush(LuOut)
        call ifill((2*(1+ncont_tot_gb(basis))),0,int_mb(k_tmp),1)
        ibs_cn2bfr(H_ibs,basis)  = h_tmp
        ibs_cn2bfr(K_ibs,basis)  = k_tmp
        ibs_cn2bfr(SZ_ibs,basis) = 2*(1+ncont_tot_gb(basis))
      endif
c
c build nprim_tot_gb, ncoef_tot_gb, nbf_tot_gb, and 
c contraction -> basis function range map 
c find nbfmax for basis (initialized in block data statement)
c
      nbf_tot_gb(basis)   = 0
      nprim_tot_gb(basis) = 0
      do 00700 i = 1,ncont_tot_gb(basis)
        iu_cont = sf_ibs_cn2ucn(i,basis)
c
        nbf = nbf_from_ucont(iu_cont,basisin)
        nbfmax_bs(basis) = max(nbfmax_bs(basis),nbf)
c
        int_mb(mb_ibs_cn2bfr(1,i,basis)) = nbf_tot_gb(basis) + 1
        int_mb(mb_ibs_cn2bfr(2,i,basis)) = nbf_tot_gb(basis) + nbf

        nbf_tot_gb(basis) = nbf_tot_gb(basis) + nbf
        nprim_tot_gb(basis) = nprim_tot_gb(basis) +
     &         infbs_cont(CONT_NPRIM,iu_cont,basis)
        ncoef_tot_gb(basis) = ncoef_tot_gb(basis) + 
     &         infbs_cont(CONT_NPRIM,iu_cont,basis)*
     &         infbs_cont(CONT_NGEN,iu_cont,basis)
00700 continue
c
c set zero elements of cn2bfr
c
      int_mb(mb_ibs_cn2bfr(1,0,basis)) = 0
      int_mb(mb_ibs_cn2bfr(2,0,basis)) = 0

c
c build high angular momentum of this loaded <basis|geom> pair
c note angular_bs(*) initialized in block data function
c
      if (.not.bas_high_angular(basisin,myang))call errquit
     &      ('bas_geobas_build: error bas_high_angular',911,
     &       BASIS_ERR)
*

      end


      subroutine cneo_int1er_init(rtdb, geom, bases)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "rtdb.fh"
#include "global.fh"
#include "bas.fh"
#include "ccneo.fh"
#include "cconstrain.fh"
#include "util.fh"

      integer rtdb
      integer geom
      integer bases(nuc_num), basis
      integer iatom,i,j
      double precision center(3)
      integer g_dens, nbf
      double precision moments(4)
      integer g_int(3)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      character*5 str_i

      do i=1,nuc_num
         iatom = nuc_cent(i)
         basis = bases(i)
         str_i = ' '
         write(str_i,'(I0)') i

         call dcopy(3,coords(1,iatom,geom),1,exp_pos(1,i),1)
         call dcopy(3, exp_pos(1,i),1,center(1),1)

         g_int1e_r(1,i)=ga_create_atom_blocked(geom,basis,
     $      'int1er:x'//trim(str_i))
         g_int1e_r(2,i)=ga_create_atom_blocked(geom,basis,
     $      'int1er:y'//trim(str_i))
         g_int1e_r(3,i)=ga_create_atom_blocked(geom,basis,
     $      'int1er:z'//trim(str_i))

         call int1e_r_ao(basis, center, g_int1e_r(1,i))

      enddo
      if (util_print('cneo_int1er', print_debug)) then
         do i=1,nuc_num
            do j=1,3
               call ga_print(g_int1e_r(j,i))
            enddo
         enddo
      endif

      end

      subroutine int1e_r_ao(basis,center,g_int)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

      integer basis
      double precision center(3)
      integer g_int(3)

      integer g_s
      integer i

      if (.not. ga_duplicate(g_int(1),g_s,'int1e:tmp_ovlp'))
     $      call errquit('int1e_r_ao:create ovlp')
      call ga_zero(g_s)
      call int_1e_ga(basis,basis,g_s,'overlap',.false.)

      call int_dip_ga(basis,basis,g_int(1),g_int(2),g_int(3))
      do i=1,3
         call ga_dadd(1.0d0,g_int(i),-center(i),g_s,g_int(i))
      enddo

      if (.not. ga_destroy(g_s)) call errquit
     $      ('int1e_r_ao:ga_destroy?',0,0)

      end

      subroutine cneo_vecs_guess(rtdb)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "ccneo.fh"
#include "cscf_cneo.fh"

      integer rtdb
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical oprint
      integer g_dens(2), g_odens
      integer nproc
      integer nxtask
      external nxtask
      integer i, iatom

      nproc = ga_nnodes()

      oprint = util_print('mo guess',print_default) .and.
     $      ga_nodeid().eq.0
      
      call scf_lindep(rtdb, geom_cneo, bas_e, g_movecs_e, nmo_e)
      call cneo_schwarz_set(bas_e)

      g_dens(1) = ga_create_atom_blocked(geom_cneo, bas_e, 
     $      'Guess density:e')
      call guess_dens(rtdb, geom_cneo, bas_e, g_dens(1))

      call fock_force_direct(rtdb)
      call rhf_dens_to_mo(rtdb, geom_cneo, bas_e, nelec_e,
     &        nbf_e, nmo_e, 1d-7,
     &        g_dens(1), g_movecs_e, .false., dbl_mb(k_eval_e))

      call movecs_swap(rtdb, 
     &     'scf', 'RHF', 
     &     g_movecs_e, dbl_mb(k_occ_e), dbl_mb(k_eval_e))


      call movecs_fix_phase(g_movecs_e)
      call rohf_dens(g_movecs_e,nbf_e,nclosed_e,nopen_e, 
     $                                    g_dens(1),g_dens(2))

      do i=1,nuc_num
         iatom = nuc_cent(i)
         call cneo_vecs_guess_nuc(rtdb, i,iatom,g_dens(1))
      enddo
      if (.not. ga_destroy(g_dens(1))) call errquit
     &        ('scf_vectors_guess: destroy of density failed',0,
     &       GA_ERR)

      end

      subroutine cneo_vecs_guess_nuc(rtdb, i,iatom,g_edens)
      implicit none 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "ccneo.fh"
#include "cscf_cneo.fh"

      integer rtdb,i, iatom
      integer g_edens

      integer g_nhcore, g_nfock, g_over, g_tmp
      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      g_nfock = ga_create_atom_blocked(geom_cneo, bas_n(i), 
     $            'Temp:n:fock')
      if (.not. ga_create(MT_DBL,nbf_n,nbf_n,'cneo_vecs_guess_n:over',
     $     0, 0, g_over)) call
     $     errquit('cneo_vecs_guess_n: ga failed for over', 0, GA_ERR)

      if (.not. ga_create(MT_DBL,nbf_n,nbf_n,'cneo_vecs_guess_n:tmp',
     $     0, 0, g_tmp)) call
     $     errquit('cneo_vecs_guess_n: ga failed for tmp', 0, GA_ERR)
      if (.not. ga_create(MT_DBL,nbf_n,nbf_n,'cneo_vecs_guess_n:over',
     $     0, 0, g_nhcore)) call
     $     errquit('cneo_vecs_guess_n: ga failed for over', 0, GA_ERR)

      call ga_zero(g_nfock)
      call get_j_dm(geom_cneo, bas_n(i), bas_e, -1.0d0, 1.0d-9,
     $       .false., g_edens, g_nfock, 'neints', i,0,.false.)
      call hcore_nuc(geom_cneo, bas_n(i), g_nhcore, iatom)
      call ga_dadd(1.0d0, g_nhcore, 1.0d0, g_nfock, g_nfock)
      if (util_print('cneo_vec_nuc', print_debug)) 
     $         call ga_print(g_nfock)

      call ga_zero(g_over)
      call int_1e_ga(bas_n(i), bas_n(i), g_over, 'overlap', .false.)

      call ga_diag(g_nfock, g_over, g_tmp, dbl_mb(k_eval_n(i)))
      call ga_copy(g_tmp, g_movecs_n(i))

      if (.not. ga_destroy(g_tmp)) call errquit
     $      ('cneo_vecs_guess_nuc: destroy ga',0,GA_ERR)
      if (.not. ga_destroy(g_over)) call errquit
     $      ('cneo_vecs_guess_nuc: destroy ga',0,GA_ERR)
      if (.not. ga_destroy(g_nfock)) call errquit
     $      ('cneo_vecs_guess_nuc: destroy ga',0,GA_ERR)
      if (.not. ga_destroy(g_nhcore)) call errquit
     $      ('cneo_vecs_guess_nuc: destroy ga',0,GA_ERR)

      if (util_print('cneo_vec_guess_nuc', print_debug)) 
     $         call ga_print(g_movecs_n(i))

      end