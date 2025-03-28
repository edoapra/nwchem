c
c     Functions are modified from basis.F to process nuclear basis
c
      subroutine cneo_bas_input(rtdb)
c $Id$
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "inp.fh"
#include "stdio.fh"
c
c   basis "nuc basis" [spherical|cartesian] [print|noprint]
c
c     tag library [tag in library] <standard set> [file <filename>]
c     
c     end basis
c
c
c::functions
      logical  bas_do_destroy
      external bas_do_destroy
c::passed
      integer rtdb              ! [input] handle to database
c::local
      integer nopt
      parameter (nopt = 4)
      character*10 opts(nopt)
      character*255 test, name, filename, standard
      character*255 ecpname, soname
      character*256 mydir
      character*1000 errmsg
      logical status, ospherical, osegment, oprint
      logical o_add_ecpname, o_add_soname
      logical oHas_Star
      logical oIs_ecp, oIs_so, oIs_rel
      integer ind, basis
      logical bas_add_ucnt_init, bas_set_spherical, bas_add_ucnt_tidy
      external bas_add_ucnt_init, bas_set_spherical
      data opts /
     &    'spherical', 'cartesian', 
     &    'print',     'noprint'/
c
      ecpname = '                                                '
      soname = '                                                '
      mydir = ' '
c
c     Check is a basis/ecp directive and read in name of the basis
c
      call inp_set_field(0)
      status = inp_a(test)
      oIs_ecp = .false.
      oIs_so = .false.
      status = status.and.
     &    (inp_compare(.false.,test,'basis') .or.
     &     oIs_ecp .or. oIs_so)
      if (.not.status) goto 10000
*:debug-s
*debug:      write(luout,*)' debug output - start'
*debug:      write(luout,*)'      test?        ',test
*debug:      write(luout,*)' oIs_ecp is        :',oIs_ecp
*debug:      write(luout,*)' debug output - end  '
*:debug-e
c
c     Parse rest of basis/ecp directive line
c
      name = ' '
      filename = ' '
      test = ' '
      standard = ' '
      ospherical = .false. ! Default is cartesian
      osegment   = .true.  ! Default is to force segmentation
      oprint     = .true.  ! Default is to print the basis on input
      o_add_ecpname = .false. ! Default is to not associate an ecp name
      o_add_soname = .false.  ! Default is to not associate a so name
      oIs_rel = .false.    ! Default is nonrelativistic
      oHas_Star = .false.  ! Default, not using star-tags 
c
 10   if (inp_a(test)) then
c
         if (.not. inp_match(nopt, .false., test, opts, ind)) then
c
c     Not a recognized option ... the name of the basis or an error
c
            if ((name.ne.' ') .or. (inp_cur_field() .ne. 2)) then
               write(LuOut,*)
     &            ' bas_input: basis/ecp/so name must be first option'
               goto 10000
            endif
            name = test
            goto 10
         endif
c
         goto (100, 200, 700, 800) ind
         goto 10000
c
c     spherical
c
 100     ospherical = .true.
         goto 10
c
c     cartesian
c     
 200     ospherical = .false.
         goto 10
c
c     print
c
 700     oprint = .true.
         goto 10
c
c     noprint
c
 800     oprint = .false.
         goto 10
c
      endif
c
c     Open a new basis set to receive the new data
c
      if (name .eq. ' ') name = 'nuc basis'
      if (.not. bas_create(basis, name))
     $     call errquit('bas_input: failed to create basis', 0,
     &       BASIS_ERR)
c
c    Now left with reading in from the input additional specifications
c    for basis functions or standard sets on specific tags
c
      if (.not.bas_add_ucnt_init(basis)) call errquit
     &    ('bas_input: failed to init add_ucnt ',911,
     &       BASIS_ERR)
      if (.not. bas_set_spherical(basis, ospherical))
     &      call errquit
     &      (' bas_set_spherical failed ',911, BASIS_ERR)
      
      call cneo_bas_input_body(basis, osegment, oIs_rel, oHas_Star)
c
c
c     Now have processed the entire basis directive.  Print out
c     info if desired, write it to the data base, tidy up and go home
c
c     Always print summary. Only print full basis when it does not
c     contain any star tags. Do store the info that the basis should be 
c     printed on rtdb, so taht it can be printed at the bas_rtdb_load 
c     stage
c
      if (oprint) then
         test = 'basisprint:'//name(1:inp_strlen(name))
         if (.not. rtdb_put(rtdb,test,mt_log,1,oHas_Star))
     $          call errquit('bas_input: rtdb_put * failed', 0,
     &       RTDB_ERR)
         if (.not. oHas_Star) then
            if (.not. bas_print(basis))
     $          call errquit('bas_input: print failed', 0, BASIS_ERR)
         endif
         if (.not.bas_summary_print(basis))
     &      call errquit('scf: basis summary print failed',911,
     &       BASIS_ERR)
      endif
c
      if (.not. bas_rtdb_store(rtdb, name, basis)) call errquit
     $     ('bas_input: failed to store basis', 0, BASIS_ERR)
c
      if (.not. bas_add_ucnt_tidy(basis))
     $     call errquit('bas_input: tidy failed',0, BASIS_ERR)
c
      if (.not. bas_do_destroy(basis)) call errquit
     $     ('bas_input: bas_destroy failed', 0, BASIS_ERR)
c
      return
c
10000 write(LuOut,10001)
10001 format(' basis [<name>] [library <standard set>] \\',/,
     $    '      [file <filename>] [segment]\\',
     $    '      [spherical|cartesian]\\',
     &    '      [print|noprint] [rel] [ecpset name_of_ecp]',
     &    ' [soset name_of_so_potential] [bse]')
      call errquit('bas_input: invalid format for basis directive', 0,
     &       INPUT_ERR)
c
      end
*.....................................................................
      subroutine cneo_bas_input_body(basis,osegment,oIs_rel,oHas_Star)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "inp.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "stdio.fh"
#include "util.fh"
c::passed
      integer basis             ! [input]
      logical osegment          ! [input]
      logical oHas_Star         ! [output] true if basis input contains star-tags
c
c     Read the body of a basis directive that describes the
c     tags/exponents/contraction coefficients 
c
c
c     tag library [<tag in library>] <standard set> [file <filename>]
c       ...
c     end basis
c     
      character*16 tag, tag_in_lib
      character*16 type_of_cont, rel_shell
      character*255 field, standard, filename
      integer nltypes           ! No. of known angular momentum types
      integer nsptypes          ! No. of known sp type shells
      integer nopts             ! No. of options
      parameter (nltypes = 8, nsptypes = 3, nopts = 2)
      character*1 ltypes(nltypes)
      character*2 sptypes(nsptypes)
      character*8 opts(nopts)
      logical oIs_ecp            
      logical oIs_so
      logical oIs_rel
      logical oshell_is_rel
      logical oIs_star
      integer spvalues(nsptypes)
      integer num_elec
      integer l_value, ind, icount
      character*(nw_max_path_len) compiled_name
cc AJL/Begin/SPIN-POLARISED ECPs
      integer channel           ! Both (Default)=0; Alpha=1; Beta=2
cc AJL/End
      data ltypes /'s','p','d','f','g','h','i','k'/
      data sptypes / 'sp', 'l ', 'ul'/
      data spvalues/  -1 , -1 ,   -1 /
      data opts / 'library', 'file' /
c
      oIs_ecp = .false.
      oIs_so  = .false.
c
c     Input a new line
c
 10   if (.not. inp_read()) call errquit
     $     ('bas_input_body: premature EOF', 0, INPUT_ERR)
c
c     Start parsing current line
c
 20   call inp_set_field(0)
      standard = ' '            ! Must reset these for every tag
      call util_nwchem_srcdir(compiled_name)
      compiled_name=compiled_name(1:inp_strlen(compiled_name))
     $     //"/cneo/basis/"
      filename = compiled_name
      oshell_is_rel = oIs_rel
      oIs_star = .false.
      if (.not. inp_a(tag)) goto 10000
c
c  Check if we have a tag that contains a '*'. If so, we have to store
c  the data "as is" in the RTDB until we know the full geometry without 
c  processing
c
      ind = inp_strlen(tag)
c
*. . . . . . . . . . . . . . . . . . . . . . ! End of basis directive
      if (inp_compare(.false., 'end', tag)) goto 9000 
c
      if (.not. inp_a(type_of_cont)) goto 10000

      call inp_prev_field
 22   if (inp_a(field)) then
         if (inp_match(nopts, .false., field, opts, ind)) then
            goto (25, 26) ind
c     
c     library
c     
 25         if (.not. inp_a(standard)) goto 10000
            tag_in_lib = tag
            if (inp_a(field)) then
               if (inp_match(nopts, .false., field, opts, ind)) then
                  call inp_prev_field
               else
                  tag_in_lib = standard
                  standard   = field
*                   write(LuOut,*) ' tag_in_lib, standard ', 
*     $                tag_in_lib, standard
               endif
            endif
            goto 22
c     
c     file
c     
 26         if (.not. inp_a(filename)) goto 10000
            goto 22
         else
            goto 10000       ! Unknown option
         endif
      endif
c
      if (standard .ne. ' ') then
         call bas_tag_lib(basis, osegment, tag, tag_in_lib, standard,
     $        filename, oshell_is_rel, .false.)
         goto 10
      else 
         write(luout,*) 'basis set name undefined'
         goto 10000
      endif
c     
 9000 return
c
10000 write(LuOut,1)
 1    format(' nuclear basis directive body format is:',/,
     $'       tag library <standard set> [file <filename>]',/,
     $       ' end basis')
      call errquit('bas_input_body: format error in the input', 0,
     &       INPUT_ERR)
c
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
