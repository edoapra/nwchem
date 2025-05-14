      subroutine cneo_nuc_init(rtdb, geom)
      implicit none
#include "errquit.fh"
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
        if(.not.rtdb_get(rtdb,tmp,mt_int,nuc_num,nuc_cent)) call 
     +      errquit('cneo_nuc_init:can not load nuc_cent',0,RTDB_ERR)
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
          if(.not.rtdb_put(rtdb,tmp,mt_int,1,ii)) call errquit
     +         ('cneo_nuc_init: save nuc_num?',0,RTDB_ERR)
          tmp(k:) = ' '
          tmp(k:) = 'idx'
          if(.not.rtdb_put(rtdb,tmp,mt_int,nuc_num,nuc_cent)) call 
     +         errquit('cneo_nuc_init:save idx?',0,RTDB_ERR)
          go to 120
        else
          write(6,*) 'Warning: no quantum idx/symbol defined'
          write(6,*) 'Hydrogen atoms are set by default'
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

      if(.not.geom_create(geom_cneo, 'geometry_cneo')) 
     +   call errquit('cneo_geom_init:geom_create?',0,0)
      if(.not.geom_rtdb_load(rtdb, geom_cneo, 'geometry'))
     +   call errquit('cneo_geom_init:geom_load?',0,0)
c  read quantum_nuc
      if(.not.geom_ncent(geom_cneo,nat)) call errquit
     +   ('cneo_geom_init:geom_ncent?',0,0)
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
#include "global.fh"

      integer rtdb, geom_cneo
      integer bases_nuc(max_quantum_nuc)
      logical status
      integer i, iatom

      bases_nuc = -999
      do i=1,nuc_num
         iatom = nuc_cent(i)
         if (ga_nodeid().eq.0) then
            write(6,*), 'cneo_bas_init: building basis for atom',iatom
         endif
         if (.not. bas_create(bases_nuc(i), 'nuc basis'))
     $      call errquit('cneo init: bas create',0,0)
         call cneo_bas_load(rtdb,geom_cneo,bases_nuc(i),'nuc basis',
     $                                       iatom)
      enddo

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
#include "ccneo_constrain.fh"
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
     $      'cneo:int1er:x'//trim(str_i))
         g_int1e_r(2,i)=ga_create_atom_blocked(geom,basis,
     $      'cneo:int1er:y'//trim(str_i))
         g_int1e_r(3,i)=ga_create_atom_blocked(geom,basis,
     $      'cneo:int1er:z'//trim(str_i))

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

      subroutine cneo_int1er_tidy
      implicit none 
#include "errquit.fh"
#include "global.fh"
#include "nwc_const.fh"
#include "ccneo.fh"
#include "ccneo_constrain.fh"

      integer i,j
      character*32 pname

      pname = 'cneo_int1er_tidy: '

      do i=1,nuc_num
         do j=1,3
            if(.not. ga_destroy(g_int1e_r(j,i))) call errquit
     +         (pname//'could not destroy int1er',0,GA_ERR)
         enddo
      enddo

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