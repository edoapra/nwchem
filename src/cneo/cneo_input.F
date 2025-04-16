      subroutine cneo_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "case.fh"
#include "nwc_const.fh"
#include "stdio.fh"
      integer rtdb
      logical s
      character*256 test
      double precision tmp_idx
      integer quantum_nuc_idx(nw_max_atom) ! indices for quantum nuclei
      character*2 quantum_nuc_smb
      logical oidx
      logical osmb
      character*256 tmp,tmp1
      integer k
      integer i,nuc_num
      character*256 nuc_bas_name
      character*255 cdft_bas_input

      do i=1,10
         quantum_nuc_idx(i) = 0
      end do

      oidx = .false.
      osmb = .false.
      i=0

 100  if (.not. inp_read()) 
     $     call errquit('cdft_input: inp_read failed',0,INPUT_ERR)
      

      if (.not. inp_a(test))
     $     call errquit('cdft_input: failed to read keyword',0,
     &       INPUT_ERR)


      if (inp_compare(.false.,'quantum_nuc',test)) then
 110     if (inp_f(tmp_idx)) then
            oidx = .true.
            i = i+1
            quantum_nuc_idx(i)=int(tmp_idx)
            go to 110
         else if (inp_a(quantum_nuc_smb)) then
            osmb = .true.
         else
            write(luout,*) 'CNEO: No quantum nucleus is defined.'
         end if

         tmp = 'cneo:quantum_nuc:'
         k = inp_strlen(tmp)+1
         if (oidx) then
            nuc_num = i
            tmp(k:) = ' '
            tmp(k:) = 'idx'
            s = rtdb_put(rtdb,tmp,mt_int,i,quantum_nuc_idx)
            tmp(k:) = ' '
            tmp(k:) = 'nuc_num'
            s = s .and. rtdb_put(rtdb,tmp,mt_int,1,nuc_num)
         else if (osmb) then
            tmp(k:) = ' '
            tmp(k:) = 'symbol'
            s = rtdb_cput(rtdb,tmp,1,quantum_nuc_smb)
         end if
         go to 100
      else if (inp_compare(.false.,'end',test)) then
         go to 900
      else if (inp_compare(.false.,'basis',test)) then
         call cneo_bas_input(rtdb)
         go to 100
      else if (inp_compare(.false.,'dft',test)) then
         call dft_input(rtdb)
         go to 100
      else if (inp_compare(.false.,'tddft',test)) then
         call tddft_input(rtdb)
         go to 100
      end if

 900  continue

      end