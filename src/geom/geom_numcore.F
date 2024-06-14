C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     geom_num_core -- return the number of orbitals considered to be
C     in the "core" in a system
C
C SYNOPSIS
C> \ingroup geom
C> @{
c
C> \brief Calculate the number of core orbitals in a geometry
C
C> Determines the number of core orbitals in a system based
C> simply on the constituent atoms and the standard general chemistry
C> ideas of core and valence orbitals.
C
C> \return Returns .true. if a valid geometry is provided, and
C> .false. otherwise.
c
      logical function geom_num_core(rtdb, geom, module, ncore)
*
* $Id$
*
      Implicit NONE
C
C DESCRIPTION
C     Determines the number of core orbitals in a system based on the
C     simply the constituent atoms and the standard general chemistry
C     ideas of core and valance.
C
C BUGS
C     Without any knowledge of the basis we have no choice but to 
C     assume the physical 5d, 7f shells, but if the calculation is
c     done with cartesian representations of these shells, this approach
c     will underestimate the actual number of core orbitals in the
c     calculation.
C     HvD: I do not believe this is a problem. The way how the physical
c          orbitals are represented does not impact on which orbitals
c          are occupied and which ones are not. The only case this fails
c          is if you limit the basis set so severely that electrons
c          are forced into the contaminant orbitals. But by that point
c          you have left any sense of physics well behind anyway...
C
C INCLUDE FILES
#include "nwc_const.fh"
#include "periodic.fh"
#include "geomP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
C 
C EXTERNAL ROUTINES
      Logical Geom_Check_Handle, Geom_Tag_to_Element
      External Geom_Check_Handle, Geom_Tag_to_Element

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer I, AtNr ! Atomic number
      Character*2 Sy  ! Atomic symbol
      Character*16 El ! Element name
      Logical Is_At   ! Is this tag an atom or something else?
      
      Integer num_elem
      Integer tagfrz(60)
      Integer numfrz(60)
      Integer rtdb         !< [Input] the RTDB handle
      Integer Geom         !< [Input] the geometry handle
      character*(*) module !< [Input] the code module
      Integer NCore        !< [Output] the number of core orbitals
      Integer tag
      character*80 key 
      character*20 space
      Integer mlen
#include "geom_numcore.fh"      
c
c     rjh: To avoid problems with ECPs and atoms that aren't really
c     .    atoms 'cos the user is messing with their charge, simply
c     .    ignore atoms that don't have a charge = atomic number.
C
     
      space = 'number tag elements'
      key = ' '
      mlen = inp_strlen(module) 
      write (key, '(a,'':'',a,'':'',a)') module(1:mlen),
     $    space(1:inp_strlen(space)), 'freeze by atoms'

      if(rtdb_get(rtdb, key, mt_int, 1, num_elem)) then
          
          space = 'tag freeze'
          key = ' '
          write (key, '(a,'':'',a,'':'',a)') module(1:mlen),
     $        space(1:inp_strlen(space)), 'freeze by atoms' 
      
          if(rtdb_get(rtdb, key, mt_int, num_elem, tagfrz)) then
            space = 'num freeze'
            write (key, '(a,'':'',a,'':'',a)') module(1:mlen),
     $        space(1:inp_strlen(space)), 'freeze by atoms'
            
             if(rtdb_get(rtdb, key, mt_int, num_elem, numfrz)) then
              do i = 1, num_elem
                tag = tagfrz(i)
                PerTbl(tag) = numfrz(i)
              enddo
            endif
          endif
       endif  

      Geom_Num_Core = Geom_Check_Handle(Geom, 'Geom_Num_Core') 
      If ( .NOT. Geom_Num_Core) Return
C
      NCore = 0
      Do I = 1, NCenter(Geom)
         Is_At = Geom_Tag_To_Element( Tags(i, Geom), Sy, El, AtNr)
         If ( Is_At ) then
            if (abs(charge(i,geom)-dble(atnr)) .lt. 1d-6) then
               NCore = NCore + PerTbl( AtNr )
            endif
         endif
      EndDo
C
      Return
      End
C> @}
      integer function geom_ncore_element(atnr)
      Implicit NONE
      integer atnr ! atomic number [in]
#include "nwc_const.fh"
#include "periodic.fh"
#include "geom_numcore.fh"
      if(atnr.lt.0.or.atnr.gt.nelements) call
     E     errquit(' geom_ncore_element: wrong z ',atnr,0)
      geom_ncore_element=pertbl(atnr)
      return
      end
