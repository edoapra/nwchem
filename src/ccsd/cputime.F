      subroutine cputime(time)
      Implicit NONE
#include "util.fh"      
C     $Id$
      real*8 time
c#include "tcgmsg.fh"
c     time=tcgtime()
      time=util_wallsec()
      return
      end
