      subroutine nwchem_libxc_read(rtdb,funname,success,
     $           fac,lfac,nlfac)
#ifdef USE_LIBXC
      use,intrinsic :: iso_c_binding
      use xc_f03_lib_m
#endif
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "cdft.fh"
#include "case.fh"
#include "inp.fh"
#ifdef USE_LIBXC
#include "libxc.fh"
      type(xc_f03_func_t) :: xcfunc
      type(xc_f03_func_info_t) :: xcinfo
      integer(c_int) :: funcid, xcflags, xckind, xcfamily
      real(c_double) :: alpha
      integer(c_int),parameter :: error = -1
#endif
      double precision fact
      double precision fac(numfunc)
      logical lfac(numfunc),nlfac(numfunc)
      integer rtdb
      character*(*) :: funname
      logical :: success
      success = .false.

#ifdef USE_LIBXC

      funcid = xc_f03_functional_get_number(trim(funname))
      if (funcid.eq.error) return

      if (.not.(inp_f(fact))) fact=1d0

      call xc_f03_func_init(xcfunc,funcid,XC_UNPOLARIZED)

      xcinfo = xc_f03_func_get_info(xcfunc)
      xckind = xc_f03_func_info_get_kind(xcinfo)
      xcflags = xc_f03_func_info_get_flags(xcinfo)
      xcfamily = xc_f03_func_info_get_family(xcinfo)

      select case(xcfamily)
      case (XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
      case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      case (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
      case default
        call errquit('nwchem_libxc_read: unknown functional family',
     $                xcfamily,input_err)
      endselect

      select case(xcfamily)
      case (XC_FAMILY_HYB_LDA, XC_FAMILY_HYB_GGA, XC_FAMILY_HYB_MGGA)
        if (iand(xcflags,xc_flags_hyb_cam).eq.xc_flags_hyb_cam) then
          call xc_f03_hyb_cam_coef(xcfunc,cam_omega,cam_alpha,cam_beta)

          cam_exch = .true.
          if (cam_alpha.eq.0d0) then
            lfac(1) = .true.
            nlfac(1) = .true.
            fac(1) = cam_beta
            cam_beta = 1d0
            cam_srhf = .true.
          else
            cam_alpha = cam_alpha + cam_beta
            cam_beta = -cam_beta
            lfac(1) = .true.
            nlfac(1) = .true.
            fac(1) = 1.0d0
            cam_srhf = .false.
          endif
          call xc_setcamparam(rtdb,cam_exch,cam_srhf,cam_omega,
     $                        cam_alpha,cam_beta)
        else
          alpha = xc_f03_hyb_exx_coef(xcfunc)
          lfac(1) = .true.
          nlfac(1) = .true.
          fac(1) = alpha
        endif
      endselect

      if (iand(xcflags,xc_flags_hyb_camy).eq.xc_flags_hyb_camy) goto 10
      if (iand(xcflags,xc_flags_hyb_lcy).eq.xc_flags_hyb_lcy) goto 10
      if (iand(xcflags,xc_flags_vv10).eq.xc_flags_vv10) goto 20

      success = .true.
      libxcon = .true.

      libxc_nfuncs = libxc_nfuncs + 1
      libxc_funcs(libxc_nfuncs) = funcid
      libxc_family(libxc_nfuncs) = xcfamily
      libxc_kind(libxc_nfuncs) = xckind
      libxc_flags(libxc_nfuncs) = xcflags
      libxc_facts(libxc_nfuncs) = fact

      call xc_f03_func_end(xcfunc)

      return

   10 call errquit('nwchem_libxc_read: Yukawa kernel is not available',
     $              0,input_err) 
   20 call errquit('nwchem_libxc_read: VV10 kernel is not available',
     $              0,input_err) 
#endif
      end subroutine


