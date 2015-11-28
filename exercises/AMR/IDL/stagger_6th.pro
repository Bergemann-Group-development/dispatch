; $Id: stagger_6th.pro,v 1.8 2008/10/12 20:54:38 aake Exp $

;***********************************************************************
;
;  6th order staggered derivatives, 5th order half zone staggering
;
;***********************************************************************
PRO init_stagger
END

;***********************************************************************
PRO init_derivs, dxa=dxa, dya=dya, dza=dza
@common_cdata
  if n_elements(dxa) gt 0 then dx=dxa
  if n_elements(dya) gt 0 then dy=dya
  if n_elements(dza) gt 0 then dz=dza
END

;***********************************************************************
FUNCTION xup, f
;
;  f is centered on (i-.5,j,k), fifth order stagger
;
      COMMON cabc, a, b, c      ; for plot title purposes only
;
;-----------------------------------------------------------------------
      s=size(f)
      c = 3./256.
      b = -25./256.
      a = .5-b-c
      if s(0) eq 1 then $
        return, a*(f          + shift(f,-1)) $
              + b*(shift(f,1) + shift(f,-2)) $
              + c*(shift(f,2) + shift(f,-3))
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, c*(shift(f,2,0,0) + shift(f,-3,0,0)) $
            + b*(shift(f,1,0,0) + shift(f,-2,0,0)) $
            + a*(f              + shift(f,-1,0,0))
END
;***********************************************************************
FUNCTION yup, f
;
;  f is centered on (i,j-.5,k)
;
;
;-----------------------------------------------------------------------
      s=size(f)
      c = 3./256.
      b = -25./256.
      a = .5-b-c
      if s(0) le 1 then return,f
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, c*(shift(f,0,2,0) + shift(f,0,-3,0)) $
            + b*(shift(f,0,1,0) + shift(f,0,-2,0)) $
            + a*(f              + shift(f,0,-1,0))
END
;***********************************************************************
FUNCTION zup, f
;
;  f is centered on (i,j,k-.5)
;
;
;-----------------------------------------------------------------------
      s=size(f)
      c = 3./256.
      b = -25./256.
      a = .5-b-c
      if s(0) le 1 then return,f
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, c*(shift(f,0,0,2) + shift(f,0,0,-3)) $
            + b*(shift(f,0,0,1) + shift(f,0,0,-2)) $
            + a*(f              + shift(f,0,0,-1))
END
;***********************************************************************
FUNCTION xdn, f
;
;  the result is centered on (i-.5,j,k)
;
;-----------------------------------------------------------------------
      s=size(f)
      c = 3./256.
      b = -25./256.
      a = .5-b-c
      if s(0) eq 0 then return,f
      if s(0) eq 1 then $
        return, a*(f           + shift(f,1)) $
              + b*(shift(f,-1) + shift(f,2)) $
              + c*(shift(f,-2) + shift(f,3))
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, c*(shift(f,-2,0,0) + shift(f,3,0,0)) $
            + b*(shift(f,-1,0,0) + shift(f,2,0,0)) $
            + a*(f               + shift(f,1,0,0))
END
;***********************************************************************
FUNCTION ydn, f
;
;  the result is centered on (i,j,k-.5)
;
;-----------------------------------------------------------------------
      s=size(f)
      c = 3./256.
      b = -25./256.
      a = .5-b-c
      if s(0) le 1 then return,f
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, c*(shift(f,0,-2,0) + shift(f,0,3,0)) $
            + b*(shift(f,0,-1,0) + shift(f,0,2,0)) $
            + a*(f               + shift(f,0,1,0))
END
;***********************************************************************
FUNCTION zdn, f
;
;  the result is centered on (i,j,k-.5)
;
;-----------------------------------------------------------------------
      s=size(f)
      c = 3./256.
      b = -25./256.
      a = .5-b-c
      if s(0) le 1 then return,f
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, c*(shift(f,0,0,-2) + shift(f,0,0,3)) $
            + b*(shift(f,0,0,-1) + shift(f,0,0,2)) $
            + a*(f               + shift(f,0,0,1))
END
;***********************************************************************
FUNCTION xup1, f
;
;  f is centered on (i-.5,j,k)
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return, 0.5*(f + shift(f,-1))
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, 0.5*(f + shift(f,-1,0,0))
END
;***********************************************************************
FUNCTION yup1, f
;
;  f is centered on (i,j-.5,k)
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return,f
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, 0.5*(f + shift(f,0,-1,0))
END
;***********************************************************************
FUNCTION zup1, f
;
;  f is centered on (i,j,k-.5)
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return,f
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      return, 0.5*(f + shift(f,0,0,-1))
END
;***********************************************************************
FUNCTION xdn1, f
;
;  the result is centered on (i-.5,j,k)
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return, 0.5*(f + shift(f,1))
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, 0.5*(f + shift(f,1,0,0))
END
;***********************************************************************
FUNCTION ydn1, f
;
;  the result is centered on (i,j-.5,k)
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return,f
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, 0.5*(f + shift(f,0,1,0))
END
;***********************************************************************
FUNCTION zdn1, f
;
;  the result is centered on (i,j,k-.5)
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return,f
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, 0.5*(f + shift(f,0,0,1))
END
;***********************************************************************
FUNCTION ddxup, f
;
;  X partial derivative
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if n_elements(dx) eq 0 then dx=1.
      c = (-1.+(3.^5-3.)/(3.^3-3.))/(5.^5-5.-5.*(3.^5-3))/dx
      b = (-1./dx-120.*c)/24.
      a = (1./dx-3.*b-5.*c)
      if s(0) eq 1 then $
        return, a *(shift(f,-1) - f) $
              + b *(shift(f,-2) - shift(f,1)) $
              + c *(shift(f,-3) - shift(f,2))
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, c *(shift(f,-3,0,0) - shift(f,2,0,0)) $
            + b *(shift(f,-2,0,0) - shift(f,1,0,0)) $
            + a *(shift(f,-1,0,0) - f)
END
;***********************************************************************
FUNCTION ddyup, f
;
;  Y partial derivative
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return, make_array(size=size(f))
      if n_elements(dy) eq 0 then dy=1.
      c = (-1.+(3.^5-3.)/(3.^3-3.))/(5.^5-5.-5.*(3.^5-3))/dy
      b = (-1./dy-120.*c)/24.
      a = (1./dy-3.*b-5.*c)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, c *(shift(f,0,-3,0) - shift(f,0,2,0)) $
            + b *(shift(f,0,-2,0) - shift(f,0,1,0)) $
            + a *(shift(f,0,-1,0) - f)
END
;***********************************************************************
FUNCTION ddzup, f
;
;  Z partial derivative
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return, make_array(size=size(f))
      if n_elements(dz) eq 0 then dz=1.
      c = (-1.+(3.^5-3.)/(3.^3-3.))/(5.^5-5.-5.*(3.^5-3))/dz
      b = (-1./dz-120.*c)/24.
      a = (1./dz-3.*b-5.*c)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, c *(shift(f,0,0,-3) - shift(f,0,0,2)) $
            + b *(shift(f,0,0,-2) - shift(f,0,0,1)) $
            + a *(shift(f,0,0,-1) - f)
END
;***********************************************************************
FUNCTION ddxdn, f
;
;  X partial derivative
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if n_elements(dx) eq 0 then dx=1.
      c = (-1.+(3.^5-3.)/(3.^3-3.))/(5.^5-5.-5.*(3.^5-3))/dx
      b = (-1./dx-120.*c)/24.
      a = (1./dx-3.*b-5.*c)
      if s(0) eq 1 then $
        return, a*(f           - shift(f,1)) $
              + b*(shift(f,-1) - shift(f,2)) $
              + c*(shift(f,-2) - shift(f,3))
      return, c*(shift(f,-2,0,0) - shift(f,3,0,0)) $
            + b*(shift(f,-1,0,0) - shift(f,2,0,0)) $
            + a*(f               - shift(f,1,0,0))
END
;***********************************************************************
FUNCTION ddydn, f
;
;  Y partial derivative
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return, make_array(size=size(f))
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if n_elements(dy) eq 0 then dy=1.
      c = (-1.+(3.^5-3.)/(3.^3-3.))/(5.^5-5.-5.*(3.^5-3))/dy
      b = (-1./dy-120.*c)/24.
      a = (1./dy-3.*b-5.*c)
      return, c*(shift(f,0,-2,0) - shift(f,0,3,0)) $
            + b*(shift(f,0,-1,0) - shift(f,0,2,0)) $
            + a*(f - shift(f,0,1,0))
END
;***********************************************************************
FUNCTION ddzdn, f
;
;  Z partial derivative
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return, make_array(size=size(f))
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if n_elements(dz) eq 0 then dz=1.
      c = (-1.+(3.^5-3.)/(3.^3-3.))/(5.^5-5.-5.*(3.^5-3))/dz
      b = (-1./dz-120.*c)/24.
      a = (1./dz-3.*b-5.*c)
      return, c*(shift(f,0,0,-2) - shift(f,0,0,3)) $
            + b*(shift(f,0,0,-1) - shift(f,0,0,2)) $
            + a*(f - shift(f,0,0,1))
END
;***********************************************************************
FUNCTION ddx, f
;
;  X partial derivative
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if n_elements(dx) eq 0 then dx=1.
      a =  3./(4.*dx)
      b = -3./(20.*dx)
      c =  1./(60.*dx)
      if s(0) eq 1 then $
        return, a*(shift(f,-1) - shift(f,1)) $
              + b*(shift(f,-2) - shift(f,2)) $
              + c*(shift(f,-3) - shift(f,3))
      return, a*(shift(f,-1,0,0) - shift(f,1,0,0)) $
            + b*(shift(f,-2,0,0) - shift(f,2,0,0)) $
            + c*(shift(f,-3,0,0) - shift(f,3,0,0))
END
;***********************************************************************
FUNCTION ddy, f
;
;  Y partial derivative
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return, make_array(size=size(f))
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if n_elements(dy) eq 0 then dy=1.
      a =  3./(4.*dy)
      b = -3./(20.*dy)
      c =  1./(60.*dy)
      return, a*(shift(f,0,-1,0) - shift(f,0,1,0)) $
            + b*(shift(f,0,-2,0) - shift(f,0,2,0)) $
            + c*(shift(f,0,-3,0) - shift(f,0,3,0))
END
;***********************************************************************
FUNCTION ddz, f
;
;  Z partial derivative
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return, make_array(size=size(f))
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if n_elements(dz) eq 0 then dz=1.
      a =  3./(4.*dz)
      b = -3./(20.*dz)
      c =  1./(60.*dz)
      return, a*(shift(f,0,0,-1) - shift(f,0,0,1)) $
            + b*(shift(f,0,0,-2) - shift(f,0,0,2)) $
            + c*(shift(f,0,0,-3) - shift(f,0,0,3))
END
;***********************************************************************
FUNCTION difx3, f
;
;  Finite difference one zone left to two zones right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if s(0) eq 1 then return, shift(f,1) - shift(f,-2)
      return, shift(f,1,0,0) - shift(f,-2,0,0)
END
;***********************************************************************
FUNCTION dify3, f
;
;  Finite difference one zone left to two zones right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return, make_array(size=size(f))
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, shift(f,0,1,0) - shift(f,0,-2,0)
END
;***********************************************************************
FUNCTION difz3, f
;
;  Finite difference one zone left to two zones right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return, make_array(size=size(f))
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, shift(f,0,0,1) - shift(f,0,0,-2)
END
;***********************************************************************
FUNCTION difx2, f
;
;  Finite difference one zone left to two zones right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return, shift(f,1,0,0) - shift(f,-1)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, shift(f,1,0,0) - shift(f,-1,0,0)
END
;***********************************************************************
FUNCTION dify2, f
;
;  Finite difference one zone left to two zones right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return, make_array(size=size(f))
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, shift(f,0,1,0) - shift(f,0,-1,0)
END
;***********************************************************************
FUNCTION difz2, f
;
;  Finite difference one zone left to two zones right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return, make_array(size=size(f))
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, shift(f,0,0,1) - shift(f,0,0,-1)
END
;***********************************************************************
FUNCTION difx1, f
;
;  Finite difference to one zone right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return, f - shift(f,-1)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, f - shift(f,-1,0,0)
END
;***********************************************************************
FUNCTION dify1, f
;
;  Finite difference to one zone right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return, make_array(size=size(f))
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, f - shift(f,0,-1,0)
END
;***********************************************************************
FUNCTION difz1, f
;
;  Finite difference to one zone right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 1 then return, make_array(size=size(f))
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      return, f - shift(f,0,0,-1)
END

;***********************************************************************
function ddxup1, f
@common_cdata

      s=size(f)
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if n_elements(dx) eq 0 then dx=1.
      a = (1./dx)
      if s(0) eq 1 then $
      return, a*(shift(f,-1) - f )
      return, a*(shift(f,-1,0,0) - f )
end
;***********************************************************************
function ddyup1, f
@common_cdata

      s=size(f)
      if s(0) eq 1 then return, 0.*f
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if n_elements(dy) eq 0 then dy=1.
      a = (1./dy)
      return, a*(shift(f,0,-1,0) - f )
end
;***********************************************************************
function ddzup1, f
@common_cdata

      s=size(f)
      if s(0) lt 3 then return, 0.*f
      if n_elements(dz) eq 0 then dz=1.
      a = (1./dz)
      return, a*(shift(f,0,0,-1) - f )
end
;***********************************************************************
function ddxdn1, f
@common_cdata

      s=size(f)
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if n_elements(dx) eq 0 then dx=1.
      a = (1./dx)
      if s(0) eq 1 then $
      return, a*(f - shift(f,1))
      return, a*(f - shift(f,1,0,0))
end
;***********************************************************************
function ddydn1, f
@common_cdata

      s=size(f)
      if s(0) eq 1 then return, 0.*f
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if n_elements(dx) eq 0 then dy=1.
      a = (1./dy)
      return, a*(f - shift(f,0,1,0))
end
;***********************************************************************
function ddzdn1, f
@common_cdata

      s=size(f)
      if s(0) lt 3 then return, 0.*f
      if n_elements(dz) eq 0 then dz=1.
      a = (1./dz)
      return, a*(f - shift(f,0,0,1))
end
;***********************************************************************

PRO stagger_6th
  print,'stagger_6th initialized'
END
