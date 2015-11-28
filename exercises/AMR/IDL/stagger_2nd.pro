;***********************************************************************
;
;  2nd order staggered derivatives, 1st order half zone staggering
;
;***********************************************************************
;***********************************************************************
PRO init_derivs, dxa=dxa, dya=dya, dza=dza
@common_cdata
  if n_elements(dxa) eq 1 then dx=dxa
  if n_elements(dya) eq 1 then dy=dya
  if n_elements(dza) eq 1 then dz=dza
END

FUNCTION xup, f
;
;  f is centered on (i-.5,j,k)
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if      s(0) eq 1 then return, 0.5*(f + shift(f,-1    )) $
      else if s(0) eq 2 then return, 0.5*(f + shift(f,-1,0  )) $
      else                   return, 0.5*(f + shift(f,-1,0,0))
END
;***********************************************************************
FUNCTION yup, f
;
;  f is centered on (i,j-.5,k)
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if      s(0) eq 1 then return, f $
      else if s(0) eq 2 then return, 0.5*(f + shift(f,0,-1  )) $
      else                   return, 0.5*(f + shift(f,0,-1,0))
END
;***********************************************************************
FUNCTION zup, f
;
;  f is centered on (i,j,k-.5)
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if      s(0) le 2 then return, f $
      else                   return, 0.5*(f + shift(f,0,0,-1))
END
;***********************************************************************
FUNCTION xdn, f
;
;  f is centered on (i-.5,j,k)
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if      s(0) eq 1 then return, 0.5*(f + shift(f,+1    )) $
      else if s(0) eq 2 then return, 0.5*(f + shift(f,+1,0  )) $
      else                   return, 0.5*(f + shift(f,+1,0,0))
END
;***********************************************************************
FUNCTION ydn, f
;
;  f is centered on (i,j,k-.5)
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if      s(0) eq 1 then return, f $
      else if s(0) eq 2 then return, 0.5*(f + shift(f,0,+1  )) $
      else begin
	ydn = shift(f,0,+1,0)
	ydn = 0.5*(f+temporary(ydn))
        return, ydn
      end
END
;***********************************************************************
FUNCTION zdn, f
;
;  f is centered on (i,j,k-.5)
;
@common_cdata
;
;-----------------------------------------------------------------------
      s=size(f)
      if      s(0) le 2 then return, f $
      else                   return, 0.5*(f + shift(f,0,0,+1))
END
;***********************************************************************
FUNCTION xup1, f
;
;  f is centered on (i-.5,j,k)
;
@common_cdata
;
;-----------------------------------------------------------------------
      return, xup(f)
END
;***********************************************************************
FUNCTION yup1, f
;
;  f is centered on (i,j-.5,k)
;
@common_cdata
;
;-----------------------------------------------------------------------
      return, yup(f)
END
;***********************************************************************
FUNCTION zup1, f
;
;  f is centered on (i,j,k-.5)
;
@common_cdata
;
;-----------------------------------------------------------------------
      return, zup(f)
END
;***********************************************************************
FUNCTION xdn1, f
;
;  f is centered on (i-.5,j,k)
;
@common_cdata
;
;-----------------------------------------------------------------------
      return, xdn(f)
END
;***********************************************************************
FUNCTION ydn1, f
;
;  f is centered on (i,j,k-.5)
;
@common_cdata
;
;-----------------------------------------------------------------------
      return, ydn(f)
END
;***********************************************************************
FUNCTION zdn1, f
;
;  f is centered on (i,j,k-.5)
;
@common_cdata
;
;-----------------------------------------------------------------------
      return, zdn(f)
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
      if      s(0) eq 1 then return, (1./dx)*(shift(f,-1    )-f) $
      else if s(0) eq 2 then return, (1./dx)*(shift(f,-1,0  )-f) $
      else                   return, (1./dx)*(shift(f,-1,0,0)-f)
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
      if      s(0) eq 1 then return, replicate(0.,s(1)) $
      else if s(0) eq 2 then return, (1./dy)*(shift(f,0,-1  )-f) $
      else                   return, (1./dy)*(shift(f,0,-1,0)-f)
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
      if      s(0) eq 1 then return, replicate(0.,s(1)) $
      else if s(0) eq 2 then return, replicate(0.,s(1),s(2)) $
      else                   return, (1./dz)*(shift(f,0,0,-1)-f)
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
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      if n_elements(dx) eq 0 then dx=1.
      return, (1./dx)*(f - shift(f,1,0,0))
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
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      if n_elements(dy) eq 0 then dy=1.
      return, (1./dy)*(f - shift(f,0,1,0))
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
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      if n_elements(dz) eq 0 then dz=1.
      return, (1./dz)*(f - shift(f,0,0,1))
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
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      if n_elements(dx) eq 0 then dx=1.
      return, (0.5/dx)*(shift(f,-1,0,0) - shift(f,1,0,0))
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
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      if n_elements(dy) eq 0 then dy=1.
      return, (0.5/dy)*(shift(f,0,-1,0) - shift(f,0,1,0))
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
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      if n_elements(dz) eq 0 then dz=1.
      return, (0.5/dz)*(shift(f,0,0,-1) - shift(f,0,0,1))
END
;***********************************************************************
FUNCTION difx3, f
;
;  Positive part of finite difference one zone left to two zones right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      return, shift(f,1,0,0) - shift(f,-2,0,0)
END
;***********************************************************************
FUNCTION dify3, f
;
;  Positive part of finite difference one zone left to two zones right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      return, shift(f,0,1,0) - shift(f,0,-2,0)
END
;***********************************************************************
FUNCTION difz3, f
;
;  Positive part of finite difference one zone left to two zones right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      return, shift(f,0,0,1) - shift(f,0,0,-2) 
END
;***********************************************************************
FUNCTION difx2, f
;
;  Positive part of finite difference one zone left to two zones right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      return, shift(f,1,0,0) - shift(f,-1,0,0)
END
;***********************************************************************
FUNCTION dify2, f
;
;  Positive part of finite difference one zone left to two zones right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      return, shift(f,0,1,0) - shift(f,0,-1,0)
END
;***********************************************************************
FUNCTION difz2, f
;
;  Positive part of finite difference one zone left to two zones right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      return, shift(f,0,0,1) - shift(f,0,0,-1) 
END
;***********************************************************************
FUNCTION difx1, f
;
;  Finite difference to one zone right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      return, f - shift(f,-1,0,0)
END
;***********************************************************************
FUNCTION dify1, f
;
;  Finite difference to one zone right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      return, f - shift(f,0,-1,0)
END
;***********************************************************************
FUNCTION difz1, f
;
;  Finite difference to one zone right.
;
;-----------------------------------------------------------------------
      s=size(f)
      if s(0) eq 2 then f=reform(f,s(1),s(2),1)
      if s(0) eq 1 then f=reform(f,s(1),1,1)
      return, f - shift(f,0,0,-1)
END
;***********************************************************************
PRO stagger_2nd
  print,'stagger_2nd initialized'
END

