;+
; Function:
;    sunitvec.
;
; Purpose:
;   I calculate unit vector of vector.
;    
; Parameters:
;   vec0, in, type = array of certain type, required.
;     For a vector, dimension is [m].
;     For array of vectors, dimension is [n,m] default.
;     can be [m,n] if keyword transpose is set.
;       
; Keywords:
;   transpose = transpose, in, type = boolean, optional.
;     Set transpose when vec0 is in [m, n].
;
; Return:
;   return, out, type = dblarr.
;     Same dimension as vec0.
;
; Example:
;   unit = sunitvec([1,2,3])
;
; Dependence:
;   none.
;
; Notes:
;   none.
;   
; Author:
;   Sheng Tian.
;
; History:
;   2012-10-01, Sheng Tian, create.
;-

function sunitvec, vec0, transpose = transpose

  compile_opt idl2
 
  on_error, 2
  
  vecsize = size(vec0)
  
  ; scalar or one m-dimension vector.
  if vecsize[0] le 1 then $
    return, double(vec0) / sqrt(total(vec0^2))
  
  ; array of m-dimension vector.
  if vecsize[0] gt 2 then $
    message, 'too many dimensions ...'
  
  ; vec0 in [n, m].  
  if n_elements(transpose) eq 0 then begin
    vec = double(vec0)
    for ii = 0, vecsize[1] - 1 do $
      vec[ii,*] = vec0[ii,*] / sqrt(total(vec0[ii,*]^2))
  ; vec0 in [n, m].
  endif else begin
    vec = double(vec0)
    for ii = 0, vecsize[2] - 1 do $
      vec[*,ii] = vec0[*,ii] / sqrt(total(vec0[*,ii]^2))
  endelse
  
  return, vec

end
