;+
; Function:
;   snorm.
;
; Purpose:
;   Calculate magnitude of vector.
;
; Parameters:
;   vec, in, type = array of certain type, required.
;     For a vector, dimension is [m].
;     For array of vectors, dimension is [n,m] default,
;     can be [m,n] if keyword transpose is set.
; 
; Keywords:
;   transpose = transpose, in, type = boolean, optional.
;     Set transpose when vec1 and vec2 are in [m,n].
;   
; Return:
;   return, out, type = array of certain type.
;     [0] for vector in [m].
;     [n] for vector in [n,m] or [m,n].
;   
; Example:
;   r = snorm(dblarr(3,100), /transpose).
;   
; Notes:
;   * I do NOT check dimensions, keep track yourself.
;
; Dependence:
;   none.
;  
; Author:
;   Sheng Tian.
; 
; History:
;   2010-07-28, Sheng Tian, create.
;   2012-11-01, Sheng Tian, move to slib.
;-

function snorm, vec, transpose = transpose

  compile_opt idl2
  
  on_error, 2
  
  sz = size(vec)
  
  ; scalar.
  if sz[0] eq 0 then return, abs(vec)
  
  ; pure vector.
  if sz[0] eq 1 then return, sqrt(total(vec^2))
  
  ; vector array.
  if sz[0] ne 2 then message, 'wrong dimension ...'
 
  ; vec in [m, n].
  if n_elements(transpose) eq 0 then $
    return, sqrt(total(vec^2, 2)) $
  else return, sqrt(total(vec^2, 1))

end