;+
; Function:
;   sdot.
;
; Purpose:
;   Calculate dot product of vec1 and vec2.
;
; Parameters:
;   vec1, in, type = array of certain type, required.
;     For a vector, dimension is [3].
;     For array of vectors, dimension is [n,3] default,
;     can be [3,n] if keyword transpose is set.
;     
;   vec2, in, type = array of certain type, required.
;     Same as vec1.
; 
; Keywords:
;   transpose = transpose, in, type = boolean, optional.
;     Set transpose when vec1 and vec2 are in [3,n].
;   
; Return:
;   return, out, type = array of certain type.
;     Has the same dimension as vec1.
;   
; Example:
;   l = sdot(r, v).
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
;   2012-09-19, Sheng Tian, add document.
;   2012-09-20, Sheng Tian, add transpose.
;-

function sdot, vec1, vec2, transpose = transpose

  compile_opt idl2

  on_error, 2

  vecsize = size(vec1)
  
  ; scalar or one m-dimension vector.
  if vecsize[0] le 1 then $
    return, total(vec1*vec2)
  
  ; vec1 and vec2 are in [n, 3].  
  if n_elements(transpose) eq 0 then begin
    vec = vec1[*,0] * vec2[*,0] + $
          vec1[*,1] * vec2[*,1] + $
          vec1[*,2] * vec2[*,2]
  ; vec1 and vec2 are in [3, n].
  endif else begin
    vec = vec1[0,*] * vec2[0,*] + $
          vec1[1,*] * vec2[1,*] + $
          vec1[2,*] * vec2[2,*]
  endelse
  
  return, vec

end