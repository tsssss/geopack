;+
; Function:
;   srmdim.
;
; Purpose:
;   I remove the selected dimension of given array.
;   Given no varyflags, I shrink unnecessary dimensions.
;    
; Parameters:
;   array0, in, type = array of certain type, required.
;     The array to be shrinked.
;     
;   varyflag, in, type = intarr, optional.
;     Set 0 on dimension you don't want.
;     e.g., array = intarr(2,1,3), varyflag = [0, 1, 1],
;     then the first dim is shrinked, return intarr(1,3).
;
; Keywords:
;   none.
;
; Return:
;   return, out, type = array of certain type.
;     The reduced array.
;
; Example:
;   newarray = srmdim(intarr(1,3)).
;   newarray = srmdim(intarr(1,3), [1,0])
; 
; Notes:
;   none.
;
; Dependence:
;   none.
;   
; Author:
;   Sheng Tian.
;
; History:
;   2011-05-15, Sheng Tian, create.
;   2012-07-02, Sheng Tian, revise.
;   2012-08-14, Sheng Tian, fix overflow.
;   2012-10-04, Sheng Tian, use product(), change name.  
;-

function srmdim, array0, varyflags

  compile_opt idl2
  
  on_error, 0
  
  varyndims = n_elements(varyflags)
  
  ; no varyflags, same as reform().
  if varyndims eq 0 then $
    return, reform(array0)
  
  ; old array's dimension info.
  oldsize = size(array0, /struct)
  oldndims = oldsize.n_dimensions
  ; scalar, nothing to be reduced.
  if oldndims eq 0 then return, array0
  if varyndims ne oldndims then $
    message, 'array dim and varyflags do not match ...'
  olddims = oldsize.dimensions[0:oldndims-1]
  
  ; new array's dimension info.
  vrydims = where(varyflags eq 1, newndims)  
  ; shrink to scalar.
  if newndims eq 0 then return, array0[0]
  ; no shrink.
  if newndims eq oldndims then return, array0
  newdims = olddims[vrydims]
  
  ; new array.
  array1 = make_array(newdims, type = oldsize.type, /nozero)
  
  ; calc the base of old and new array.
  oldbase = product([1,olddims[0:oldndims-2]], $
    /integer, /cumulative)
  newbase = oldbase[vrydims]
  
  ; loop through the new array.
  for i1 = 0ULL, n_elements(array1) - 1ULL do begin
    ; find indices of element i1 in the new array.
    idx1 = ulon64arr(newndims)
    tmpt = i1
    for jj = 0ULL, newndims - 1ULL do begin
      idx1[jj] = tmpt mod newdims[jj]
      tmpt = tmpt / newdims[jj]
    endfor
    i0 = 0ULL   ; 1-d index in old array.
    for jj = 0ULL, n_elements(idx1) - 1ULL do $
      i0 = i0 + newbase[jj] * idx1[jj]
    array1[i1] = array0[i0]
  endfor

  return, array1

end

arr = make_array([4,1,3,5], /nozero)
flg = [1,1,1,1]
help, arr
help, srmdim(arr, flg)
help, sshrinkmatrix(arr, flg)
end

