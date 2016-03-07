;+
; Type: function.
; Purpose: Interpolate given data of at abscissa x to new x.
; Parameters:
;   data, in, [*], req. For scalar array, dimension is [n0].
;     For m dimension vector array, dimension is [m, n0] or [n0, m].
;   oldabs, in, [n0], req. Old abscissa.
;   newabs, in, [n1], req. New abscissa.
; Keywords:
;   _extra = extra, in, struct, opt. Keywords for interpol, see idl help.
;     /LSQuadratic, /NaN, /Quatradic, /Spline.
; Return: [*]. If data is array of [n0], then return [n1];
;     if data is array of [m, n0], then return, [m, n1].
;     if data is array of [n0, m], then return, [n1, m].
; Notes: Do not work if given data has more than 2 dimensions. Do NOT set /NaN 
;   automatically for interpol().
; History:
;   2011-07-20, Sheng Tian, create.
;   2011-08-16, Sheng Tian, add throw exterpolate.
;   2012-09-17, Sheng Tian, auto deal with dims.
;-

function sinterpol, data, oldabs, newabs, _extra = extra
  compile_opt idl2 & on_error, 0
  
  oldyy = data
  oldxx = oldabs
  newxx = newabs
  
  oldsize = size(oldyy)
  oldndims = oldsize[0]
  
  if oldndims gt 2 then $
    message, '# of dimension greater than 2 ...'
      
  if oldndims eq 0 then begin
    message, 'data is scalar, return ...', /continue
    return, data
  endif
  
  ; scalar array.
  if oldndims eq 1 then $
    return, interpol(oldyy, oldxx, newxx, _extra = extra)
  
  ; vector array.
  ; data in [n0, m].
  if oldsize[1] eq n_elements(oldxx) then begin
    newdims = [n_elements(newxx), oldsize[2]]
    newyy = replicate(oldyy[0], newdims)    ; keep type.
    for ii = 0, newdims[1]-1 do $
      newyy[*,ii] = interpol(oldyy[*,ii], oldxx, newxx, _extra = extra)
  ; data in [m, n0].
  endif else begin
    if oldsize[2] ne n_elements(oldxx) then $
      message, 'data incorrect dimension ...'
    newdims = [oldsize[1], n_elements(newxx)]
    newyy = replicate(oldyy[0], newdims)
    for ii = 0, newdims[0]-1 do $
      newyy[ii,*] = interpol(oldyy[ii,*], oldxx, newxx, _extra = extra)
  endelse
    
  return, newyy
end
