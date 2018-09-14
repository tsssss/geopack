;+
; Function: sdatarate.
; Purpose: Get data rate of a time series.
; Parameters:
;   t0s, in, dblarr(n), req. Time series in sec (e.g., unix time).
; Keywords: none.
; Return: out, double. Data rate in sec.
; Notes: The input time series needs to be quasi-uniform, but can has gaps.
; Dependence: slib.
; History:
;   2012-09-28, Sheng Tian, create.
;-

function sdatarate, t0s
  compile_opt idl2
  on_error, 2
  
  ; uniq.
  t1s = suniq(t0s)
  
  ; check time series.
  sz = size(t1s)
  if sz[0] ne 1 then $                ; sz[0] is # of dim.
    message, 'time series needs to be in 1-d ...'
  if sz[2] ne 4 and sz[2] ne 5 then $  ; sz[2] is type code.
    message, 'time series needs to be float or double ...'
  
  
  ; get difference.
  diff = t1s[1:sz[1]-1] - t1s[0:sz[1]-2]
  
  ; check monotonically increase.
  idx = where(diff lt 0D)
  if idx[0] ne -1 then $
    message, 'time series is not monotonically increasing ...'
    
  ; compute predominant data rate.
  minrate = min(diff)
  rate = mean(diff[where(diff lt minrate*1.2)])
  return, rate
end