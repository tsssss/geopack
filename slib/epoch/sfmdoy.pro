;+
; Type: function.
; Purpose: Get year, month, day from year and day of year.
; Parameters:
;   yr, in, int, req. Year. If omit, treat to be non-leap year.
;   doy, in, int, req. Day of year.
; Keywords: none.
; Return: intarr[2]. [month, day].
; Notes: none.
; Dependence: none.
; History:
;   2012-07-06, Sheng Tian, create.
;-
function sfmdoy, yr, doy
    days = [0,31,59,90,120,151,181,212,243,273,304,334]
    if yr mod 400 eq 0 or yr mod 100 ne 0 and yr mod 4 eq 0 then days[2:*] += 1
    mo = (where(doy le days))[0]
    dy = doy - days[mo-1]
    return, [mo,dy]
end
