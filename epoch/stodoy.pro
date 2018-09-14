;+
; Type: function.
; Purpose: Get day of year from year, month, day.
; Parameters:
;   yr, in, int, req. Year, used to tell leap year.
;   mo, in, int, req. Month.
;   dy, in, int, req. Day.
; Keywords: none.
; Return: int. Day of year.
; Notes: none.
; Dependence: none.
; History:
;   2012-07-06, Sheng Tian, create.
;-
function stodoy, yr, mo, dy
    days = [0,31,59,90,120,151,181,212,243,273,304,334]
    if yr mod 400 eq 0 or yr mod 100 ne 0 and yr mod 4 eq 0 then days[2:*] += 1
    return, days[mo-1]+dy
end
