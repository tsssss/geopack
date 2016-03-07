;+
; Type:
;   function.
;
; Name:
;   sleap.
;
; Purpose:
;   Tell if given year(s) is leap year(s).
;
; Parameters:
;   yr, in, type = int/intarr[n], required.
;
; Keywords:
;   none.
;
; Return:
;   return, out, type = byte/bytarr[n].
;       1: leap year, 0: non-leap year.
;
; Example:
;   print, isleap([1998,2000,2004]).
;
; Notes:
;   none.
; 
; Dependence:
; 	none.
; 
; Author:
; 	Sheng Tian.
; 
; History:
; 	2012-07-06, Sheng Tian, create.
;-

function sleap, yr

    return, yr mod 400 eq 0 or yr mod 100 ne 0 and yr mod 4 eq 0

end
