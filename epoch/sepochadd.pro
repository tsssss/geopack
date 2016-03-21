;+
; Type: function.
; Purpose: Add certain time to given epoch.
; Parameters:
;   et0, in, double, required. Given epoch in double.
;   num, in, int, required. The quantity to add. Negative is substract.
;   type, in, string, optional. Can be 'yr','mo','dy','hr','mi','sc','msc'.
; Keywords: none.
; Return: out, double. Output epoch.
; Notes: none.
; Dependence: none.
; History:
;   2013-06-18, Sheng Tian, create.
;-
function sepochadd, et0, num, type
    compile_opt idl2
    cdf_epoch, et0, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
    if n_elements(type) eq 0 then type = 'dy'
    case type of
        'yr': yr += num
        'mo': mo += num
        'dy': dy += num
        'hr': hr += num
        'mi': mi += num
        'sc': sc += num
        'msc': msc += num
    endcase
    cdf_epoch, et, yr, mo, dy, hr, mi, sc, msc, /compute_epoch
    return, et
end
