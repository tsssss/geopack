;+
; Type: function.
; Purpose: Ceil epoch to year, month, day, etc.
; Parameters: 
;   et0, in, double, required. Input epoch.
;   type, in, string, optional. Can be 'yr','mo','dy','hr','mi','sc'.
; Keywords: none.
; Return: out, double. Output epoch.
; Notes: Do not work for array.
; Dependence: none.
; History:
;   2013-06-19, Sheng Tian, create.
;-
function sepochceil, et0, type
    compile_opt idl2
    cdf_epoch, et0, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
    if n_elements(type) eq 0 then type = 'dy'
    switch type of      ; order matters!
        'yr': mo = 0
        'mo': dy = 0
        'dy': hr = 0
        'hr': mi = 0
        'mi': sc = 0
        'sc': msc = 0
    endswitch
    case type of
        'yr': yr += 1
        'mo': mo += 1
        'dy': dy += 1
        'hr': hr += 1
        'mi': mi += 1
        'sc': sc += 1
        else: return, t0
    endcase
    cdf_epoch, et, yr, mo, dy, hr, mi, sc, msc, /compute_epoch
    return, et
end
