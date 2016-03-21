;+
; Type: pro.
; Purpose: Set current device to index color (decomposed=0).
; Parameters:
;   dev0, in, string, opt. The device to be set to index color mode. Default is
;       the current device.
; Keywords:
;   ct, in, int, opt. The color table. Do not load color table if no id is set.
;   ct2, in, int/boolean, opt. Set to use loadct2.
; Notes: none.
; Dependence: none.
; History:
;   2014-04-08, Sheng Tian, create.
;-
pro sgindexcolor, dev0, ct = ct, ct2 = ct2, _extra = ex
    on_error, 2
    
    case size(dev0,/type) of
        0: dev = !d.name        ; no device.
        7: dev = dev0           ; dev0 is device string.
        else: begin             ; dev0 is color table id.
            dev = !d.name
            ct = dev0
        end
    endcase

    ; change device and set to index color.
    set_plot, dev
    dev = strlowcase(dev)
    case dev of
        'ps': device, decomposed = 0, /color, bits_per_pixel = 8
        'z': device, decomposed = 0, set_pixel_depth = 8
        'x': device, decomposed = 0
        'win': device, decomposed = 0
    endcase
    
    ; load color table.
    if keyword_set(ct2) then begin
        if n_elements(ct) ne 0 then loadct2, ct, _extra = ex $
        else loadct2, ct2, _extra = ex
    endif else if keyword_set(ct) then loadct, ct, _extra = ex

end
