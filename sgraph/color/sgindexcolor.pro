;+
; Type: pro.
; Purpose: Set current device to index color (decomposed=0).
; Parameters:
;   dev0, in, string, opt. The device to be set to index color mode. Default is
;       the current device.
; Keywords:
;   ct, in, int, opt. The color table. Do not load color table if no id is set.
;   file, in, string, opt. Set the file name of color table files, then ct is
;       the color table id in that file.
; Notes: none.
; Dependence: none.
; History:
;   2014-04-08, Sheng Tian, create.
;-
pro sgindexcolor, dev0, ct = ct, file = file, _extra = ex
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
    if n_elements(ct) eq 0 then ct = 0
    if n_elements(file) ne 0 then begin
        ctfn = srootdir()+'/'+file+'.tbl'
        loadct, ct, file = ctfn, _extra = ex
    endif else loadct, ct, _extra = ex

end
