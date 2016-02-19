;+
; Type: pro.
; Purpose: Set current device to true color (decomposed=1).
; Parameters: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   2014-04-08, Sheng Tian, create.
;-
pro sgindexcolor, dev0, ct = ct, _extra = ex
    on_error, 2
    
    if n_elements(dev0) eq 0 then dev = !d.name
    
    case size(dev0,/type) of
        0: dev = !d.name
        7: dev = dev0
        else: begin
            dev = !d.name
            ct = dev0
        end
    endcase
    
    set_plot, dev
    
    dev = strlowcase(dev)
    case dev of
        'ps': device, decomposed = 0, /color, bits_per_pixel = 8
        'z': device, decomposed = 0, set_pixel_depth = 8
        'x': device, decomposed = 0
        'win': device, decomposed = 0
    endcase
    
    if keyword_set(ct) then begin
        if ct gt 41 then loadct2, ct, _extra = ex else loadct, ct, _extra = ex
    endif

end