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
pro sgtruecolor, dev
    on_error, 2
    
    if n_elements(dev) then set_plot, dev
    
    dev = strlowcase(!d.name)
    case dev of
        'ps': device, decomposed = 1, /color, bits_per_pixel = 8
        'z': device, decomposed = 1, set_pixel_depth = 24
        'x': device, decomposed = 1
        'win': device, decomposed = 1
    endcase

end