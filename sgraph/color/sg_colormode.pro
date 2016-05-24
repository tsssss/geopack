;+
; Type: procedure.
; Purpose: Get or set decompose state and color depth to certain device.
; Parameters:
;   dev, in, string, opt. Device name whose color mode will be returned.
;       Default is the current device's name.
; Keywords:
;   decomposed, in/out, boolean. 1: true color, 0: index color.
;   depth, in/out, int, opt. Current device's color depth.
;   set, in, boolean, opt. Set to modify current device's color mode.
; Notes: Depth may have different meanings in different devices.
; Dependence: none.
; History:
;   2015-11-06, Sheng Tian, create.
;-

pro sg_colormode, dev, decomposed = dec, depth = depth, set = set
    on_error, 2

    ; get the original device and target device.
    dev0 = !d.name
    dev1 = (n_elements(dev) eq 0)? dev0: dev

    dev0 = strlowcase(dev0)
    dev1 = strlowcase(dev1)
    
    set_plot, dev1

    if keyword_set(set) then begin
        if n_elements(dec) eq 0 then dec = 0
        if n_elements(depth) eq 0 then depth = (dec eq 0)? 8: 24
        case dev1 of
            'ps': device, decomposed = dec, /color, bits_per_pixel = depth
            'z': device, decomposed = dec, set_pixel_depth = depth
            'x': device, decomposed = dec
            'win': device, decomposed = dec
            else: message, 'unknown device: '+dev1+' ...'
        endcase
    endif else begin
        case dev1 of
            'ps': begin
                if float(!versin.release) gt 7.1 then begin
                    device, get_decomposed = dec
                endif else begin
                    help, /device, output = tmp
                    dec = strpos(strupcase(tmp[4]), 'DECOMPOSED')) ne -1
                endelse
                if dec then depth = 24 else depth = 8
            end
            'z': begin
                device, get_decomposed = dec, get_pixel_depth = depth
            end
            'x': device, get_decomposed = dec, get_visual_depth = depth
            'win': device, get_decomposed = dec, get_visual_depth = depth
            else: begin
                dec = 0
                depth = 8
            end
        endcase
    endelse

    if dev1 ne dev0 then set_plot, dev0

end
