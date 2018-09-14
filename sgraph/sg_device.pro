;+
; Type: procedure.
; Purpose: Change device and set properties.
; Parameters:
;   dev, in, string, opt. The device to switch to. Do not switch device
;       if omitted.
; Keywords:
;   info, in, struct, opt. This is the modified !d structure, see sgopen.
;   options, in, struct, opt. This is the keywords accepted by device.
;   set, in, boolean, opt. Set this to change device property.
; Notes: none.
; Dependence: none.
; History:
;   2015-11-06, Sheng Tian, create.
;   2016-03-21, Sheng Tian, document.
;-

pro sg_device, dev, info = d0, options = opt, set = set

    on_error, 0

    ; all known devices.
    devp = 'ps'
    devz = 'z'
    devw = (!version.os_family eq 'unix')? 'x': 'win'

    ; get the original device and target device.
    dev0 = !d.name
    dev1 = (n_elements(dev) eq 0)? dev0: dev

    dev0 = strlowcase(dev0)
    dev1 = strlowcase(dev1)
    
    px2cm = 0.025d

    set_plot, dev1

    if keyword_set(set) then begin
        if n_elements(opt) eq 0 then begin      ; optins not set, use info.
            if n_elements(d0) eq 0 then message, 'no info ...'
            if size(d0,/type) ne 8 then message, 'info must be a struct ...'
            case dev1 of
                devp: opt = {set_character_size:[d0.x_ch_size,d0.y_ch_size], $
                    xsize:d0.x_size*px2cm, ysize:d0.y_size*px2cm, inch:0}
                devz: opt = {set_character_size: [d0.x_ch_size,d0.y_ch_size], $
                    set_resolution: [d0.x_size,d0.y_size]}
                devw: opt = {set_character_size: [d0.x_ch_size,d0.y_ch_size]}
            endcase
        endif
        set_plot, dev1
        device, _extra = opt
        if n_elements(d0) ne 0 then $
            sg_colormode, decomposed = d0.decomposed, depth = d0.depth, /set
    endif else begin
        sg_colormode, decomposed = dec, depth = dep
        d0 = create_struct(!d, 'decomposed',dec, 'depth',dep)
    endelse

    if dev1 ne dev0 then set_plot, dev0
    
end
