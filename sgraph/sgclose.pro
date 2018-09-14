;+
; Type: procedure.
; Purpose: Cloase a canvas.
; Parameters: none.
; Keywords:
;   wdelete, in, boolean, opt. Close window in 'x' or 'win'. Default
;       behaviour is to leave the window opened.
; Notes: none.
; Dependence: slib.
; History:
;   2015-11-06, Sheng Tian, create.
;   2016-03-21, Sheng Tian, create.
;-


pro sgclose, wdelete = wdelete

    ; all known devices.
    devp = 'ps'
    devz = 'z'
    devw = (!version.os_family eq 'unix')? 'x': 'win'

    ; save file or close window.
    dev1 = strlowcase(!d.name)
    case dev1 of
        devw: mode1 = !sgraph.wmode
        devp: mode1 = !sgraph.pmode
        devz: mode1 = !sgraph.zmode
    endcase
    ; filename, extension.
    fn = mode1.sgid
    if size(fn,/type) eq 7 then begin
        extpos = strpos(fn,'.',/reverse_search)
        ext = strmid(fn,extpos+1)
    endif

    ; close current device, save to file.
    case dev1 of

        devw: begin
            if keyword_set(wdelete) then wdelete, !d.window
        end

        devp: begin
            device, /close
            if strlowcase(ext) eq 'pdf' then begin
                tmp = strmid(fn,0,extpos)+'.eps'
                file_copy, fn, tmp, /overwrite
                sps2pdf, tmp, /rm
            endif
        end

        devz: begin
            device, get_pixel_depth = true
            true = (true eq 8)? 0: 1
            case strlowcase(ext) of
                'png': begin
                    if true then write_png, fn, tvrd(/true) else begin
                        tvlct, r, g, b, /get
                        write_png, fn, tvrd(), r, g, b
                    endelse & end
                'jpg': begin tvlct, r, g, b, /get
                    write_jpeg, fn, tvrd(true=3), true=3 & end
                'jpeg': begin tvlct, r, g, b, /get
                    write_jpeg, fn, tvrd(true=3), true=3 & end
                else: message, 'unsupported raster format ...'
            endcase
        end
    endcase

    ; restore !p.
    !p = !sgraph.p0

    ; restore current device's settings.
    sg_device, info = !sgraph.d1, /set
    
    ; restore original device's settings.
    if !sgraph.d0.name eq !sgraph.d1.name then return
    set_plot, !sgraph.d0.name
    sg_device, info = !sgraph.d0, /set

    ; restore color.
    common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr
    r_curr = !sgraph.d0.r0
    g_curr = !sgraph.d0.g0
    b_curr = !sgraph.d0.b0

end
