;+
; Type: <+++>.
; Purpose: <+++>.
; Parameters: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   <+yyyy-mm-dd+>, Sheng Tian, create.
;-


pro sgopen, fn0, xsize = xsize, ysize = ysize, cm = cm, inch = inch, $
    magnify = magc0, background = bgcolor, _extra = extra
    
    compile_opt idl2
    
    if n_elements(fn0) eq 0 then fn0 = 0    ; w-mode, window 0.
    fn = fn0[0]

    ; all known devices.
    devp = 'ps'
    devz = 'z'
    devw = (!version.os_family eq 'unix')? 'x': 'win'

    ; target device.
    if size(fn,/type) eq 7 then begin       ; string.
        dir = file_dirname(fn)
        if ~file_test(dir,/directory) then file_mkdir, dir
        extpos = strpos(fn,'.',/reverse_search)
        ext = strmid(fn,extpos+1)
        case strlowcase(ext) of
            'pdf': dev0 = devp
            'ps':  dev0 = devp
            'eps': dev0 = devp
            'png': dev0 = devz
            'jpg': dev0 = devz
            'jpeg':dev0 = devz
            else: message, 'does not support the extension '+ext+'...'
        endcase
    endif else dev0 = devw                  ; assume to be window id.


; **** save original status.
    
    ; save status for original device and target device.
    sg_device, !d.name, info = d0
    sg_device, dev0, info = d1

    ; save plot settings and colors.
    common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr
    c0 = { $
        pr_orig:ptr_new(r_orig), pr_curr:ptr_new(r_curr), $
        pg_orig:ptr_new(g_orig), pg_curr:ptr_new(g_curr), $
        pb_orig:ptr_new(b_orig), pb_curr:ptr_new(b_curr)}
    p0 = !p


; **** setup !sgraph.

    ; magnify?
    if n_elements(magc0) ne 0 then magc = double(magc0) else magc = 1d
    
    ; window size, in double, in pixel. Use default size if needed.
    if n_elements(xsize) eq 0 then xsz = !d.x_size else xsz = xsize
    if n_elements(ysize) eq 0 then ysz = !d.y_size else ysz = ysize
    xsz = double(xsz)*magc
    ysz = double(ysz)*magc

    ; convert from cm or inch.
    cm2px = 40d & in2cm = 2.54d & in2px = 101.6d    ; in2px = in2cm*cm2px.
    px2cm = 0.025d
    if keyword_set(cm) then begin
        if keyword_set(inch) then message, 'both inch and cm are set ...'
        xsz*= cm2px & ysz*= cm2px
    endif else if keyword_set(inch) then begin
        xsz*= in2px & ysz*= in2px
    endif
    
    ; char size, in double, in pixel.
    pchsz = [216d,360]*magc & zchsz = [9d,15]*magc & wchsz = [9d,15]*magc
    
    
    ; prepare window id, filename for 'ps' and 'z'.
    if dev0 ne devw then begin
        set_plot, devw
        window, /free, /pixmap
        wid = !d.window
        wdelete, wid
        if dev0 eq devp then begin
            pfn = fn & zfn = strmid(fn,0,extpos)+'png'
        endif else begin
            zfn = fn & pfn = strmid(fn,0,extpos)+'eps'
        endelse
    endif else begin
        wid = fn
        pfn = filepath('idl.eps', root_dir = shomedir())
        zfn = filepath('idl.png', root_dir = shomedir())
    endelse
    
    ; prepare sgraph structure.
    pmode = {sgid:pfn, device:devp, char:pchsz, thick:2d, area:[xsz,ysz]*px2cm}
    zmode = {sgid:zfn, device:devz, char:zchsz, thick:1d, area:[xsz,ysz]}
    wmode = {sgid:fix(wid), device:devw, char:zchsz, thick:1d, area:[xsz,ysz]}

    ; set up sgraph structure.
    defsysv, '!sgraph', exists = flag
    if flag eq 0 then begin     ; undefined.
        sgraph = {d0:d0, d1:d1, p0:p0, c0:c0, $
            pmode:pmode, zmode:zmode, wmode:wmode}
        defsysv, '!sgraph', sgraph
    endif else begin
        !sgraph.pmode = pmode
        !sgraph.zmode = zmode
        !sgraph.wmode = wmode
    endelse
    
; **** setup new device.

    case dev0 of
        devp: begin
            mode = !sgraph.pmode
            opt = {filename:mode.sgid, set_character_size:mode.char, $
                xsize:mode.area[0], ysize:mode.area[1], inches:0, $
                encapsulate:(ext ne 'ps'), copy:1, interpolate:1}
        end
        devz: begin
            mode = !sgraph.zmode
            opt = {set_character_size:mode.char, set_resolution:mode.area}
        end
        devw: begin
            mode = !sgraph.wmode
            opt = {set_character_size:mode.char, retain:2, true_color:24}
        end
    endcase

    ; set new device.
    sg_device, mode.device, options = opt, /set
    set_plot, mode.device
    
    !p.font = 1
    !p.thick = mode.thick
    if dev0 ne devp then begin
        !p.background = sgcolor('white')
        !p.color = sgcolor('black')
    endif

    ; set color.
    sgtruecolor

    ; draw window for window mode.
    if dev0 eq devw then $
        window, mode.sgid, xsize = mode.area[0], ysize = mode.area[1]

    ; erase.
    if n_elements(bgcolor) eq 0 then bgcolor = sgcolor('white')
    polyfill, [0,1,1,0,0], [0,0,1,1,0], /normal, color = bgcolor

end

fn = shomedir()+'/idl.png'
fn = 0
sgopen, fn, xsize = 400, ysize = 700, background = sgcolor('yellow')
xyouts, 10, 200, /device, 'Hello World!', color = sgcolor('red')
end
