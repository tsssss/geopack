
pro sgopen, id0, xsize = xsize, ysize = ysize, cm = cm, inch = inch, $
    magnify = magc0, background = bgc0, _extra = extra
    
    compile_opt idl2
    
    if n_elements(id0) eq 0 then id0 = 0    ; w-mode, window 0.
    fn = id0[0]
    
    if n_elements(bgc0) eq 0 then bgc0 = sgcolor('white')
    
    
    ; save the old device's config.
    ; device and plot structure.
    d0 = !d & p0 = !p & x0 = !x & y0 = !y & z0 = !z
    sg_colormode, d0.name, decomposed = dec, depth = depth
    tvlct, r0, g0, b0, /get
    d0 = create_struct(d0, 'decomposed', dec, 'depth', long(depth), $
        'r0',r0, 'g0',g0, 'b0',b0)

    
    ; magnification coef.
    if n_elements(magc0) ne 0 then magc = double(magc0) else magc = 1d
    
    
    ; canvas size, in double, in pixel.
    if n_elements(xsize) eq 0 then xsz = !d.x_size else xsz = double(xsize)
    if n_elements(ysize) eq 0 then ysz = !d.y_size else ysz = double(ysize)
    
    ; convert to pixel from cm or inch.
;    cm2px = 40d & in2cm = 2.54d & in2px = 101.6d    ; in2px = in2cm*cm2px.
    px2cm = 0.025d
    if keyword_set(cm) and keyword_set(inch) then $
        message, 'cannot set incn and cm at the same time ...'
    case 1 of
        keyword_set(cm): tmp = 40d      ; cm2px.
        keyword_set(inch): tmp = 101.6d ; in2px = in2cm*2.54.
        else: tmp = 1d
    endcase
    
    ; apply magnfication to canvas.
    tmp *= magc
    xsz = double(xsz)*tmp
    ysz = double(ysz)*tmp

    
    ; char size, in double, in pixel.
    pchsz = [216d,360]*magc
    zchsz = [9d,15]*magc
    wchsz = [9d,15]*magc
    
    
    ; device names.
    devw = (!version.os_family eq 'unix')? 'x': 'win'
    devp = 'ps'
    devz = 'z'
    
    ; determine the current device.
    if size(fn,/type) eq 7 then begin       ; string.
        dir = file_dirname(fn)
        if dir eq '' then cd, current = dir
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
    endif else dev0 = devw                  ; window.
    
    
    ; set the default file name for all devices.
    if dev0 ne devw then begin
        set_plot, devw
        window, /free, /pixmap
        wid = !d.window
        wdelete, wid
        if dev0 eq devp then begin
            pfn = fn
            zfn = strmid(fn,0,extpos)+'png'
        endif else begin
            zfn = fn
            pfn = strmid(fn,0,extpos)+'eps'
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


    ; set up the sgraph structure.
    defsysv, '!sgraph', exists = flag
    if flag eq 0 then begin
;        sgraph = {d0:d0, d1:d0, p0:!p, x0:!x, y0:!y, z0:!z, $
;            r0:r0, g0:g0, b0:b0, pmode:pmode, zmode:zmode, wmode:wmode}
        sgraph = {d0:d0, d1:d0, p0:!p, x0:!x, y0:!y, z0:!z, $
            pmode:pmode, zmode:zmode, wmode:wmode}
        defsysv, '!sgraph', sgraph
    endif else begin
        !sgraph.pmode = pmode
        !sgraph.zmode = zmode
        !sgraph.wmode = wmode
    endelse
    
    ; set to the new device.
    case dev0 of
        devw: mode = !sgraph.wmode
        devp: mode = !sgraph.pmode
        devz: mode = !sgraph.zmode
    endcase

    if dev0 eq devp then $     ; encapsulate for 'pdf' and 'eps'.
        opt = {filename:mode.sgid, set_character_size:mode.char, $
            xsize:mode.area[0], ysize:mode.area[1], inches:0, $
            encapsulate:(ext ne 'ps'), copy:1, interpolate:1}
    if dev0 eq devz then $
        opt = {set_character_size:mode.char, set_resolution:mode.area}
    if dev0 eq devw then $
        opt = {set_character_size:mode.char, retain:2, true_color:24}
    
    set_plot, mode.device
    
    
    ; save the new device's original config.
    d1 = !d
    sg_colormode, d1.name, decomposed = dec, depth = depth
    tvlct, r0, g0, b0, /get
    d1 = create_struct(d1, 'decomposed', dec, 'depth', long(depth), $
        'r0',r0, 'g0',g0, 'b0',b0)
    !sgraph.d1 = d1


    ; set the new device to sgraph config.
    device, _extra = opt
    !p.font = 1
    !p.thick = mode.thick
    !p.charthick = mode.thick
    !x.thick = mode.thick & !y.thick = mode.thick & !z.thick = mode.thick
    sgtruecolor
    !p.background = bgc0
    !p.color = sgcolor('black')

    
    if dev0 eq devw then begin
        window, mode.sgid, xsize = mode.area[0], ysize = mode.area[1]
        polyfill, [0,1,1,0,0], [0,0,1,1,0], /normal, color = !p.background
    endif
end

fn = shomedir()+'/idl.png'
fn = 0
sgopen, fn, xsize = 400, ysize = 700
sgclose
end
