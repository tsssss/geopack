;+
; Type: procedure.
; Purpose: Draw image.
; Parameters:
;   img0, in, intarr[m,n]/intarr[m,n,3]/intarr[m,3,n]/intarr[3,m,n], req.
;       The image. The size of image determines the mode to draw image.
;       If img0 is in [m,n], then emulate index color mode in true color mode.
;       Otherwise img0 is already in true color mode.
; Keywords:
;   position, in, dblarr[4], opt. The position in normalized coord. Default is
;       calculated by sgcalcpos.
;   ct, in, int, opt. The color table id. Omitted if img0 is in true color
;       mode. Default value is 0.
;   file, in, string, opt. Set the file name of color table files, then ct is
;       the color table id in that file.
; Notes: none.
; Dependence: slib.
; History:
;   2016-03-22, Sheng Tian, rewrite.
;-
pro sgtv, img0, position = pos0, ct = ct0, file = file, _extra = extra
    
    ; store current color.
    common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

;    r0 = r_orig & g0 = g_orig & b0 = b_orig
;    r1 = r_curr & g1 = g_curr & b1 = b_curr
    tvlct, rgb0, /get
    
    ; determine index color mode or true color mode.
    imgsz = size(img0,/dimensions)
    case n_elements(imgsz) of
        2: cmode = 0
        3: cmode = 1
        else: message, 'invalid image dimension ...'
    endcase
    
    if cmode eq 1 then begin    ; true color input.
        img = img0
    endif else begin            ; index color input.
        ; emulate true color mode.
        img = bytarr([3,imgsz]) ; [r,g,b], true = 1.
        ; determine and load color table.
        ncolor = 255
        if n_elements(ct0) eq 0 then tvlct, rgb1, /get $
        else rgb1 = sgcolor(bindgen(ncolor), ct = ct, file = file, /triplet)
        
        img[0,*,*] = rgb1[img0 mod ncolor,0]
        img[1,*,*] = rgb1[img0 mod ncolor,1]
        img[2,*,*] = rgb1[img0 mod ncolor,2]
    endelse

    
    ; determine interlace mode.
    imgsz = size(img,/dimensions)
    true = where(imgsz eq 3)
    imgsz = imgsz(where(imgsz ne 3))
    img1 = transpose(img,shift([0,1,2],-true)) ; convert to [3,imgsz].
    true = 1
    
    ; sgtv should operate in true color mode always.
;    sgtruecolor
    
    ; deal with position.
    if n_elements(pos0) eq 0 then pos1 = sgcalcpos()

    ; centering image, if shrink is needed.
    pos1 = convert_coord(pos0[[0,2]],pos0[[1,3]],/normal,/to_device)
    pos1 = pos1[[0,1,3,4]]      ; pos in pixel, i.e., device coord.
    ysz1 = fix(pos1[3]-pos1[1]) ; target ysize in pixel.
    xsz1 = fix(pos1[2]-pos1[0]) ; target xsize in pixel.
    asp1 = double(ysz1)/xsz1
    ysz0 = imgsz[1]             ; image ysize in pixel.
    xsz0 = imgsz[0]             ; image xsize in pixel.
    asp0 = double(ysz0)/xsz0
    if asp1 gt asp0 then begin  ; fit x.
        xsz0 = xsz1
        ysz0 = xsz0*asp0
        x0 = pos1[0]
        y0 = pos1[1]+(ysz1-ysz0)/2
    endif else begin
        ysz0 = ysz1
        xsz0 = ysz0/asp0
        x0 = pos1[0]+(xsz1-xsz0)/2
        y0 = pos1[1]
    endelse 
    
    img = bytarr(3,xsz0,ysz0)
    for i = 0, 2 do $
        img[i,*,*] = congrid(reform(img1[i,*,*]), xsz0, ysz0)
    if strlowcase(!d.name) eq 'ps' then loadct, 0, /silent
    tv, img, true = true, x0, y0, xsize = xsz0, ysize = ysz0, _extra = extra

    ; restore color.
;    r_orig = r0 & g_orig = g0 & b_orig = b0
;    r_curr = r1 & g_curr = g1 & b_curr = b1
    tvlct, rgb0

end

dir = shomedir()
sgtruecolor
sz = [500,300]
img0 = dist(sz[0],sz[1])
img0 = bytscl(img0,min=0,top=255)
ct = 43

sgopen, dir+'/test.pdf', xsize = 500, ysize = 300
erase, color = sgcolor('white')
;polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color = sgcolor('grey')
;plot, findgen(10), color = sgcolor('lime'), /noerase, /nodata
;oplot, findgen(10), color = sgcolor('yellow')
sgtv, img0, position = [0.1,0.1,0.9,0.9], ct = ct
sgclose

sgopen, dir+'/test.png', xsize = 500, ysize = 300
sgtruecolor
erase, color = sgcolor('white')
;plot, findgen(10), color = sgcolor('lime'), /noerase, /nodata
;oplot, findgen(10), color = sgcolor('yellow')
sgtv, img0, position = [0.1,0.1,0.9,0.9], ct = ct, file = 'ct2'
sgclose

end
