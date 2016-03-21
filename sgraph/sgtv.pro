;+
; wrap IDL's tv.
; if color mode is true color, using it.
; if color mode is index color, or if image is 2-d array, 
; convert indexed colors to true color.
; for 3-d array image, set true to avoid ambiguity.
; position must be normal coord.
;-
pro sgtv, img0, position = pos0, _extra = extra
    common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr
    
    if n_elements(r_curr) eq 0 then loadct, 0, /silent
    r0 = r_curr & g0 = g_curr & b0 = b_curr
    
    device, get_decomposed = cmode0 ; 0 for index, 1 for true.
    
    imgsz = size(img0,/dimensions)
    case n_elements(imgsz) of
        2: cmode = 0
        3: cmode = 1
        else: message, 'invalid image dimension ...'
    endcase
    
    ; cmode is determined by image size.
    if cmode eq 1 then img = img0 else begin
        img = bytarr([3,imgsz])
        true = 1
        ncolor = n_elements(r_curr)
        if ncolor eq 0 then loadct, 0, /silent
        ncolor = n_elements(r_curr)
        img[0,*,*] = r0[img0 mod ncolor]
        img[1,*,*] = g0[img0 mod ncolor]
        img[2,*,*] = b0[img0 mod ncolor]
    endelse
    
    ; determine true.
    imgsz = size(img,/dimensions)
    if n_elements(true) eq 0 then true = where(imgsz eq 3)+1
    imgsz = imgsz(where(imgsz ne 3))
    
    ; set true color.
    sgtruecolor
    
    ; deal with position.
    if n_elements(pos0) eq 0 then begin
        tv, img, true = true, _extra = extra
    endif else begin
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
        
        if strlowcase(!d.name) eq 'ps' then begin
            loadct, 0, /silent
            tv, img, true = true, x0, y0, xsize = xsz0, ysize = ysz0, _extra = extra
            tvlct, r0, g0, b0
            r_curr = r0 & g_curr = g0 & b_curr = b0
            r_orig = r0 & g_orig = g0 & b_orig = b0
        endif else begin
            tmp = bytarr(3,xsz0,ysz0)
            for i = 0, 2 do $
                tmp[i,*,*] = congrid(reform(img[i,*,*]), xsz0, ysz0)
            tv, tmp, true = true, x0, y0, xsize = xsz0, ysize = ysz0, _extra = extra
        endelse
    endelse

    device, decomposed = cmode0

end

sgtruecolor
sz = [500,300]
img0 = dist(sz[0],sz[1])
;get_data, 'asi', tmp, img0
ct = 43

sgpsopen, '~/test.eps', xsize = 500, ysize = 300
loadct2, ct
erase, color = sgcolor('white')
;polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color = sgcolor('grey')
;plot, findgen(10), color = sgcolor('lime'), /noerase, /nodata
;oplot, findgen(10), color = sgcolor('yellow')
sgtv, img0, position = [0.1,0.1,0.9,0.9]
sgpsclose, /pdf

sgzopen, '~/test.png', xsize = 500, ysize = 300
loadct2, ct
sgtruecolor
erase, color = sgcolor('white')
;plot, findgen(10), color = sgcolor('lime'), /noerase, /nodata
;oplot, findgen(10), color = sgcolor('yellow')
sgtv, img0, position = [0.1,0.1,0.9,0.9]
sgzclose

end
