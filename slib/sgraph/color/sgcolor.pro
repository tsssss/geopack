;+
; Type: function.
; Purpose: Get 24-bit color or emulate 8-bit color table.
; Parameters:
;   c0, in, strarr/intarr, req. Input color(s), string for true color, 
;       int for index color.
; Keywords:
;   ct, in, int, opt. Color table number, set if color is in int.
;   triplet, in, boolean, opt. Set to return [n,3] color in rrggbb.
;   names, in, boolean, opt. Set to return all colors IDL knows.
; Return: bytarr[n,3]/tripletarr[n]. Color in true color regime, 
;   in [n,3] by default, in [n] if triplet is set.
; Notes: none.
; Dependence: none.
; History:
;   2014-04-05, Sheng Tian, create.
;-

function sgcolor, c0, ct = ct, triplet = triplet, names = names
    common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr
    
    colors = tag_names(!color)
    if keyword_set(names) then return, colors

    ncolor = n_elements(c0)
    if ncolor eq 0 then return, -1
    
    ; [r,g,b].
    rgb = bytarr(ncolor,3)
    if size(c0, /type) eq 7 then begin
        for i = 0, ncolor-1 do begin
            idx = where(strupcase(c0[i]) eq colors, cnt)
            if cnt eq 0 then message, c0+', no such color ...'
            rgb[i,*] = !color.(idx)
        endfor
    endif else begin    ; c0 is in number, interpret as index color.
        if ~keyword_set(ct) then ct = 0
        loadct, ct, rgb_table = rgb0, /silent
        nc = n_elements(rgb0)/3
        rgb = rgb0[c0 mod nc,*]
    endelse

    if keyword_set(triplet) then return, rgb
    rgb = 256L*(256L*rgb[*,2]+rgb[*,1])+rgb[*,0] ; bbggrr.
    if n_elements(rgb) eq 1 then return, rgb[0]
    return, rgb
end

print, sgcolor([6,4,2])
print, sgcolor(['red','green','blue'])
end
