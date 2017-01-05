;+
; Type: function.
; Purpose: Calculate panel positions in normal coordinate.
; Parameters:
;   nypanel, in, int, opt. # of panels in y direction. Default is 1.
; Keywords:
;   ypad, in, int, opt. Space b/w y-panels, in # of ycharsize. Default is 0.4.
;   lpad, in, int, opt. Line skip in # of ycharsize. Default is 0.1.
;   position, in, dlbarr[4], opt. In normal coord, set the area for panels.
;   region, in, dblarr[4], opt. In normal coord, set area including margins.
;   lmargin, in, int, opt. # of xcharsize, left margin. Default is 10.
;   rmargin, in, int, opt. # of xcharsize, right margin. Default is 10.
;   tmargin, in, int, opt. # of ycharsize, top margin. Default is 4.
;   bmargin, in, int, opt. # of bcharsize, bottom margin. Default is 4.
;   margins, in, int/intarr[2]/intarr[4], opt. # of charsize. 
;       intarr[4]:[l,b,r,t](in position sense), inarr[2]:[lr,tb], int:[lrtb].
; Return:
;   dblarr[4,n]. The position of panels.
; Notes: In this code, "region" is the area contains margins and plotting area,
;   "positon" is the plotting area that excludes margins, it is also exactly 
;   the boundary the y-panels.
;       Effectively, the difference between region and position is the margins,
;   but region and position are in normal coordinate, margins are in charsize,
;   i.e., we provide 2 ways to set margin. [lrtb]margin overwrite margins, 
;   position overwrites [lrtb]margin.
;       Our panels stack in y-direction, so there is no nxpanel, xpad, etc.
; Dependence: none.
; History:
;   2014-04-08, Sheng Tian, create.
;-
function sgcalcpos, nypanel, nxpanel, ypad = ypad0, xpad = xpad0, $
    position = pos0, region = region, margins = margins, lpad = lpad0, $
    lmargin = lmg0, rmargin = rmg0, bmargin = bmg0, tmargin = tmg0

    ; x- and y-charsize in normal coord.
    xchsz = double(!d.x_ch_size)/double(!d.x_size)
    ychsz = double(!d.y_ch_size)/double(!d.y_size)

    ; # of x- and y-panels.
    nypan = (n_elements(nypanel) eq 0)? 1: nypanel
    nxpan = (n_elements(nxpanel) eq 0)? 1: nxpanel

    ; x- and y-panel skip, line skip in ycharsize.
    ypad = (n_elements(ypad0) eq 0)? 0.4: ypad0
    xpad = (n_elements(xpad0) eq 0)? 4: xpad0   ; b/c title is vertical.
    lpad = (n_elements(lpad0) eq 0)? 0.1: lpad0

    ; margins in ycharsize, mgs in [l,r,t,b].
    case n_elements(margins) of
        ; no settings.
        0: mgs = [10,8,4,6]
        ; [lrtb].
        1: mgs = replicate(margins,4)
        ; [lr],[tb].
        2: mgs = [replicate(margins[0],2),replicate(margins[1],2)]
        ; [l,b,r,t], same as position and region.
        4: mgs = [margins[0],margins[2],margins[3],margins[1]]
        else: message, 'wrong margin ...'
    endcase
    ; overwrite specific margin, still in ycharsize.
    lmg = (n_elements(lmg0) eq 0)? mgs[0]: lmg0
    rmg = (n_elements(rmg0) eq 0)? mgs[1]: rmg0
    tmg = (n_elements(tmg0) eq 0)? mgs[2]: tmg0
    bmg = (n_elements(bmg0) eq 0)? mgs[3]: bmg0

    ; convert margins and skips to normal coord.
    lmg *= ychsz & rmg *= ychsz & tmg *= ychsz & bmg *= ychsz
    ypad *= ychsz & xpad *= ychsz & lpad *= ychsz

    if n_elements(region) eq 0 then region = [0d,0,1,1]
    if n_elements(region) ne 4 then message, 'region must have 4 elements ...'
    x0 = region[0] & y0 = region[1] & x1 = region[2] & y1 = region[3]

    if n_elements(pos0) ne 0 then begin
        if n_elements(pos0) ne 4 then message, 'pos must have 4 elements ...'
        lmg = pos0[0]-region[0] & rmg = region[2]-pos0[2]
        bmg = pos0[1]-region[1] & tmg = region[3]-pos0[3]
    endif

    ; **** consider elements that have size, title, [xy]title, label, etc.
    ; to be implemented
    ; y-direction
    ; title.
    ytls = dblarr(nypan)
    ; xtitle (xticks).
    yxtls = dblarr(nypan)
    
    ; x-direction
    ; ytitle lines.
    xytl = ychsz+lpad
    ; check ztitle lines.
    xztl = 0
    ; check ylabel width.
    xylb = 0

    ; **** calc position ****
    if nypan eq 1 and nxpan eq 1 then $
        return, [x0+lmg, y0+bmg, x1-rmg, y1-tmg]

    ; will shrink unnecessary dimensions in the end.
    pos = dblarr(4,nxpan,nypan)

    ; calc inter panel size.
    if nxpan le 1 then xpads = 0 else xpads = dblarr(nxpan-1)+xpad
    if nypan le 1 then ypads = 0 else ypads = dblarr(nypan-1)+ypad


    ; size of panel's [xy]size. 
    xpan = (x1-x0-lmg-rmg-total(xpads))/nxpan
    ypan = (y1-y0-bmg-tmg-total(ypads))/nypan
    pos[0,*,*] = x0+lmg
    pos[2,*,*] = x1-rmg
    
    ; panel size.
    xpans = dblarr(nxpan)+xpan
    ypans = dblarr(nypan)+ypan
    
    ; upper left position.
;    pos[0,0,*] = lmg
    for i = 0, nxpan-2 do pos[0,i+1,*] = pos[0,i,*]+xpads[i]+xpans[i]
    pos[3,*,0] = y1-tmg
    for i = 0, nypan-2 do pos[3,*,i+1] = pos[3,*,i]-ypads[i]-ypans[i]
    ; lower right position.
    for i = 0, nxpan-1 do pos[2,i,*] = pos[0,i,*]+xpans[i]
    for i = 0, nypan-1 do pos[1,*,i] = pos[3,*,i]-ypans[i]

    return, reform(pos)
end
