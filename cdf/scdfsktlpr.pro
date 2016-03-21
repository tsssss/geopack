;+
; Type: procedure.
; Purpose: Print given skeleton variable to file or console.
; Parameters:
;   skeleton, in, struct, req. Skeleton variable. See scdfskt.pro.
;   filename, in, string, opt. Output filename, omit to output to console.
; Keywords: none.
; Notes: none.
; Dependence: none.
; History:
;   2013-03-06, Sheng Tian, create.
;   2013-04-04, Sheng Tian, combine scdfskt2file and scdfskt2console.
;-
pro scdfsktlpr, skt, fn
    compile_opt idl2 & on_error, 0

    tenter = ''     ; this works better than string(10B) or string(13B).

    if n_elements(fn) ne 0 then openw, lun, fn, /get_lun $
    else lun = -1       ; console's lun is -1.

    ; =============
    ; print header.
    ; =============
    printf, lun, '! Skeleton table for the "'+skt.name+'" CDF.'
    printf, lun, '! Generated: ' +systime()
    printf, lun, '! CDF version: ' +skt.header.version
    printf, lun, tenter
    printf, lun, '#header'
    printf, lun, tenter
    tformat = '(A, T16, A)'
    printf, lun, 'CDF NAME:', skt.name, format = tformat
    printf, lun, 'DATA ENCODING:', skt.header.encoding, format = tformat
    printf, lun, 'DATA DECODING:', skt.header.decoding, format = tformat
    printf, lun, 'MAJORITY:', skt.header.majority, format = tformat
    printf, lun, 'FORMAT:', skt.header.cdfformat, format = tformat
    printf, lun, tenter
    tformat = '(A, T6, A, T16, A, T26, A, T36, A)'
    printf, lun, '!', 'R.Vars', 'Z.Vars', 'G.Atts', 'V.Atts', format = tformat
    printf, lun, '!', '------', '------', '------', '------', format = tformat
    tformat = '(T6, I0, T16, I0, T26, I0, T36, I0)'
    
    ngatt = skt.header.ngatt
    nvatt = skt.header.nvatt
    nrvar = skt.header.nrvar
    nzvar = skt.header.nzvar
    nvar = nrvar + nzvar
    
    printf, lun, nrvar, nzvar, ngatt, nvatt, format = tformat
    printf, lun, tenter
    
    ; =================
    ; global attribute.
    ; =================
    if ngatt gt 0 then begin
        printf, lun, tenter
        printf, lun, '#GLOBALattributes'
        printf, lun, tenter
        gatt = skt.att
        for ii = 0, ngatt-1 do begin
            tatt = gatt.(ii)
            printf, lun, string(ii, format = '(I4)')+' '+$
                tatt.name, ': ', tatt.value
        endfor
        printf, lun, tenter
    endif

    ; ==========
    ; variables.
    ; ==========
    if nvar gt 0 then begin
        printf, lun, tenter
        printf, lun, '#Variables'
        printf, lun, tenter
        var = skt.var
        for ii = 0, nvar-1 do begin
            tvar = var.(ii)
            recvary = tvar.recvary? 'VARY': 'NOVARY'
            dimvary = strjoin(string(tvar.dimvary,format='(I0)'),', ')
            printf, lun, ii, tvar.name, format = '(I4, " ", A)'
            printf, lun, tenter
            printf, lun, format = '("!", T6, "CDF Type", T20, "# Elem", T30,'+$
                '"Max Rec", T42, "Rec Vary", T55, "Dimensions", T70,"Dim Vary")'
            printf, lun, format = '("!", T6, "--------", T20, "------", T30,'+$
                '"-------", T42, "--------", T55, "----------", T70,"--------")'
            printf, lun, tvar.cdftype, tvar.nelem, tvar.maxrec, recvary, $
                strjoin(string(tvar.dims, format = '(I0)'), ', '), dimvary, $
                format = '(T6, A, T20, I0, T30, I0, T42, A, T55, A, T70, A)'
            printf, lun, tenter
            ; variable attribute.
            if tvar.natts ne 0 then begin
                vatt = tvar.att
                for jj = 0, tvar.natts-1 do begin
                    tatt = vatt.(jj)
                    printf, lun, tatt.name, string(tatt.value), $
                        format = '(T6, A, T40, A)' 
                endfor
                printf, lun, tenter
            endif
        endfor
    endif

    if lun ne -1 then free_lun, lun
end
