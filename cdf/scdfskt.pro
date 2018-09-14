;+
; Type: procedure.
; Purpose: Print cdf skeleton table to file.
; Parameters:
;   cdf0, in, string/long, req. Input cdf file's full file name or cdfid.
;   skeleton, out, struct, opt. Skeleton info structure.
; Keywords:
;   filename, in, string, opt. Write skeleton info to this file. If not 
;       specified, then print to console.
; Notes: If skt.header.ngatt = 0, then no skt.att.
;   If skt.header.nzvar+skt.header.nrvar = 0, then no skt.var.
; Dependence: slib.
; History:
;   2013-03-06, Sheng Tian, create.
;-
pro scdfskt, cdf0, skeleton, filename = fn
    compile_opt idl2
    on_error, 0

    ; get cdf id.
    if size(cdf0, /type) eq 7 then begin
        if ~file_test(cdf0) then $
            message, 'file ' + cdf0 + ' does not exist ...'
        cdfid = cdf_open(cdf0)
    endif else cdfid = cdf0

    ; header info.
    cdfinq = cdf_inquire(cdfid)
    cdf_doc, cdfid, vsn, rls, cpy, increment = inc
    vsn = string(vsn, rls, inc, format = '(I0,".",I0,".",I0)')
    cdf_control, cdfid, get_filename = fn0, $
        get_format = fmt, get_numattrs = natts
    ngatt = natts[0]
    nvatt = natts[1]
    nrvar = cdfinq.nvars
    nzvar = cdfinq.nzvars

    header = {$
        copyright: cpy, $
        cdfformat: fmt, $
        decoding: cdfinq.decoding, $
        encoding: cdfinq.encoding, $
        filename: fn0, $
        majority: cdfinq.majority, $
        version: vsn, $
        ngatt: ngatt, $
        nvatt: nvatt, $
        nzvar: nzvar, $
        nrvar: nrvar }

    ; attribute.
    natts = cdfinq.natts
    if natts ne 0 then begin
        attnames = strarr(natts)
        scpflags = bytarr(natts)
        ; check scope.
        for ii = 0, natts-1 do begin
            cdf_attinq, cdfid, ii, attname, scope, dummy
            attnames[ii] = attname
            case strmid(scope, 0, 1) of
                'G': scpflags[ii] = 1B          ; gatt.
                'V': scpflags[ii] = 0B          ; vatt.
                else: message, 'invalid scope: '+scope+' ...'
            endcase
        endfor
    endif

    ; global attribute.
    gatts = {}           ; idl 8.0 or higher.
    if ngatt gt 0 then begin
        attids = where(scpflags eq 1B)
        for ii = 0, ngatt-1 do begin
            ; deal with name.
            attname0 = attnames[attids[ii]]
            attname  = idl_validname(attname0, /convert_all)
            ; read entry info.
            cdf_control, cdfid, attribute = attname0, get_attr_info = attinfo
            nentrys = attinfo.numgentries
            maxentrys = attinfo.maxgentry
            ; read entry values.
            if nentrys le 0 then begin
                gatts = create_struct(gatts, attname, $
                    {name: attname0, value: ''})
            endif else begin
                value = []          ; need idl 8.0 or higher.
                for jj = 0, maxentrys do begin
                    if not cdf_attexists(cdfid, attname0, jj) then continue
                    cdf_attget, cdfid, attname0, jj, tvalue
                    value = [value, tvalue]
                endfor
                gatts = create_struct(gatts, attname, $
                    {name: attname0, value: value})
            endelse
        endfor
    endif

    ; variable attribute.
    if nvatt gt 0 then begin
        attids = where(scpflags eq 0B)
        vattnames0 = attnames[attids]
        vattnames  = idl_validname(vattnames0, /convert_all)
    endif

    ; variables.
    nvar = nrvar+nzvar
    var = {}        ; idl 8,0 or higher.

    ; rvariable.
    if nrvar gt 0 then begin
        for jj = 0, nrvar-1 do begin
            ; variable info.
            varinq = cdf_varinq(cdfid, jj)
            cdf_control, cdfid, variable = jj, $
                get_var_info = varinfo
            recvary = varinq.recvar eq 'VARY'
            info = { $
                natts: 0, $
                name: varinq.name, $
                cdftype: varinq.datatype, $
                nelem: varinq.numelem, $
                iszvar: 0, $
                recvary: recvary, $
                maxrec: varinfo.maxrec+1, $
                dims: cdfinq.dim, $
                dimvary: varinq.dimvar }
            vname0 = info.name
            vname  = idl_validname(strtrim(vname0,2), /convert_all)
            if n_elements(var) ne 0 then begin
                idx = where(tag_names(var) eq strupcase(vname), cnt)
                if cnt ne 0 then vname = vname+'x'  ; fix vars: +15V & -15V.
            endif
            ; vatt.
            vatts = {}
            for kk = 0, nvatt-1 do begin
                if ~cdf_attexists(cdfid, vattnames0[kk], vname0) then continue
                cdf_attget, cdfid, vattnames0[kk], vname0, value
                vatts = create_struct(vatts, vattnames[kk], $
                    {name: vattnames0[kk], value: value})
                info.natts++
            endfor
            ; add to info.
            if n_elements(vatts) ne 0 then $
                info = create_struct(info, 'att', vatts)
            var = create_struct(var, vname, info)
        endfor
    endif

    ; zvariable.
    if nzvar gt 0 then begin
        for jj = 0, nzvar-1 do begin
            ; varialbe info.
            varinq = cdf_varinq(cdfid, jj, /zvariable)
            cdf_control, cdfid, variable = jj, $
                get_var_info = varinfo, /zvariable
            recvary = varinq.recvar eq 'VARY'
            info = { $
                natts: 0, $
                name: varinq.name, $
                cdftype: varinq.datatype, $
                nelem: varinq.numelem, $
                iszvar: 1, $
                recvary: recvary, $
                maxrec: varinfo.maxrec+1, $
                dims: varinq.dim, $
                dimvary: varinq.dimvar }
            vname0 = info.name
            vname  = idl_validname(strtrim(vname0,2), /convert_all)
            if n_elements(var) ne 0 then begin
                idx = where(tag_names(var) eq strupcase(vname), cnt)
                if cnt ne 0 then vname = vname+'x'  ; fix vars: +15V & -15V.
            endif
            ; vatt.
            vatts = {}
            for kk = 0, nvatt-1 do begin
                if ~cdf_attexists(cdfid, vattnames0[kk], vname0) then continue
                cdf_attget, cdfid, vattnames0[kk], vname0, value
                vatts = create_struct(vatts, vattnames[kk], $
                    {name: vattnames0[kk], value: value})
                info.natts++
            endfor
            ; add to info.
            if n_elements(vatts) ne 0 then $
                info = create_struct(info, 'att', vatts)
            var = create_struct(var, vname, info)
        endfor
    endif

    ; free cdfid.
    if size(cdf0, /type) eq 7 then cdf_close, cdfid

    skeleton = { $
        name: header.filename, $
        header: header }

    if n_elements(gatts) ne 0 then $
        skeleton = create_struct(skeleton, 'att', gatts)

    if n_elements(var) ne 0 then $
        skeleton = create_struct(skeleton, 'var', var)

    ; output.
    if n_params() eq 1 then scdfsktlpr, skeleton, fn

end
