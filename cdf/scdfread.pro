;+
; Type: function.
; Purpose: Read data in given cdf file.
; Parameters:
;   cdf0, in, string/long, req. Full file name or cdfid.
;   vnames, in, strarr, opt. Omit to load all vars.
;   recs, in, int/intarr[2], opt. Omit to locate all records for each var.
;       If it is int then it is the record we load, if intarr[2] it is the
;       record range. To load all records, [-1,-1] works too. Negative record
;       is treated as counting from tail, e.g., [2,-2] means [2,nrec-2], -1
;       means to load the last record.
; Keywords:
;   drec, in, int, opt. Record interval to down sample data.
;   skt, in/out, struct, opt. Set to avoid loading it every time.
;   silent, in, boolean, opt. Set to suppress printing messages.
; Return: struct. (Quantities in paranthesis may not be presented.)
;     
;     return = {
;       __name: 'cdf'
;       header: struct
;       att: struct/gatt
;       var: struct_array }
;     
;     header = {
;       __name: 'cdf.header'
;       filename: string
;       cdfformat: string
;       copyright: string
;       decoding: string
;       encoding: string
;       majority: string
;       version: string }
;       
;     gatt = {
;       natts: long
;       attname1: value1
;       attname2: value2
;         ...
;       attnameN: valueN }
;     
;     var = {
;       nvar: long
;       varname1: struct/var1
;       varname2: struct/var2
;         ...
;       varnameN: struct/varN }
;     
;     varN = {
;       name: string
;       value: array of certain type
;       nrecs: long
;       dims: longarr
;       att: struct/att }
; Notes: The data of arrays is in [n, d1, d2, ..., dm], where n is number of 
;   records, [d1, d2, ..., dm] is the dimension of data at each record. For 
;   example, an array of epoch will be in [n]; an array of electric field will 
;   be in [n,3]; an array of aurora image will be in [n,256,256].
; Dependence: slib.
; History:
;   2012-07-31, Sheng Tian, create.
;   2012-10-03, Sheng Tian, documented.
;   2013-03-06, Sheng Tian, revised.
;-
function scdfread, cdf0, vnames, recs0, drec = drec, skt = skt, silent = silent
    
    compile_opt idl2
    on_error, 0
    quiet0 = !quiet
    !quiet = 1
    
    if n_elements(drec) eq 0 then drec = 1

    ; get cdf id.
    if size(cdf0, /type) eq 7 then begin
        if ~file_test(cdf0) then $
            message, 'file ' + cdf0 + ' does not exist ...'
        cdfid = cdf_open(cdf0)
    endif else cdfid = cdf0

    ; read skeleton.
    if n_elements(skt) eq 0 then scdfskt, cdf0, skt

    ; original variable names.
    novar = n_tags(skt.var)
    ovnames = strarr(novar)
    for i = 0, novar-1 do ovnames[i] = skt.var.(i).name

    ; claimed variable names.
    vnames = (n_elements(vnames) ne 0)? vnames: ovnames
    nvar = n_elements(vnames)
    case n_elements(recs0) of
        0: recs = lonarr(nvar,2)-1
        1: recs = [[replicate(recs0,nvar)],[replicate(recs0,nvar)]]
        2: recs = [[replicate(recs0[0],nvar)],[replicate(recs0[1],nvar)]]
    endcase
    
    
    ; read vars.
    vars = []
    if ~keyword_set(silent) then print, 'reading '+skt.name+' ...'
    for i = 0, nvar-1 do begin
        idx = where(strtrim(ovnames,2) eq vnames[i], cnt)
        if cnt eq 0 then continue
        vinfo = skt.var.(idx)
        if ~keyword_set(silent) then print, 'reading '+vinfo.name+' ...'
        ; record info.
        if vinfo.recvary eq 0 then begin
            ; no vary, read one record.
            rec0 = 0 & nrec = 1
        endif else begin
            ; deal with negative recs. [-1,-1] read all, -1 to read last.
            if recs[i,0] eq -1 and recs[i,1] eq -1 then begin
                if n_elements(recs0) eq 1 then recs[i,*]+= vinfo.maxrec $
                else recs[i,*] = [0,vinfo.maxrec]
            endif else begin
                idx = where(recs[i,*] lt 0, cnt)
                if cnt gt 0 then recs[i,idx]+= vinfo.maxrec
            endelse
            rec0 = recs[i,0]>0 & recs[i,1]<=vinfo.maxrec
            nrec = recs[i,1]-rec0
            if nrec le 0 then nrec = 1      ; read one record.
            rec0 <= vinfo.maxrec-1
            if rec0 le 0 then rec0 = 0      ; fix vinfo.maxrec = 0 
            nrec /= drec
        endelse
        if nrec le 0 then begin
            if ~keyword_set(silent) then print, vinfo.name+' no data ...'
            tvar = {name:vinfo.name, value:ptr_new(), nrec:0}
            vars = [vars, tvar]
            continue
        endif
        ; check for nonvary dimensions.
        shrink = total(vinfo.dimvary eq 0) gt 0
        if vinfo.dims[0] eq 0 then shrink = 0  ; scalar element.
        ; read variable.
        if shrink then begin
            cdf_varget, cdfid, vinfo.name, tval, /string, rec_start = 0
            tmp = [nrec,vinfo.dims]
            vals = make_array(type = size(tval,/type), $
                tmp[where([1,vinfo.dimvary] eq 1)])
            for j = 0, nrec-1 do begin
                cdf_varget, cdfid, vinfo.name, tval, /string, $
                    rec_start = rec0+j*drec
                vals[j,*,*,*,*,*,*,*] = srmdim(tval, vinfo.dimvary)
            endfor
        endif else begin
            cdf_varget, cdfid, vinfo.name, vals, /string, $
                rec_start = rec0, rec_interval = drec, rec_count = nrec
            vals = reform(vals)
            ; permute dimensions.
            if nrec ne 1 and size(vals,/n_dimensions) gt 1 then $
                vals = transpose(vals,shift(indgen(n_elements(vinfo.dims)+1),1))
        endelse
        tvar = {name: vinfo.name, value: ptr_new(vals), nrec:long64(nrec)}
        vars = [vars, tvar]
    endfor
    
    ; free cdfid.
    if size(cdf0, /type) eq 7 then cdf_close, cdfid
    
    !quiet = quiet0

    return, vars
end
