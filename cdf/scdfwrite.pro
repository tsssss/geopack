
pro scdfwrite, cdf0, vname, recs, skt = skt, value = val, attributes = attinfo, reset = reset

    compile_opt idl2
    on_error, 0
    quiet0 = !quiet
    !quiet = 1

    ; get cdf id.
    if size(cdf0, /type) eq 7 then begin
        if ~file_test(cdf0) then begin  ; creat the file.
            message, 'file ' + cdf0 + ' does not exist ...', /continue
            cdfid = cdf_create(cdf0)
            newfile = 1
        endif else newfile = 0
        if newfile eq 0 then cdfid = cdf_open(cdf0)
    endif else cdfid = cdf0


    ; read skeleton.
    if n_elements(skt) eq 0 then scdfskt, cdfid, skt

    ; original variable names.
    novar = skt.header.nrvar+skt.header.nzvar
    if novar eq 0 then newvar = 1 else begin
        ovnames = strarr(novar)
        for i = 0, novar-1 do ovnames[i] = skt.var.(i).name
        idx = where(ovnames eq vname, cnt)
        if cnt eq 0 then newvar = 1 else newvar = 0
    endelse
    
    if newvar eq 0 and keyword_set(reset) then begin
        newvar = 1
        for i = 0, novar-1 do if skt.var.(i).name eq vname then break
        vattinfo0 = skt.var.(i).att
        cdf_vardelete, cdfid, vname, zvariable = skt.var.(i).iszvar
    endif
    
    if newvar eq 1 then begin
        if n_elements(dimvary) eq 0 then begin
            varid = cdf_varcreate(cdfid, vname, _extra = extra, /zvariable)
        endif else begin
            varid = cdf_varcreate(cdfid, vname, dimvary, _extra = extra, /zvariable)
        endelse
    endif

    cdf_varput, cdfid, vname, val
    
    ; variable attributions.
    nvatt = skt.header.nvatt
    if nvatt gt 0 then vatts = tag_names(skt.att)
    if n_elements(vattinfo) ne 0 then begin
        vattnames = tag_names(vattinfo)
        nvattname = n_elements(vattnames)
        for i = 0, nvattname-1 do begin
            if nvatt eq 0 then newatt = 1 else begin
                idx = where(vatts eq vattnames[i], cnt)
                if cnt eq 0 then newatt = 1 else newatt = 0
            endelse
            if newatt eq 1 then $
                attid = cdf_attcreate(cdfid, vattnames[i], /variable_scope)
            cdf_attput, cdfid, vattnames[i], vname, vattinfo.(i)
        endfor
    endif
    
    if n_elements(vattinfo0) ne 0 then begin
        vattnames = tag_names(vattinfo0)
        nvattname = n_elements(vattnames)
        for i = 0, nvattname-1 do begin
            cdf_attput, cdfid, vattnames[i], vname, vattinfo0.(i).value
        endfor
    endif
    
    cdf_close, cdfid

end

fn = shomedir()+'/test.cdf'
vname = 'angle'
val = 1
scdfwrite, fn, vname, value = val
cdf = scdfread(fn, vname)
print, *cdf[0].value

val = [2,3,4]
scdfwrite, fn, vname, value = val
cdf = scdfread(fn, vname)
print, *cdf[0].value

file_delete, fn
end