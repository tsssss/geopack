;+
; Type: procedure (internal).
; Purpose: Return norminal filenames, for given pattern and time.
;   Pure string operations, do not check file existence.
; Parameters:
;   tr0, in, double/string or dblarr[2]/strarr[2], optional. If in double or
;       string, set the time; if in dblarr[2] or strarr[2], set the time range.
;       For double or dblarr[2], it's the unix time or UTC. For string or
;       strarr[2], it's the formatted string accepted by stoepoch, e.g.,
;       'YYYY-MM-DD/hh:mm'. Omitting tr0 means no time replacement.
; Keywords:
;   dt, in, double, optional. If tr0 is a time range, dt sets the steps of
;       the wanted times. In micro second, e.g. 3600000 for hour. Default 
;       time step is day.
;   paths, in, strarr[n], required. The paths. Suggested organization is
;       [root directory, paths, basename].
;   flags, in, bytarr[n], optional. The flag to tell which paths need time
;       replacement. 1 for time replacement, 0 for no. By default, no time 
;       replacement if tr0 is undefined. When tr0 is defined, do not do time
;       replacement for the root directory.
; Notes: none.
; Dependence: slib.
; History:
;   2014-04-28, Sheng Tian, create.
;   2016-03-05, Sheng Tian, updated.
;-

function sprepfile, tr0, dt = dt, paths = paths, flags = flags

    ; prepare ets, paths, flags, files.
    npath = n_elements(paths)
    if npath eq 0 then message, 'no input path ...'
    
    ; 1 file, several or no times.
    if n_elements(tr0) eq 0 then begin
        flags = bytarr(npath)
    endif else begin
        ; times.
        if size(tr0,/type) eq 7 then tformat = '' else tformat = 'unix'
        ets = stoepoch(tr0,tformat)
        if n_elements(tr0) eq 2 then ets = sbreaktr(ets, dt)
        ; time flags.
        if n_elements(flags) eq 0 then begin
            flags = bytarr(npath)+1
            if npath gt 1 then flags[0] = 0   ; protect root directory.
        endif
    endelse
    
    ; # of files.
    nfile = n_elements(ets)/2 > 1
    
    files = strarr(nfile)
    
    for i = 0, nfile-1 do begin
        idx = where(flags eq 1, cnt)
        tpaths = paths
        if cnt ne 0 then begin
            tmp = strjoin(tpaths[idx],'%')
            tmp = sfmepoch(ets[0,i],tmp)
            tpaths[idx] = strsplit(tmp,'%',/extract)
        endif
        files[i] = strjoin(tpaths,'/')
    endfor
    
    return, files[uniq(files)]
end

locroot = spreproot('themis')
remroot = 'http://themis.ssl.berkeley.edu/data/themis'
sites = ['atha']
baseptn = 'thg_l1_asc_'+sites+'_yyyyMMdd_v01.cdf'
rempaths = [remroot,'thg','l1','asi',sites,'yyyy','MM',baseptn]
locpaths = [locroot,'thg','l1','asi',sites,'yyyy','MM',baseptn]

tr0 = '2013-04-05'

locfns = sprepfile(tr0, paths = locpaths)
remfns = sprepfile(tr0, paths = rempaths)

end