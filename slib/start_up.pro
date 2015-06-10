;+
; Type: procedure.
;
; Purpose: Initiate IDL !path. A simplified version of sidlpath.pro.
;
; Parameters:
;   file, in, string, optional. The full filename of a formatted file, 
;       which contains path entries that will be added to IDL !path. The format
;       of the path entry is explained in Notes. The default file is 'idl_path'
;       in the same directory with this code, ie, start_up.pro.
;
; Keywords:
;   rootdir, in/out, string, optional. The root directory for relative path.
;       Default is the directory contains this code, ie, start_up.pro; if
;       file is set, then is the directory contains file.
;
; Notes: IDL path entry can use the following notations. (1) '~' for home
;   directory, eg, '~/codes/idl'; (2) '..' for parent directory, eg, '../lib';
;   (3) '*' for root directory; (4) '+' before path entry means to include all
;   subdirectories (that contain *.pro or *.sav files).
;
; Dependence: none.
;
; History:
;   2012-05-17, Sheng Tian, create.
;-

pro start_up, file, rootdir = rootdir

    sep = path_sep()
    sep1 = path_sep(/search_path)
    plun = -1   ; output to console.

    ; need to check rootdir before checking file.
    if ~keyword_set(rootdir) then rootdir = ''
    if n_elements(file) ne 0 then begin
        ifn = file_basename(file)
        rootdir = file_dirname(file)
    endif else begin
        ifn = 'idl_path'
        rootdir = ''    ; find parent directory.
    endelse

    ; ensure rootdir exists.
    if file_test(rootdir) eq 0 then begin
        call = scope_traceback(/structure)
        rootdir = file_dirname(call[n_elements(call)-1].filename)
    endif

    ifn = rootdir+sep+ifn
    ifn = file_search(ifn, /expand_tilde)  ; ~ not work in Windows.
    ifn = ifn[0]
    if ifn eq '' then return    ; cannot find file.
    print, 'reading input file: ', ifn

    ; read paths.
    npath = file_lines(ifn)
    if npath le 0 then return       ; no path.
    paths = strarr(npath)
    openr, lun, ifn, /get_lun
    readf, lun, paths
    free_lun, lun

    tpath = strjoin(paths, sep1)
    paths = strsplit(tpath, ';:,', /extract)    ; all possible path separators.
    npath = n_elements(paths)

    ; prepare the paths, replace / to \ for Windows.
    if !version.os_family eq 'Windows' then $
        for i = 0, npath-1 do $
            paths[i] = strjoin(strsplit(paths[i],'/',/extract),sep)

    ; deal with homedir, '~'.
    sysvar = (!version.os_family eq 'Windows')? 'UserProfile': 'HOME'
    homedir = getenv(sysvar)
    for i = 0, npath-1 do begin
        pos = strpos(paths[i],'~')
        if pos lt 0 then continue
        paths[i] = strmid(paths[i],0,pos)+homedir+strmid(paths[i],pos+1)
    endfor

    ; deal with rootdir, '*'.
    if file_test(rootdir, /directory) eq 0 then cd, current = rootdir
    for i = 0, npath-1 do begin
        pos = strpos(paths[i],'*')
        if pos lt 0 then continue
        paths[i] = strmid(paths[i],0,pos)+rootdir+strmid(paths[i],pos+1)
    endfor

    ; deal with rootdir's parent directory, '..'.
    for i = 0, npath-1 do begin
        pos = strpos(paths[i],'..')
        if pos lt 0 then continue
        cnt = 0     ; how many levels to go up.
        repeat begin
            cnt+= 1
            paths[i] = strmid(paths[i],0,pos)+strmid(paths[i],pos+3)
            pos = strpos(paths[i],'..')
        endrep until pos lt 0
        tmp = strsplit(rootdir,sep,/extract)
        if strmid(rootdir,0,1) eq sep then tmp = ['',tmp]
        cnt = n_elements(tmp)-1-cnt > 0
        paths[i] = strjoin([tmp[0:cnt],''],sep)+paths[i]
        pos = strpos(paths[i],'+')
        if pos lt 0 then continue
        paths[i] = '+'+strmid(paths[i],0,pos)+strmid(paths[i],pos+1)
        paths[i] = strjoin(strsplit(paths[i],sep,/extract),sep)
    endfor

    ; deal with '+'.
    for i = 0, n_elements(paths)-1 do $
        paths[i] = expand_path(paths[i])    ; include all_dirs keyword?

    ; check <IDL_DEFUALT>.
    idx = where(paths eq '<IDL_DEFAULT>', cnt)
    if cnt eq 0 then paths = [paths,'<IDL_DEFAULT>']
    idx = where(paths eq '<IDL_DEFAULT>')
    paths[idx] = expand_path('<IDL_DEFAULT>')

    ; check uniqueness.
    tmp = reverse(paths)        ; reverse b/c uniq keeps the last duplicate.
    idx = uniq(tmp, sort(tmp))  ; get unique index.
    idx = idx[sort(idx)]        ; sort to recover the ordering.
    tmp = tmp[idx[sort(idx)]]
    paths = reverse(tmp)        ; we want to keep the first duplicate.

    ; check existence.
    paths = strsplit(strjoin(paths,sep1),sep1,/extract)
    npath = n_elements(paths)
    flags = bytarr(npath)
    for i = 0, npath-1 do flags[i] = file_test(paths[i],/directory)
    idx = where(flags eq 1b, cnt)   ; cnt must be >0, from <IDL_DEFAULT>.
    paths = paths[idx]

    ; commit the changes to IDL path.
    tpath = strjoin(paths,sep1)
    pref_set, 'IDL_PATH', tpath, /commit
    printf, plun, 'IDL path is updated ...'

end
