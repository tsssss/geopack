;+
; Type: procedure.
;
; Purpose: Interface for IDL path management.
;
; Parameters:
;   path, in, string or strarr[n], optional. The paths to be modified. May 
;       contain an array of paths, or a string with paths separated by : or ;.
;       If not set, try to read paths from filename.
;
; Keywords:
;   filename, in, string, optional. Specify a full filename that contain paths.
;       Will be omitted if path is set.
;   reset, in, boolean, optional. Reset IDL path to original status. Used 
;       separately because setting it will omit all other functionalities.
;
; Notes: none.
; Dependence: none.
;
; History:
;   2015-06-08, Sheng Tian, create.
;-

pro sidlpath, path, filename = file, rootdir = rootdir, reset = reset

    ; reset IDL !path.
    if keyword_set(reset) then begin
        tpath = expand_path('<IDL_DEFAULT>')
        pref_set, 'IDL_PATH', tpath, /commit
        return
    endif

    sep = path_sep()
    sep1 = path_sep(/search_path)
    plun = -1   ; output to console.

    ; prepare paths.
    if n_elements(path) eq 0 then begin
        path = ''   ; if no file or cannot find it.
        if n_elements(file) ne 0 then begin
            if file_test(file) eq 1 then begin
                nlines = file_lines(file)
                path = strarr(nlines)
                openr, lun, file, /get_lun
                readf, lun, path
                free_lun, lun
            endif
        endif
    endif

    tpath = strjoin(path, sep1)
    paths = strsplit(tpath, sep1, /extract)
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
    if tpath eq '' then begin
        call = scope_traceback(/structure)
        rootdir = file_dirname(call[n_elements(call)-1].filename)
    endif else rootdir = file_dirname(file)
    for i = 0, npath-1 do begin
        pos = strpos(paths[i],'*')
        if pos lt 0 then continue
        paths[i] = strmid(paths[i],0,pos)+rootdir+strmid(paths[i],pos+1)
    endfor

    ; deal with rootdir's parent directory, '..'.
    for i = 0, npath-1 do begin
        pos = strpos(paths[i],'..')
        if pos lt 0 then continue
        cnt = 0
        repeat begin
            cnt+= 1
            paths[i] = strmid(paths[i],0,pos)+strmid(paths[i],pos+3)
            pos = strpos(paths[i],'..')
        endrep until pos lt 0
        tmp = strsplit(rootdir,path_sep(),/extract)
        if strmid(rootdir,0,1) eq path_sep() then tmp = ['',tmp]
        cnt = n_elements(tmp)-1-cnt > 0
        paths[i] = strjoin(tmp[0:cnt],path_sep())+paths[i]
        pos = strpos(paths[i],'+')
        if pos lt 0 then continue
        paths[i] = '+'+strmid(paths[i],0,pos)+strmid(paths[i],pos+1)
    endfor
    
    ; deal with '+'.
    for i = 0, n_elements(paths)-1 do $
        paths[i] = expand_path(paths[i])    ; include all_dirs keyword?

    ; check <IDL_DEFUALT>.
    idx = where(paths eq '<IDL_DEFAULT>', cnt)
    if cnt eq 0 then paths = [paths,'<IDL_DEFAULT>']
    idx = where(paths eq '<IDL_DEFAULT>')
    paths[idx] = expand_path('<IDL_DEFAULT>')
    
    ; check path existence.
    paths = strsplit(strjoin(paths,sep1),sep1,/extract)
    npath = n_elements(paths)
    flags = bytarr(npath)
    for i = 0, npath-1 do flags[i] = file_test(paths[i],/directory)
    idx = where(flags eq 1B, cnt)   ; cnt must be >0, from <IDL_DEFAULT>.
    tpath = strjoin(paths[idx],sep1)

    ; commit the changes to IDL path.
    pref_set, 'IDL_PATH', paths, /commit
    printf, plun, 'IDL path is updated ...'

end
