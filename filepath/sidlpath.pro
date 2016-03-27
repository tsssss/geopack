;+
; Type: procedure.
;
; Purpose: Interface for IDL path management. Can reset to default; set to
;   specified paths; add in front of existing paths; append to existing paths.
;   Specify paths by array or a contatenated string or a formatted file.
;
; Parameters:
;   info, in, string or strarr[n], optional. The paths to be modified. May 
;       contain an array of paths, or a string with paths separated by : or ;.
;       If not set, try to read paths from filename.
;
; Keywords:
;   filename, in, string, optional. Specify a full filename that contain paths.
;       Will be omitted if path is set.
;   rootdir, in/out, string, optional. Root directory for '*' notation. If set,
;       should be full path or relative to home directory, ie, '~/<path>'.
;       If omitted, find it through filename if available, then if the given
;       path is presented and contain 1 element, use it. Otherwise, try the 
;       root directory of the routine that calls sidlpath. Finally, test the
;       existence of rootdir, if not, set it to current directory.
;   reset, in, boolean, optional. Reset IDL path to original status. Used 
;       separately because setting it will omit all other functionalities.
;   current, in, boolean, optional. Save current !path to user preference file.
;   add, in, boolean, optional. Set to add to existing paths (in front of).
;   append, in, boolean, optional. Set to append to existing paths (after).
;
; Notes: IDL path entry can use the following notations. (1) '~' for home
;   directory, eg, '~/codes/idl'; (2) '..' for parent directory, eg, '../lib';
;   (3) '*' for root directory; (4) '+' before path entry means to include all
;   subdirectories (that contain *.pro or *.sav files).
;
; Dependence: none.
;
; History:
;   2015-06-08, Sheng Tian, create.
;-

pro sidlpath, infos, filename = file, rootdir = rootdir, $
    reset = reset, add = add, append = append, current = current

    ; reset IDL !path.
    if keyword_set(reset) then begin
        tinfo = expand_path('<IDL_DEFAULT>')
        pref_set, 'IDL_PATH', tinfo, /commit
        return
    endif
    
    if keyword_set(current) then begin
        pref_set, 'IDL_PATH', !path, /commit
        return
    endif

    sep = path_sep()
    sep1 = path_sep(/search_path)
    plun = -1   ; output to console.

    ; prepare info.
    if n_elements(infos) eq 0 then begin
        infos = ''   ; if no file or cannot find it.
        if n_elements(file) ne 0 then begin
            if file_test(file) eq 1 then begin
                ninfo = file_lines(file)
                infos = strarr(ninfo)
                openr, lun, file, /get_lun
                readf, lun, infos
                free_lun, lun
            endif
        endif
    endif

    tinfo = strjoin(infos, sep1)
    infos = strsplit(tinfo, ';:,', /extract)    ; all possible path separators.
    ninfo = n_elements(infos)

    ; deal with $Rs.
    pos = stregex(infos,'=')
    idx = where(pos ne -1, cnt)
    
    ; replace each $R.
    for i = 0, cnt-1 do begin
        ti = idx[i]
        tvar = strtrim(strmid(infos[ti],0,pos[ti]),2)   ; var name.
        root = strtrim(strmid(infos[ti],pos[ti]+1),2)   ; root value.
        vlen = strlen(tvar)
        ; replace the $Rs.
        for j = 0, ninfo-1 do begin
            tinfo = infos[j]
            if j eq ti then continue
            tpos = strpos(tinfo, tvar)
            if tpos eq -1 then continue
            infos[j] = strmid(tinfo,0,tpos)+root+strmid(tinfo,tpos+vlen)
        endfor
    endfor
    
    ; remove those lines define root dirs.
    if cnt ne 0 then begin
        infos = infos[where(pos eq -1)]
        ninfo = n_elements(infos)
    endif


;    ; prepare the paths, replace / to \ for Windows.
;    if !version.os_family eq 'Windows' then $
;        for i = 0, npath-1 do $
;            paths[i] = strjoin(strsplit(paths[i],'/',/extract),sep)


    ; deal with home dir '~'.
    pos = stregex(infos,'\~')
    idx = where(pos ne -1, cnt)
    
    tvar = (!version.os_family eq 'unix')? 'HOME': 'UserProfile'
    home = getenv(tvar)
    for j = 0, cnt-1 do $
        infos[j] = strmid(infos[j],0,pos[j])+home+strmid(infos[j],pos[j]+1)


    ; deal with '.'
    pos = stregex(infos,'\.')
    idx = where(pos ne -1, cnt)

    if file_test(f0[0]) eq 1 then curr = file_dirname(f0[0]) else begin
        calls = scope_traceback(/struct)
        ncall = n_elements(calls)
        curr = file_dirname(calls[ncall-2].filename)
    endelse
    if n_elements(rootdir) ne 0 then curr = rootdir

    for j = 0, cnt-1 do $
        infos[j] = strmid(infos[j],0,pos[j])+curr+strmid(infos[j],pos[j]+1)

    
    ; deal with <IDL_DEFAULT>.
    idx = where(infos eq '<IDL_DEFAULT>', cnt)
    if cnt eq 0 and ~keyword_set(add) and ~keyword_set(append) then $
        infos = [infos,'<IDL_DEFAULT>']     ; don't add for add & append modes.
    idx = where(infos eq '<IDL_DEFAULT>')
    infos[idx] = expand_path('<IDL_DEFAULT>')


    ; deal with '+', expand all sub-directories.
    for i = 0, n_elements(infos)-1 do $
        infos[i] = expand_path(infos[i])    ; include all_dirs keyword?

    
    ; add or append.
    if keyword_set(add) then infos = [infos, strsplit(!path, sep1, /extract)]
    if keyword_set(append) then infos = [strsplit(!path, sep1, /extract), infos]
    ninfo = n_elements(infos)

    
    ; check uniqueness, for duplication, only keep the 1st appearance.
    tmp = reverse(infos)        ; reverse b/c uniq keeps the last duplicate.
    idx = uniq(tmp, sort(tmp))  ; get unique index.
    idx = idx[sort(idx)]        ; sort to recover the ordering.
    tmp = tmp[idx[sort(idx)]]
    infos = reverse(tmp)        ; we want to keep the first duplicate.

    ; check existence.
    infos = strsplit(strjoin(infos,sep1),sep1,/extract)
    ninfo = n_elements(infos)
    flags = bytarr(ninfo)
    for i = 0, ninfo-1 do flags[i] = file_test(infos[i],/directory)
    idx = where(flags eq 1b, cnt)   ; cnt must be >0, from <IDL_DEFAULT>.
    infos = infos[idx]

    ; commit the changes to IDL path.
    tinfo = strjoin(infos,sep1)
    pref_set, 'IDL_PATH', tinfo, /commit
    printf, plun, 'IDL path is updated ...'

end
