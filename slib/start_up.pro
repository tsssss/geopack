;+
; Type: procedure.
;
; Purpose: Precisely initiate IDL !path.
;
; Parameters:
;   infile, in, string, optional. The full filename of a formatted file, 
;       which specifies the paths will be added to IDL !path. The format is 
;       explained in 'Notes'. The default infile is 'idl_path' in the same
;       directory with this code, ie, start_up.pro.
;
; Keywords:
;   rootdir, in/out, string, optional. The root directory for relative path.
;       Default is the directory contains this code, ie, start_up.pro; if
;       infile is set, then is the directory contains infile.
;
; Notes: IDL paths are set in the input file. There are several general rules:
;   each line contains one path; path separator is '/', regardless of OS;
;   <IDL_DEFAULT> is added implicitly if omitted.
;       Both absolute path and relative path are OK. For relative path,
;   notations are: (1) use '~' for home directory, eg '~/codes/'; (2) '..' for
;   parent directory, eg '/../lib/'; (3) '*' for root directory.
;       The order of paths matters. If there are name confliction in codes,
;   put the path contains the wanted code before other competing paths.
;       '+' before paths means to include all subdirectories to IDL !path.
;
; Dependence: none.
;
; History:
;   2012-05-17, Sheng Tian, create.
;   2013-10-03, Sheng Tian, use scope_traceback.
;-

pro start_up, infile, rootdir = rootdir

    sep = '/'   ; path separator, works fine with Win, Linux/Unix, Mac OS.

    ; need to check rootdir before checking infile.
    if ~keyword_set(rootdir) then rootdir = ''
    if n_elements(infile) ne 0 then begin
        ifn = file_basename(infile)
        rootdir = file_dirname(infile)
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
    if ifn eq '' then return    ; cannot find infile.
    print, 'reading input file: ', ifn

    ; read paths.
    npath = file_lines(ifn)
    if npath le 0 then return       ; no path.
    paths = strarr(npath) & tpath = ''
    openr, lun, ifn, /get_lun
    for i = 0, npath-1 do begin
        readf, lun, tpath
        paths[i] = tpath
    endfor
    free_lun, lun

    ; deal with "..".
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

    ; deal with "*".
    for i = 0, npath-1 do begin
        pos = strpos(paths[i],'*')
        if pos lt 0 then continue
        paths[i] = strmid(paths[i],0,pos)+rootdir+strmid(paths[i],pos+1)
    endfor

    ; deal with "~".
    sysvar = (!version.os_family eq 'Windows')? 'UserProfile': 'HOME'
    homedir = getenv(sysvar)
    for i = 0, npath-1 do begin
        pos = strpos(paths[i],'~')
        if pos lt 0 then continue
        paths[i] = strmid(paths[i],0,pos)+homedir+strmid(paths[i],pos+1)
    endfor

    ; check <IDL_DEFUALT>.
    idx = where(paths eq '<IDL_DEFAULT>', cnt)
    if cnt eq 0 then paths = [paths,'<IDL_DEFAULT>']
    idx = where(paths eq '<IDL_DEFAULT>')
    paths[idx] = expand_path('<IDL_DEFAULT>')

    ; deal with "+".
    for i = 0, n_elements(paths)-1 do paths[i] = expand_path(paths[i],/all_dirs)

    ; check path existence.
    sep1 = path_sep(/search_path)
    paths = strsplit(strjoin(paths,sep1),sep1,/extract)
    npath = n_elements(paths)
    flags = bytarr(npath)
    for i = 0, npath-1 do flags[i] = file_test(paths[i],/directory)
    idx = where(flags eq 1B, cnt)   ; cnt must be >0, from <IDL_DEFAULT>.
    paths = strjoin(paths[idx],sep1)
    pref_set, 'IDL_PATH', paths, /commit
    
end
