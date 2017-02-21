;+
; mode 1a: find local file.
; mode 1b: find local file, but check remote file to decide download or not.
; mode 2a: fail to find local file, if found remote file, download directly.
; mode 2b: fail to find local file, check remote index file, find the wanted remote file.
;
; provide basefn, locpath to find the file locally, basefn can be in wildcard.
; provide rempath to enable online search, and check for mtime and fsize.
;-

function sgetfile, basefn, locpath0, rempath0, $
    locidx = locidx0, remidx = remidx0, $
    check_mtime = check_mtime, check_fsize = check_fsize, $
    force_download = force_download, use_local_index = use_local_index

    ; check local directory.
    if n_elements(locpath0) eq 0 then message, 'no local path ...'
    locpath = strimpath(locpath0)
    
    ; find local file.
    flag1 = 1   ; 0: local file doesn't exit.
    ext = strmid(basefn, strpos(basefn,'.',/reverse_search))
    lines = file_search(locpath+'/*'+ext)
    idx = where(stregex(lines, basefn) ne -1, cnt)
    if cnt ne 0 then locfn = (reverse(lines[idx]))[0] else flag1 = 0
    
    if flag1 eq 0 then locfn = ''
    
    ; check remote directory, if no remote directory, then do not go online.
    if n_elements(rempath0) eq 0 then return, locfn
    
    ; if local file is found and no check_* keyword is set, then finish.
    if ~keyword_set(check_mtime) and ~keyword_set(check_fsize) then $
        if flag1 eq 1 and ~keyword_set(force_download) then return, locfn
        
    ; prepare index files.
    locidx = (n_elements(locidx0) eq 0)? '.remote-index.html': locidx0
    remidx = (n_elements(remidx0) eq 0)? '': remidx0
    
    ; use remote directory to expand functionality.
    rempath = strimpath(rempath0)
    
    ; find the remote file directly.
    remfn = rempath+'/'+basefn
    reminfo = surlinfo(remfn)
    
    ; check if hit the remote file.
    flag = 0    ; 1 for not hit.
    if size(reminfo,/type) ne 8 then flag = 1 else begin
        if reminfo.size eq 0 then flag = 1
        if reminfo.mtime eq 0 then flag = 1
    endelse
    
    ; find the remote file using remote folder or index file.
    if flag eq 1 then begin
        
        remidx = rempath+'/'+remidx     ; also applies to remote directory.
        locidx = locpath+'/'+locidx

        ; check if need to download remote index.
        flag1 = 1   ; 1 for download remote index.
        if keyword_set(use_local_index) then flag1 = 0
;        if file_test(locidx) eq 0 then flag1 = 1 else begin
;            idxinfo = file_info(locidx)
;            ; do not update locidx if it's older than 1 month.
;            if systime(1)-idxinfo.mtime ge 2.6e6 then flag1 = 0 else flag1 = 1
;        endelse
        if flag1 then begin
            idxinfo = surlinfo(remidx)
            if size(idxinfo,/type) ne 8 then $
                return, '' ; no remote directory or index.
            
            ; download remote folder or index file to local directory.
;            if file_test(locidx) eq 1 then file_delete, locidx
;            s_curl, remidx, locidx, idxinfo
            scurl, remidx, locidx
        endif
        
        ; read local index file.
        nline = file_lines(locidx)
        if nline eq 0 then return, ''
        lines = strarr(nline)
        openr, lun, locidx, /get_lun
        readf, lun, lines
        free_lun, lun
        
        ; check and find the file, with highest version.
        flag2 = 1   ; 1: for remote file exists.
        tmp = stregex(lines, basefn)
        idx = where(tmp ne -1, cnt)
        if cnt eq 0 then flag2 = 0 else begin
            flag2 = 1
            lines = lines[idx]
            tline = (lines[sort(lines)])[cnt-1] ; want the highest version file.
            tline = strmid(tline, stregex(tline, basefn))
            tfn = strmid(tline, 0, stregex(tline, ext))+ext
            remfn = rempath+'/'+tfn
            locfn = locpath+'/'+tfn
        endelse
        
        reminfo = surlinfo(remfn)
    endif
    
    ; no remote file.
    if flag2 eq 0 then return, locfn
    
    ; decide download or not.
    download = 1
    
    if flag1 eq 1 then begin
        ; compare file size or mtime.
        locinfo = file_info(locfn)
        if keyword_set(check_mtime) then $
            if locinfo.mtime ge reminfo.mtime then download = 0
        if keyword_set(check_fsize) then $
            if locinfo.size ge reminfo.size then download = 0
    endif
    if keyword_set(force_download) then download = 1
    
    if download eq 0 then return, locfn
        
    ; download data.
;    s_curl, remfn, locfn, reminfo
    scurl, remfn, locfn
    return, locfn
end

basefn = 'rbspb_l1_esvy_20160119_v[0-9]{2}.cdf'
locdir = shomedir()
remdir = 'http://themis.ssl.berkeley.edu/data/rbsp/rbspb/l1/esvy/2016/'
; fn = sgetfile(basefn, locdir)
fn = sgetfile(basefn, locdir, remdir, locidx = 'idx.html', /force_download)
end

