
pro s_curl, remfn, locfn, reminfo
    ; used inside other slib routines.
    ; remfn and locfn are separated with '/'.
    ; reminfo contains remote file info, can be omitted.
    
    compile_opt idl2, hidden

    console = -1
    spc = '    '
    mb = 2d^20
    cmd0 = 'Percent'+spc+' Mean '+spc+'  Time  '+spc+'  Time  '
    cmd1 = '  (%)  '+spc+'(MB/s)'+spc+'  Used  '+spc+'  Left  '
    
    
    if n_elements(reminfo) eq 0 then reminfo = surlinfo(remfn)
    if size(reminfo,/type) ne 8 then return     ; no remote directory.
    
    
    ; prepare remote info.
    server = reminfo.server
    port = reminfo.port
    path = reminfo.path
    fsize = reminfo.size
    mtime = reminfo.mtime

    ; prepare local file.
    locpath = file_dirname(locfn)
    if file_test(locpath,/directory) eq 0 then file_mkdir, locpath
    openw, loclun, locfn, /get_lun

    ; connect to server.
    socket, remlun, server, port, /get_lun, error = errmsg
   
    printf, remlun, 'GET '+path+' HTTP/1.0'
    printf, remlun, 'Host: '+server
    printf, remlun, ''
    
    ; read header.
    header = ''
    tline = 'xxx'
    while tline ne '' do begin
        readf, remlun, tline
        header = [header,tline]
    endwhile
    
    
;    ; treat code 301. !!!Does not work.
;    idx = stregex(header, '301 Moved Permanently')
;    idx = where(idx ne -1, cnt)
;    if cnt ne 0 then begin
;        idx = stregex(header, 'Location:')
;        idx = where(idx ne -1)  ; should contain the new location.
;        tline = header[idx[0]]
;        idx = strpos(tline,':')
;        tremfn = strtrim(strmid(tline, idx+1),2)
;        if tremfn ne remfn then s_curl, tremfn, locfn
;    endif
    

    
    ; download file.
    console = -1
    printf, console, 'downloading '+remfn+' ...'
    printf, console, 'file size: '+string(float(fsize)/mb,format='(F0.1)')+' MB'

    if fsize eq 0 then begin    ; download as text file.

        while eof(remlun) eq 0 do begin
            readf, remlun, tline
            printf, loclun, tline
        endwhile
        
    endif else begin            ; download as data file.
        
        ; prepare buffer.
        bufs = 2l^20
        nbuf = (floor(fsize/bufs))[0]
        tbuf = bytarr(bufs,/nozero)
        
        printf, console, cmd0
        printf, console, cmd1
        
        t0 = systime(1)
        t1 = systime(1)
        
        ; download.
        for i = 0, nbuf-1 do begin
            readu, remlun, tbuf
            writeu, loclun, tbuf
            
            if systime(1)-t1 le 5 then continue
            
            t1 = systime(1)
            db = float(i*bufs)
            dt = t1-t0
            cmd = ''
;            cmd+= string(float(fsize)/mb,format='(I-6)')+spc
            cmd+= '  '+string(db/fsize*100,format='(F-4.1)')+' '+spc
            cmd+= string(db/dt/mb,format='(F5.1)')+' '+spc
            cmd+= stodate(dt,'%H:%M:%S')+spc
            cmd+= stodate((float(fsize)/db-1)*dt,'%H:%M:%S')
            printf, console, cmd
        endfor
        
        tbuf = bytarr(fsize mod bufs,/nozero)
        readu, remlun, tbuf
        writeu, loclun, tbuf
        
    endelse
    
    free_lun, remlun
    free_lun, loclun
    
    printf, console, 'saved to '+locfn+' ...'

    ; fix mtime.
    if mtime eq 0 then return   ; no remote mtime.
    stouch, locfn, mtime = mtime

end

remfn = 'http://themis.ssl.berkeley.edu/data/rbsp/rbspb/l1/vb1/2015/'
;remfn = 'http://themis.ssl.berkeley.edu/data/rbsp/rbspa/l1/spec/2013/
remfn = 'https://cdaweb.gsfc.nasa.gov/pub/data/rbsp/rbspa/l2/efw/esvy_despun/2012/rbspa_efw-l2_esvy_despun_20121120_v01.cdf'
remfn = 'http://rbsp.space.umn.edu/data/rbsp/rbspa/l2/esvy_despun/2012/rbspa_efw-l2_esvy_despun_20121120_v01.cdf'
locfn = shomedir()+'/test_download.cdf'
s_curl, remfn, locfn

end