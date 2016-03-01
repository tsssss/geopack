;+
; Type: <+++>.
; Purpose: <+++>.
; Parameters: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   <+yyyy-mm-dd+>, Sheng Tian, create.
;-

pro setidlpath, f0, idl_dir = idldir, idlde = idlde, outfile = ofn0

; The format rules are:
;
;   0,  Each line contains one entry of directory. For example, 
;
;           /Users/sheng/Dropbox/slib/
;
;       You have two ways to input the path info.
;
;       a)  Create an array of string, each element contains one line of path.
;           Pass the array to setidlpath, as in the example at the end.
;
;       b)  Save the infos into a file, pass the file name to setidlpath.
;
;       Either / or \ can be used as path separator, but / is preferred. The 
;       trailing / or \ is optional.
;
;   1,  You can include subdirectories recursively using '+'. For example
;
;           +/Users/sheng/Dropbox/slib/
;
;       This will include all its subdirectories.
;
;   2,  Relative path is allowed. There are 3 ways to do it.
;
;       a)  Refer to certain directory. For example, I want to include
;           both path/ and math/ folders in /Users/sheng/Dropbox/slib/. I write
;               
;               $R = /Users/sheng/Dropbox/slib/
;               $R/path/
;               $R/math/
;
;           This is a way of saving typing and keep paths clear.
;
;           To set multiply root directories, number them as $R1, $R2, etc.
;
;               $R1 = /Users/sheng/Dropbox/slib/
;               $R2 = /Users/sheng/sdt/
;
;       b)  Refer to your home directory, using '~' (tilde). For example,
;           The above example can be simplified to be
;
;               $R1 = ~/Dropbox/slib/
;               $R2 = ~/sdt/
;
;       c)  Refer to the current directory, using '.' (dot). In the example
;           in the end, I actually wrote
;
;               +.
;
;           In this example, the current directory is understood as the
;           directory contains the procedure setidlpath.pro. Since I put
;           setidlpath.pro under slib, $R1 will be whereever the user saves
;           slib, say C:\Projects\slib\.
;
;           However, if you input a file contains the path info, the current
;           directory will be the directory of that file. For example, if
;           the file passed in is ~/codes/sdt_path, then he current directory
;           is ~/codes/. Therefore, you can easily switch between different
;           libraries, as long as you put a file ontains the path info under
;           each library's root directory.
             
;    if n_elements(idldir) eq 0 then $
;        idldir = (!version.os_family eq 'unix')? '': $
;            'C:/Program Files/Exelis/IDL82/bin/bin.x86_64/'

    ; find where idl is installed.
    if n_elements(idldir) eq 0 then begin
        
        if !version.os_family eq 'unix' then begin
            dirs = ['/usr/local/bin/idl']
            for i = 0, n_elements(dirs)-1 do begin
                spawn, 'ls -l '+dirs[i], msg, errmsg
                if errmsg[0] ne '' then continue
                pos = strpos(msg,'->')
                if pos eq -1 then continue
                idldir = strtrim(strmid(msg,pos+2),2)
            endfor
            idldir = file_dirname(idldir)
        endif else begin
            dirs = ['C:/Program Files/Exelis/', $
                'C:/Program Files (x86)/Exelis/', $
                'C:/Program Files/ITT/', $
                'C:/Program Files (x86)/ITT/']
            for i = 0, n_elements(dirs)-1 do begin
                msg = file_search(dirs[i]+'IDL??', count = cnt)
                if cnt eq 0 then continue
                idldir = msg[cnt-1]+'/bin/bin.x86_64'
                if file_test(idldir+'/idl.exe') eq 0 then idldir = ''
            endfor
        endelse
        
        
    endif
    
    if n_elements(idldir) eq 0 then idldir = ''
    if file_test(idldir,/directory) eq 0 then $
        read, idldir, prompt = 'Please enter the directory of idl: '
    
    tmp = strmid(idldir,0,1,/reverse_offset)
    if tmp ne '/' and tmp ne '\' then idldir+= '/'



    case n_elements(f0) of

        ; no input, set to IDL default path.
        0: infos = ['<IDL_DEFAULT>']

        ; input is a file name contains path info or a string of one line.
        1: begin
            ; if the file exists then read the path info.
            if file_test(f0) eq 1 then begin
                ninfo = file_infos(f0)
                infos = strarr(ninfo)
                openr, lun, f0, /get_lun
                readf, lun, infos
                free_lun, lun
            endif else infos = f0   ; otherwise view it as a path.
        end

        ; input is an array of path info.
        else: infos = f0
    endcase

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

    for j = 0, cnt-1 do $
        infos[j] = strmid(infos[j],0,pos[j])+curr+strmid(infos[j],pos[j]+1)

    
    ; deal with <IDL_DEFAULT>.
    pos = stregex(infos,'<IDL_DEFAULT>')
    idx = where(pos ne -1, cnt)
    if cnt eq 0 then begin
        infos = [infos,'<IDL_DEFAULT>']
        ninfo = n_elements(infos)
    endif


    sep1 = path_sep()
    for i = 0, ninfo-1 do begin
        infos[i] = strjoin(strsplit(infos[i],'/',/extract),sep1)
        infos[i] = strjoin(strsplit(infos[i],'\',/extract),sep1)
    endfor

    sep2 = path_sep(/search_path)
    info = strjoin(infos,sep2)

    cmd0 = (keyword_set(idlde) eq 0)? 'idl': 'idlde'
    cmd = '"'+idldir+cmd0+'" -IDL_PATH "'+info+'"'
    
    ext = (!version.os_family eq 'unix')? 'sh': 'bat'
    ofn = cmd0
    if file_test(f0[0]) eq 1 then begin
        ofn = file_basename(f0[0])
        pos = strpos(f0[0],'.',/reverse_search)
        if pos ne -1 then ofn = strmid(f0[0],0,pos)
    endif
    if keyword_set(ofn0) then ofn = ofn0
    ofn = home+'/'+ofn+'.'+ext
    ofn = strjoin(strsplit(ofn,'\',/extract),sep1)

    openw, lun, ofn, /get_lun
    printf, lun, cmd
    free_lun, lun

    print, 'The following command is generated:'
    print, cmd
    print, 'Command exported to '+ofn+' ...'

end

f0 = ['+.']
;f0 = ['$R = hahah','+$R/hehe']
setidlpath, f0, outfile = 'slib'
end
