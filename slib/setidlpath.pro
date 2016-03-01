;+
; Type: procedure.
; Purpose: Use the given path info to set an executable file, which the
;       user can use to start IDL or IDLDE with the proper IDL path settings.
; Parameters:
;   f0, in, string or strarr[n], opt. The IDL path info, see below.
;       It can be either a list of IDL path info entries or a file that 
;       contains the list. If f0 is not a valid file on disk, it will be
;       viewed as entry list.
; Keywords:
;   outfile, in, string, opt. The output executable file's file name.
;       Default file name idl.exe for Windows,idl.sh for Mac/Unix/Linux.
;       If idlde is set, then the file name will be idlde.exe or idlde.sh.
;       If is a full file name that contains a valid path, then the path will
;       be used, otherwise use home directory. The outfile should not contain
;       extension, in order to be operating system independent.
;   idlexe, in, string, opt. The executable file for idl/idlde. For example
;       /usr/local/exelis/idl/idl83/bin/idlde. Set it to explicitly use
;       certain version of IDL, or when your IDL is not installed at the
;       default location. Setting this keyword will overwrite
;       idlde, idlversion, idlexe.
;   idlde, in, boolean, opt. Set this keyword to use IDLDE instead of IDL.
;   idlversion, in, string, opt. Specify the IDL version to be used. Set
;       it when you have multiple versions of IDL installed. Default is
;       to use the highest version. The format to set this keyword is
;       'idl82','idl71', etc, case insensitive.
;   idldir, in, string, opt. The directory IDL is installed. For example,
;       /usr/local/exelis/ in Linux, /Applications/exelis/ in Mac, 
;       C:\Program Files\Exelis\ in Windows. Set it when your IDL is not
;       installed at the default location.
; Notes: The format rules of the path info entry are.
;       0,  Each line contains one entry of directory. For example, 
;               /Users/sheng/Dropbox/slib/
;           Either / or \ can be used as path separator, but / is preferred.
;           The trailing / or \ is optional.
;       1,  You can include subdirectories recursively using '+'. For example
;               +/Users/sheng/Dropbox/slib/
;       2,  Relative path is allowed. There are 3 ways to do it.
;           a)  Refer to certain directory. For example, here's how to include
;               both path/ and math/ folders in the root directory
;                   $R = /Users/sheng/Dropbox/slib/
;                   $R/path/
;                   $R/math/
;               This is a way of saving typing and keep paths clear.
;
;               You can define multiply root directories, naming rule is simply
;               $<name>, as long as the names are not duplicated, for example
;                   $R1 = /Users/sheng/Dropbox/slib/
;                   $R2 = /Users/sheng/sdt/
;           b)  Refer to your home directory, using '~' (tilde). For example,
;               The above example can be simplified to be
;                   $R1 = ~/Dropbox/slib/
;                   $R2 = ~/sdt/
;           c)  Refer to the current directory, using '.' (dot). In the example
;               in the end, I actually wrote
;                   +.
;               Here the current directory is understood as the directory which
;               contains the procedure setidlpath.pro. Since setidlpath.pro is
;               located in slib/, $R will be where the user saves slib, say
;               C:\Projects\slib\.
;
;               However, when f0 is a file contains the path info, the current
;               directory will be the directory of that file. For example, if
;               f0 is at ~/codes/sdt_path, then the current directory is
;               ~/codes/. This design allows an easy switch between different
;               libraries, simply put f0 with the desired path info under
;               each library's root directory.
;       3,  Output executable file can be set as an entry in the format of
;                   $OUTFILE = ~/idl
; Dependence: none.
; History:
;   2016-02-26, Sheng Tian, create.
;-

pro setidlpath, f0, outfile = ofn0, idlexe = idlexe, $
    idlde = idlde, idlversion = ver0, idldir = idldir

    
; **** find idl executable file's full file name.

    if n_elements(idlexe) eq 0 then begin
        
        ; directory where IDL is installed.
        if n_elements(idldir) eq 0 then begin
            case !version.os of

                'darwin': dirs = '/Applications'

                'Win32': dirs = ['C:/Program Files','C:/Program Files (x86)']
                else: dirs = '/usr/local'
            endcase
            dirs = dirs+'/'
            for i = 0, n_elements(dirs)-1 do begin
                tdirs = dirs[i]+['exelis','itt']
                for j = 0, n_elements(tdirs)-1 do begin
                    if ~file_test(tdirs[i],/directory) then continue
                    idldir = tdirs[i]
                endfor

            endfor
            if n_elements(idldir) eq 0 then idldir = ''
            if file_test(idldir,/directory) eq 0 then $
                read, idldir, prompt = 'Please enter IDL installed directory: '

        endif
        idldir = strjoin(strsplit(idldir,'\',/extract),'/')
        if strmid(idldir,0,1,/reverse_offset) ne '/' then idldir+= '/'
        
        ; idl version, eg 'idl82'.
        ver = (keyword_set(ver0) eq 1)? strlowcase(ver0): 'idl??'

        ; idl executable file.
        dirs = file_search(idldir+ver, count = cnt)
        if cnt eq 0 then begin
            read, idlexe, prompt = 'Please enter IDL execuatalbe path: '
        endif else begin
            case !version.os of
                'Win32': idlexe = dirs[cnt-1]+'/bin/bin.x86_64/'
                else: idlexe = dirs[cnt-1]+'/bin/'
            endcase
            tmp = (keyword_set(idlde) eq 1)? 'idlde': 'idl'
            idlexe = idlexe+tmp
        endelse
    endif
    
    idlexe = strjoin(strsplit(idlexe,'\',/extract),'/')
    if !version.os eq 'Win32' then idlexe = idlexe+'.exe'

    if file_test(idlexe) eq 0 then $
        read, idlexe, prompt = 'Please enter IDL execuatalbe path: '


; **** get the list of IDL path info entries.

    if n_elements(f0) eq 0 then f0 = ['+.']
    if file_test(f0[0]) eq 0 then begin     ; list of entries.
        infos = f0
    endif else begin                        ; a file contains entry list.
        ninfo = file_infos(f0)
        infos = strarr(ninfo)
        openr, lun, f0, /get_lun
        readf, lun, infos
        free_lun, lun
    endelse
    ninfo = n_elements(infos)


; **** parse the entry list.

    ; extract output file.
    pos = stregex(infos,'\$OUTFILE')
    idx = where(pos ne -1, cnt)
    if cnt ne 0 then begin
        tinfo = infos[idx[0]]
        ofn0 = strtrim(strmid(tinfo,strpos(tinfo,'=')+1),2)
    endif
    idx = where(pos eq -1, ninfo)
    infos = infos[idx]


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


; **** assemble the command.

    sep1 = path_sep()
    for i = 0, ninfo-1 do $
        infos[i] = strjoin(strsplit(infos[i],'\',/extract),sep1)

    sep2 = path_sep(/search_path)
    info = strjoin(infos,sep2)

    cmd = '"'+idlexe+'" -IDL_PATH "'+info+'"'
    
    ext = (!version.os_family eq 'unix')? 'sh': 'bat'
    ofn = (keyword_set(idlde) eq 1)? 'myidlde': 'myidl'
    if file_test(f0[0]) eq 1 then begin
        ofn = file_basename(f0[0])
        pos = strpos(f0[0],'.',/reverse_search)
        if pos ne -1 then ofn = strmid(f0[0],0,pos)
    endif
    if n_elements(ofn0) ne 0 then begin
        ofn = ofn0
        dir = file_dirname(ofn0)
    endif else dir = home
    if file_test(dir,/directory) eq 0 then dir = home
    ofn = dir+'/'+ofn+'.'+ext
    ofn = strjoin(strsplit(ofn,'\',/extract),sep1)


    openw, lun, ofn, /get_lun
    printf, lun, cmd
    free_lun, lun

    print, 'The following command is generated:'
    print, cmd
    print, 'Command exported to '+ofn+' ...'

end

f0 = ['+.','$OUTFILE=myidl']
;f0 = ['$R = hahah','+$R/hehe']
setidlpath, f0
end