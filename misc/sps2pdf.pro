;+
; Type: procedure.
; Purpose: Wrap ps2pdf in Unix.
; Parameters: fns, in, string, req. File names end with .ps or .eps.
; Keywords: rm, in, boolean, opt. Set to remove the old file.
; Return: none.
; Notes: Need to install GhostScript for Windows, Unix/Linux/Mac OS should 
;   already have ps2pdf command.
; Dependence: none.
; History:
;   2014-11-10, Sheng Tian, create.
;-

pro sps2pdf, fns, rm = rm
    on_error, 2

    nfn = n_elements(fns)
    sep = path_sep()
    
    for i = 0, nfn-1 do begin
        
        tfn = fns[i]
        if file_test(tfn) eq 0 then continue
        
        dir = file_dirname(tfn)     ; directory.
        tfn = file_basename(tfn)
        idx = strpos(tfn,'.',/reverse_search)   ; find the last '.'.
        ext = strmid(tfn,idx+1)     ; extension, ps, eps, pdf.
        tfn = strmid(tfn,0,idx)     ; filename without extention.
        ifn = dir+sep+tfn+'.'+ext
        ofn = dir+sep+tfn+'.pdf'

        ; construct command.
        cmd = 'ps2pdf '
        ; ps or eps?
        if ext eq 'eps' then cmd+= '-dEPSCrop "'
        cmd+= ifn+'" "'+ofn+'"'

        ; do conversion.
        spawn, cmd, ostr, errmsg

        ; find ps2pdf, add it to path.
        if errmsg ne '' then begin
            ; for Windows.
            if !version.os_family eq 'Windows' then $
                message, 'Please install GhostScript printer to device ...'
            ; for Unix/Linux/Mac OS.
            paths = ['/usr/local/bin','/usr/bin','/bin','/usr/sbin', $
                '/sbin','/opt/X11/bin','/usr/texbin','/opt/local/bin']
            spawn, 'which ps2pdf', ostr
            if ostr ne '' then paths = [file_dirname(ostr),paths]
            for j = 0, n_elements(paths)-1 do begin
                if file_test(paths[j]+'/ps2pdf') eq 0 then continue
                setenv, 'PATH='+getenv('PATH')+':'+paths[j]
                break
            endfor
            ; do conversion.
            spawn, cmd, ostr, errmsg
            
            if errmsg[0] ne '' then begin
                cmd = 'pstopdf '
                cmd+= '"'+ifn+'" "'+ofn+'"'
            
                ; for Unix/Linux/Mac OS.
                paths = ['/usr/local/bin','/usr/bin','/bin','/usr/sbin', $
                    '/sbin','/opt/X11/bin','/usr/texbin','/opt/local/bin']
                spawn, 'which pstopdf', ostr
                if ostr ne '' then paths = [file_dirname(ostr),paths]
                for j = 0, n_elements(paths)-1 do begin
                    if file_test(paths[j]+'/pstopdf') eq 0 then continue
                    setenv, 'PATH='+getenv('PATH')+':'+paths[j]
                    break
                endfor
                ; do conversion.
                spawn, cmd, ostr, errmsg
            endif
            if errmsg[0] ne '' then message, 'cannot find ps2pdf or pstopdf ...'
        endif

        ; remove ps/eps files.
        if keyword_set(rm) then begin
            if !version.os_family eq 'Windows' then cmd = 'del "'+ifn+'"' $
            else cmd = 'rm -f "'+ifn+'"'
            spawn, cmd
        endif
    endfor
end