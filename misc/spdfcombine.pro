;+
; Type: procedure.
; Purpose: Combine pdfs into 1 pdf.
; Parameters:
;   in0, in, string/strarr[n], req. If in string, it's the directory of the
;       pdfs. If in strarr[n], it's the full filenames of the pdfs.
; Keywords:
;   ofn, in, string, opt. Full file name of the combined pdf.
; Notes: Need GhostScript installed, work for Mac/Unix/Linux.
; Dependence: none.
; History:
;   2015-11-04, Sheng Tian, create.
;-

pro spdfcombine, in0, filename = ofn

    if n_elements(in0) eq 0 then message, 'no input ...'
    if n_elements(in0) eq 1 then begin
        dir = in0
        if file_test(dir) eq 0 then $
            message, 'directory does not exist ...'
        fns = file_search(in0+'/*.pdf', count = nfn)
        if nfn eq 0 then return     ; no file.
    endif else begin
        fns = in0
        nfn = n_elements(in0)
    endelse

    if n_elements(ofn) eq 0 then ofn = shomedir()+'/binder.pdf'

    cmd = '/usr/local/bin/gs '
    cmd+= '-dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sPDFSETTING=/prepress -sOutputFile='
    cmd+= '"'+ofn+'"'
    for i = 0, nfn-1 do cmd+= ' "'+fns[i]+'"'

    spawn, cmd
end
