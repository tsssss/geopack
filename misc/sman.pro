;+
; Type: procedure.
; Purpose: Print manual between ';+' and ';-'.
; Parameters:
;   module, in, string, req. Module name. 
;       It's assumed that each file contain one module,
;       and the module's name is the filename.
; Keywords: none.
; Notes:
;   * Adopted from man.pro by John Dombeck.
; Dependence: none.
; History:
;   2012-07-02, Sheng Tian, create.
;-

pro sman, module, extension

    ; check extension.
    if n_elements(extention) eq 0 then extension = '.pro'

    ; turn off the compile message.
    oldquiet = !quiet
    !quiet = 1

    ; find the module.
    fn0 = strlowcase(module)+extension
    paths = ['',strsplit(!path, path_sep(/search_path), /extract)]
    npath = n_elements(paths)
    sep = path_sep()
    for ii = 0, npath-1 do begin
        fn = paths[ii]+sep+fn0
        if file_test(fn) then break
    endfor

    if not file_test(fn) then message, 'file not found ...'

    ; read manual.
    nline = file_lines(fn)
    tline = ''
    tenter = string(13B)

    openr, lun, fn, /get_lun

    ; find the head.
    for ii = 0, nline-1 do begin
        readf, lun, tline
        if strmid(tline, 0, 2) ne ';+' then continue
        break
    endfor
    if ii eq nline-1 then return    ; no manual.

    ; find the tail.
    for jj = ii+1, nline-1 do begin
        readf, lun, tline
        if strmid(tline, 0, 2) eq ';-' then break   ; the tail.
        if strmid(tline, 0, 1) eq ';' then strput, tline, ' ', 0
        print, tline        ; manual.
    endfor

    free_lun, lun

    !quiet = oldquiet

end
