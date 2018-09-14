;+
; Type: function.
;
; Purpose: Return the parent directory of the code that calls me.
;
; Parameters: none.
;
; Keywords:
;   trailing_slash, in, boolean, optional. Set to add a trailing slash.
;
; Return:
;   string. The aboslute path of the parent directory.
;
; Notes: Returns '.' if called from IDL console, ie, IDL> print, srootdir().
;
; Dependence: none.
;
; History:
;   2012-05-18, Sheng Tian, create.
;   2012-10-03, Sheng Tian, use scope_traceback().
;-

function srootdir, trailing_slash = trailing_slash

    sep = path_sep()

    calls = scope_traceback(/struct)
    ncall = n_elements(calls)
    
    if ncall lt 2 then rootdir = !dir $     ; calls from $MAIN.
    else rootdir = file_dirname(calls[ncall-2].filename)
    
    if keyword_set(trailing_slash) then rootdir+= sep
    return, rootdir

end