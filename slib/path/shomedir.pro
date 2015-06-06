;+
; Type: function.
;
; Purpose: Return absolute path of home directory for current user.
;
; Parameters: none.
;
; Keywords:
;   trailing_slash, in, boolean, optional. Set to add a trailing slash.
;
; Return:
;   string. The absolute path of home directory.
;
; Notes: none.
;
; Dependence: none.
;
; History:
;   2012-07-30, Sheng Tian, create.
;-

function shomedir, trailing_slash = trailing_slash

    sep = path_sep()
    
    case !version.os_family of
        'unix'      : homedir = getenv('HOME')
        'Windows'   : homedir = getenv('UserProfile')
        else        : message, 'unknown OS ...'
    endcase
    
    if keyword_set(trailing_slash) then homedir+= sep
    return, homedir
    
end
