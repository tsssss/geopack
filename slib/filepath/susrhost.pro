;+
; Type: function.
; Purpose: Return host and user name.
; Parameters: none.
; Keywords:
;   usr, in, boolean, opt. Set to return username only.
;   host, in, boolean, opt. Set to return hostname only.
; Return: string. Info in usr@host format.
; Notes: none.
; Dependence: none.
; History:
;   2014-05-22, Sheng Tian, create.
;-
function susrhost, usr = usr, host = host
    ; get host.
    spawn, 'hostname', host0
    if keyword_set(host) then return, host0
    ; get usr.
    case !version.os_family of
        'unix': usr0 = getenv('USER')
        'Windows': usr0 = getenv('username')
        else: message, 'unknown os family ...'
    endcase
    if keyword_set(usr) then return, usr0

    return, usr0+'@'+host0
end