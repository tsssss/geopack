;+
; Type: function.
;
; Purpose: Return absolute path of given hard disk drive name.
;
; Parameters:
;   disk, in, string, required. Name of the hard disk drive.
;
; Keywords:
;   trailing_slash, in, boolean, optional. Set to add a trailing slash.
;
; Return:
;   string. The aboslute path of the hard disk.
;
; Notes: none.
;
; Dependence: none.
;
; History:
;   2012-07-30, Sheng Tian, create.
;   2013-03-15, Sheng Tian, use spawn in Windows.
;-

function sdiskdir, disk, trailing_slash = trailing_slash

    if n_elements(disk) eq 0 then message, 'no hard disk name ...'
    sep = path_sep()    ; OS dependent.
    
    case !version.os of
        'linux': diskdir = strjoin(['','media',disk], sep)
        'darwin': diskdir = strjoin(['','Volumes',disk], sep)
        'Win32': begin
            spawn, 'powershell Get-WmiObject Win32_LogicalDisk', outputs
            idx = where(stregex(outputs, disk) ne -1, cnt)  ; hitted disk info.
            if cnt eq 0 then message, 'no such disk: '+disk+' ...'
            if cnt gt 1 then message, 'more than 1 disk found ...', /continue
            idx1 = idx[0]               ; use the first disk info.
            idx = where(outputs eq '')  ; find disk info blocks.
            tmp = (where(idx ge idx1))[0]
            outputs = outputs[idx[tmp-1]:idx[tmp]]    ; locate the block.
            idx = where(stregex(outputs, 'DeviceID') ne -1)
            diskdir = strmid(strtrim(outputs[idx],2),1,/reverse_offset)
            end
        else: message, 'unkown OS ...'
    endcase
    
    if not file_test(diskdir, /directory) then message, 'no such drive ...'
    if keyword_set(trailing_slash) then diskdir+= sep
    return, diskdir

end