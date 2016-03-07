;+
; Type: function.
; Purpose: Convert time in given format to utc.
; Parameters:
;   t0s, in, [n], req. Time in given format.
;   format, in, string, opt. Specify the format of time. Can be
;       'epoch', epoch time; 'epoch16', high res epoch time; 'unix', unix time; 
;       'sdt', sdt time; 'jd', Julian day; 'mjd', modified Julian day. For them,
;       format is case insensitive.
;           Default format is 'YYYY-MM-DD hh:mm:ss'. Only case sensitive for MM.
;       YYYY or YY for year, MM for month, DD for date, DOY for day of year,
;       hh for hour, mm for minute, ss or ss.s[*] for second, lll for millisec,
;       ccc for microsec, ppp for picosec.
;       Here are some predefined short hand notations:
;           'g' : 'mm/dd/yyyy hh:mi:ss'    's' : 'yyyy-mm-ddThh:mi:ss'
;           'u' : 'yyyy-mm-dd hh:mi:ss'    'tz': 'yyyy-mm-ddThh:mi:ssZ'
; Keywords: none.
; Return: return, out, double/dblarr(n)/dcomplex/dcomplex(n). ut time(s).
; Notes: none.
; Dependence: none.
; History:
;   2014-05-20, Sheng Tian, create.
;-

function stotime, t0s, format

    ; check format, case insensitive and scalar.
    fmt = (n_elements(format) eq 0)? 'null': format[0]
    fmt0 = strlowcase(fmt)  ; case insensitive version.
    
    et4ut0 = 62167219200d
    
    ; 1-element array.
    if n_elements(t0s) eq 1 then ts = t0s[0] else ts = t0s

    ; epoch.
    if fmt0 eq 'epoch' then return, 0.001D*ts-et4ut0

    ; epoch16.
    if fmt0 eq 'epoch16' then return, real_part(ts)+imaginary(ts)*1d-12-et4ut0

    ; unix time, sec from 1970-01-01 00:00:00.000 UTC.
    ; epoch0 = 62167219200000D
    if fmt0 eq 'unix' then return, ts
    ; 1000D*ts+62167219200000D
        
    ; Julian day, days from 4713BC-01-01 12:00:00.000 UTC.
    if fmt0 eq 'jd' then return, 86400d*(ts-1721059.5d)-et4ut0
      
    ; modified Julian day, Julian day - 2400000.5.
    if fmt0 eq 'mjd' then return, 86400d*(ts+678941d)-et4ut0

    ; sdt epoch, sec from 1968-05-24 00:00:00.000 UTC.
    ; epoch0 = 62116502400000D
    if fmt0 eq 'sdt' then return, ts-50716800d

    ; string format, listed are some conventional formats.
    case fmt0 of
        'g' : fmt = 'MM/dd/yyyy hh:mm:ss'
        's' : fmt = 'yyyy-MM-ddThh:mm:ss'
        'u' : fmt = 'yyyy-MM-dd hh:mm:ss'
        'tz': fmt = 'yyyy-MM-ddThh:mm:ssZ'
        ''  : fmt = 'yyyy-MM-dd hh:mm:ss'
        else: ; do nothing
    endcase
    tmp = size(ts,/dimensions)
    if tmp eq 0 then uts = 0d else uts = make_array(tmp,/double)
    fmt0 = fmt  ; save original copy.
    yr0 = fix(strmid(systime(),3,4,/reverse_offset))   ; this year.
    for i = 0, n_elements(uts)-1 do begin
        ; start over.
        fmt = fmt0
        ; parse predictively or definitively.
        if fmt ne 'null' then begin
            ; deal with month first, it's case sensitive.
            pos = strpos(fmt,'MM')
            if pos ne -1 then mo = fix(strmid(ts[i],pos,2)) else mo = -1
            ; now case insensitive is fine.
            fmt = strlowcase(fmt)
            ; year.
            pos = stregex(fmt,'yy+', length=len)
            if pos ne -1 then yr = fix(strmid(ts[i],pos,len)) else yr = yr0
            ; date.
            pos = strpos(fmt, 'dd')
            if pos ne -1 then dy = fix(strmid(ts[i],pos,2)) else dy = -1
            ; day of year.
            pos = strpos(fmt, 'doy')
            if pos ne -1 then doy = fix(strmid(ts[i],pos,3)) else doy = -1
            ; hour.
            pos = strpos(fmt, 'hh')
            if pos ne -1 then hr = fix(strmid(ts[i],pos,2)) else hr = 0
            ; minute.
            pos = strpos(fmt, 'mm')
            if pos ne -1 then mi = fix(strmid(ts[i],pos,2)) else mi = 0
            ; sec.
            pos = strpos(fmt, 'ss')
            if pos ne -1 then sc = fix(strmid(ts[i],pos,2)) else sc = 0
            pos = stregex(fmt, '.s+',length=len)
            if pos ne -1 then begin
                tmp = string(double('0'+strmid(ts[i],pos,len)),format='(F11.9)')
                ms = fix(strmid(tmp,2,3))
                mc = fix(strmid(tmp,5,3))
                pc = fix(strmid(tmp,8,3))
            endif
            ; milli sec.
            pos = strpos(fmt, 'lll')
            if pos ne -1 then ms = fix(strmid(ts[i],pos,3)) else ms = 0
            ; micro sec.
            pos = strpos(fmt, 'ccc')
            if pos ne -1 then mc = fix(strmid(ts[i],pos,3)) else mc = 0
            ; pico sec.
            pos = strpos(fmt, 'ppp')
            if pos ne -1 then pc = fix(strmid(ts[i],pos,3)) else pc = 0
        endif else begin
            ; '/' and ' ' and tab separate date & time.
            fmt = strsplit(strcompress(ts[i]),'/ ',/extract)
            ; '-' divide year & date.
            tmp = fix(strsplit(fmt[0],'-',/extract))
            case n_elements(tmp) of
                2: begin yr = tmp[0] & doy = tmp[1] & end
                3: begin yr = tmp[0] & mo = tmp[1] & dy = tmp[2] & end
            endcase
            ; ':' and '.' separate hour, minute, second, etc.
            tmp = fix(strsplit(fmt[1],':.',/extract))
            tmp2 = n_elements(tmp)
            if tmp2 lt 6 then tmp = [tmp,dblarr(6-tmp2)] else tmp = tmp[0:5]
            hr = tmp[0] & mi = tmp[1] & sc = tmp[2]
            ms = tmp[3] & mc = tmp[4] & pc = tmp[5]           
        endelse
        ; round to nearest year.
        if yr le 1000 then begin
            tmp = 10^(floor(alog10(yr))+1)
            yr = ((yr0-yr)/tmp)*tmp+yr
            if yr0-yr ge tmp*0.5 then yr = yr+tmp
        endif
        ; decompose day of year.
        if n_elements(doy) ne 0 then begin
            if doy ne -1 then begin
                tmp = sfmdoy(yr,doy)
                mo = tmp[0] & dy = tmp[1]
            endif
        endif
        
        ; finally.
        cdf_epoch, tmp, yr, mo, dy, hr, mi, /compute_epoch
        uts[i] = 0.001D*tmp-et4ut0+sc+1d-3*(ms+1d-3*(mc+1d-3*pc))
    endfor
    return, uts
end

print, 'Test for Epoch ...'
cdf_epoch, et, 2014, 05, 19, 21, 06, 07, 123, /compute_epoch
;print, sfmepoch(stoepoch('2014-05-19/21:06'))
;print, sfmepoch(stoepoch('2014-05-19/21:06:07'),'hh:mm:ss')
print, sfmepoch(stoepoch('2014-05-19/21:06:07.123.456.789',/epoch16),$
    'hh:mm:ss.lll.ccc.ppp')
end