;+
; Type: function.
; Purpose: Convert time in given format to epoch in UTC (?).
; Parameters:
;   t0s, in, [n], req. Time in given format.
;   format, in, string, opt. Specify the format of time. Can be
;       'epoch', epoch time; 'epoch16', high res epoch time; 'unix', unix time; 
;       'sdt', sdt time; 'jd', Julian day; 'mjd', modified Julian day. For them,
;       format is case insensitive.
;           Default format is 'YYYY-MM-DD hh:mm:ss'. Only case sensitive for MM.
;       YYYY or YY for year, MM for month, DD for date, DOY for day of year,
;       hh for hour, mm for minute, ss or ss.f[*] for second, lll for millisec,
;       ccc for microsec, ppp for picosec.
;       Here are some predefined short hand notations:
;           'g' : 'mm/dd/yyyy hh:mi:ss'    's' : 'yyyy-mm-ddThh:mi:ss'
;           'u' : 'yyyy-mm-dd hh:mi:ss'    'tz': 'yyyy-mm-ddThh:mi:ssZ'
; Keywords:
;   epoch16, in, boolean, opt. Set to return epoch16.
; Return: return, out, double/dblarr(n)/dcomplex/dcomplex(n).
; Notes:
;   Use stoepoch(et,'epoch16') to convert epoch to epoch16.
; Dependence: none.
; History:
;   2010-07-06, Sheng Tian, create.
;   2014-03-12, Sheng Tian, add epoch16.
;-

function stoepoch, t0s, format, epoch16 = epoch16, tt2000 = tt2000

    ; check format, case insensitive and scalar.
    fmt = (n_elements(format) eq 0)? 'null': format[0]
    fmt0 = strlowcase(fmt)  ; case insensitive version.
    
    ; 1-element array.
    if n_elements(t0s) eq 1 then ts = t0s[0] else ts = t0s

    ; epoch, trivial.
    if fmt0 eq 'epoch' then return, ts
    
    ; epoch16.
    if fmt0 eq 'epoch16' then return, real_part(ts)*1d3+imaginary(ts)*1d-9
    
    ; tt2000.
    if fmt0 eq 'tt2000' then begin
        cdf_tt2000, t0s[0], yr, mo, dy, hr, mi, sc, mil, mic, nan, /breakdown_epoch
        cdf_epoch, tmp, yr, mo, dy, hr, mi, sc, double(mil)+1d-3*(mic+1d-3*nan), /compute_epoch
        return, (t0s-t0s[0])*1d-6+tmp
    endif

    ; unix time, sec from 1970-01-01 00:00:00.000 UTC.
    ; epoch0 = 62167219200000D
    if fmt0 eq 'unix' then return, 1000D*ts+62167219200000D
        
    ; Julian day, days from 4713BC-01-01 12:00:00.000 UTC.
    if fmt0 eq 'jd' then return, 86400000D*(ts-1721059.5D)
      
    ; modified Julian day, Julian day - 2400000.5.
    if fmt0 eq 'mjd' then return, 86400000D*(ts+678941D)

    ; sdt epoch, sec from 1968-05-24 00:00:00.000 UTC.
    ; epoch0 = 62116502400000D
    if fmt0 eq 'sdt' then return, 1000D*ts+62116502400000D
    
    ; double input without format, assume epoch.
    if size(t0s,/type) eq size(0d,/type) then return, ts

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
    tmp2 = keyword_set(epoch16)? 9:5
    if tmp eq 0 then ets = keyword_set(epoch16)? dcomplex(0): double(0) $
    else ets = make_array(tmp,type=tmp2)
    fmt0 = fmt  ; save original copy.
    yr0 = fix(strmid(systime(),3,4,/reverse_offset))   ; this year.
    for i = 0, n_elements(ets)-1 do begin
        ; start over.
        fmt = fmt0
        ; parse predictively or definitively.
        if fmt ne 'null' then begin
            ; deal with month first, it's case sensitive.
            pos = strpos(fmt,'MM')
            if pos ne -1 then mo = fix(strmid(ts[i],pos,2)) else mo = -1
            ; minute.
            pos = strpos(fmt, 'mm')
            if pos ne -1 then mi = fix(strmid(ts[i],pos,2)) else mi = 0
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
            ; sec.
            pos = strpos(fmt, 'ss')
            if pos ne -1 then sc = fix(strmid(ts[i],pos,2)) else sc = 0
            ; milli sec.
            pos = strpos(fmt, 'lll')
            if pos ne -1 then ms = fix(strmid(ts[i],pos,3)) else ms = 0
            ; micro sec.
            pos = strpos(fmt, 'ccc')
            if pos ne -1 then mc = fix(strmid(ts[i],pos,3)) else mc = 0
            ; pico sec.
            pos = strpos(fmt, 'ppp')
            if pos ne -1 then pc = fix(strmid(ts[i],pos,3)) else pc = 0
            ; fractional sec.
            pos = stregex(fmt, '\.f+',length=len)
            if pos ne -1 then begin
                tmp = string(double('0'+strmid(ts[i],pos,len)),format='(F11.9)')
                ms = fix(strmid(tmp,2,3))
                mc = fix(strmid(tmp,5,3))
                pc = fix(strmid(tmp,8,3))
            endif
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
            if n_elements(fmt) eq 1 then begin
                tmp = dblarr(6)
            endif else begin
                tmp = fix(strsplit(fmt[1],':.',/extract))
                tmp2 = n_elements(tmp)
                if tmp2 lt 6 then tmp = [tmp,dblarr(6-tmp2)] else tmp = tmp[0:5]
            endelse
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
        if keyword_set(epoch16) then begin
            cdf_epoch16, tmp, yr, mo, dy, hr, mi, sc, ms, mc, pc, /compute_epoch
            ets[i] = tmp
        endif else begin
            cdf_epoch, tmp, yr, mo, dy, hr, mi, sc, ms, /compute_epoch
            ets[i] = tmp
        endelse
    endfor
    return, ets
end

print, 'Test for Epoch ...'
fmt = 'yyyy-MM-dd/hh:mm:ss.fff'
print, sfmepoch(stoepoch('2014-05-19/21:06:07.123',fmt),fmt)
cdf_epoch, et, 2014, 05, 19, 21, 06, 07, 123, /compute_epoch
print, sfmepoch(stoepoch('2014-05-19'))
print, sfmepoch(stoepoch('2014-05-19/21'))
print, sfmepoch(stoepoch('2014-05-19/21:06'))
print, sfmepoch(stoepoch('2014-05-19/21:06:07'),'hh:mm:ss')
print, sfmepoch(stoepoch('2014-05-19/21:06:07.123.456.789',/epoch16),$
    'hh:mm:ss.lll.ccc.ppp')

print, 'Test for arrays ...'
et = stoepoch('2014-05-19/21:06:07.123.456.789',/epoch16)
ets = [et,et]
print, sfmepoch(ets,'YY-MM-DD/hh:mm:ss.lll.ccc.ppp')
end