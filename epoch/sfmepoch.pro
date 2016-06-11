;+
; Type: function.
; Purpose: Convert epoch in UTC (?) to various formats.
; Parameters:
;   et0s, in, double/dblarr(n), req. Epoch in UTC (?).
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
;   epoch16, in, boolean, opt. No need to set, just in case.
; Return: return, [*]. Array of specified type.
; Notes:
;   Use sfmepoch(et16,'epoch') to convert epoch16 to epoch.
;       The pattern is recursively substituted. Read once is ok for stoepoch.
; Dependence: none.
; History:
;   2010-07-06, Sheng Tian, create.
;   2014-02-12, Sheng Tian, add epoch16.
;-
function sfmepoch, et0s, format, epoch16 = epoch16, tt2000 = tt2000

    ; check format, case insensitive and scalar.
    fmt = (n_elements(format) eq 0)? '': format[0]
    fmt0 = strlowcase(fmt)  ; case insensitive version.
    fmt1 = fmt              ; cas sensitive version.
    
    ; check epoch16.
    if size(et0s,/type) eq 9 then epoch16 = 1
    if keyword_set(epoch16) then begin
        if fmt0 eq 'epoch16' then return, et0s
        ets = real_part(et0s)*1d3+imaginary(et0s)*1d-9
    endif else if keyword_set(tt2000) then begin
        if fmt0 eq 'tt2000' then return, et0s
        cdf_tt2000, et0s[0], yr, mo, dy, hr, mi, sc, mil, mic, nan, /breakdown_epoch
        cdf_epoch, tmp, yr, mo, dy, hr, mi, sc, double(mil)+1d-3*(mic+1d-3*nan), /compute_epoch
        ets = (et0s-et0s[0])*1d-6+tmp
    endif else ets = et0s
    
    ; 1-element array.
    if n_elements(ets) eq 1 then ets = ets[0]

    ; epoch, trivial.
    if fmt0 eq 'epoch' then return, ets
    
    ; epoch16.
    if fmt0 eq 'epoch16' then begin
        tmp = ets*1d-3 & tmp2 = tmp mod 1
        return, dcomplex(tmp-tmp2, tmp2*1d12)
    endif

    ; unix time, sec from 1970-01-01 00:00:00.000 UTC.
    ; epoch0 = 62167219200000D
    if fmt0 eq 'unix' then return, 0.001D*ets-62167219200D

    ; Julian day, days from 4713BC-01-01 12:00:00.000 UTC.
    if fmt0 eq 'jd' then return, (1D/86400000D)*ets+1721059.5D
      
    ; modified Julian day, Julian day - 2400000.5.
    if fmt0 eq 'mjd' then return, (1D/86400000D)*ets-678941D

    ; sdt epoch, sec from 1968-05-24 00:00:00.000 UTC.
    ; epoch0 = 62116502400000D
    if fmt0 eq 'sdt' then return, 0.001D*ets-62116502400D
    
    ; string input without format, assume string epoch.
    if size(et0s,/type) eq size('',/type) then return, ets

    ; string format, listed are some conventional formats.
    case fmt1 of
        'g' : fmt1 = 'MM/dd/yyyy hh:mm:ss'
        's' : fmt1 = 'yyyy-MM-ddThh:mm:ss'
        'u' : fmt1 = 'yyyy-MM-dd hh:mm:ss'
        'tz': fmt1 = 'yyyy-MM-ddThh:mm:ssZ'
        ''  : fmt1 = 'yyyy-MM-dd/hh:mm:ss'
        else: ; do nothing
    endcase
    tmp = size(ets,/dimensions)
    if tmp eq 0 then etstrs = '' else etstrs = make_array(tmp,/string)
;    fmt0 = fmt  ; save original copy.
    for i = 0, n_elements(ets)-1 do begin
        ; start over.
        fmt = fmt1
        ; break down epoch.
        if keyword_set(epoch16) then begin
            cdf_epoch16, et0s[i], yr, mo, dy, hr, mi, sc, $
                ms, mc, pc, /breakdown_epoch
        endif else begin
            cdf_epoch, ets[i], yr, mo, dy, hr, mi, sc, ms, /breakdown_epoch
            mc = 0 & pc = 0
        endelse
        ; deal with month first, it's case sensitive.
        pos = strpos(fmt,'MM')
        while pos ne -1 do begin
            strput, fmt, string(mo,format='(I02)'), pos
            pos = strpos(fmt,'MM')
        endwhile
        ; now case insensitive is fine.
;        fmt = strlowcase(fmt)
        ; year.
        pos = stregex(fmt,'yy+', length=len, /fold_case)
        while pos ne -1 do begin
            strput, fmt, strmid(string(yr,format='(I04)'),4-len,len), pos
            pos = stregex(fmt,'yy+', length=len, /fold_case)
        endwhile
        ; date.
        pos = stregex(fmt,'dd', /fold_case)
        while pos ne -1 do begin
            strput, fmt, string(dy,format='(I02)'), pos
            pos = stregex(fmt,'dd', /fold_case)
        endwhile
        ; day of year.
        pos = stregex(fmt,'doy', /fold_case)
        while pos ne -1 do begin
            strput, fmt, string(stodoy(yr,mo,dy),format='(I03)'), pos
            pos = stregex(fmt,'doy', /fold_case)
        endwhile

        ; hour.
        pos = stregex(fmt,'hh', /fold_case)
        while pos ne -1 do begin
            strput, fmt, string(hr,format='(I02)'), pos
            pos = stregex(fmt,'hh', /fold_case)
        endwhile
        ; minute.
        pos = strpos(fmt,'mm')
        while pos ne -1 do begin
            strput, fmt, string(mi,format='(I02)'), pos
            pos = strpos(fmt,'mm')
        endwhile
        ; sec.
        pos = strpos(fmt,'ss')
        while pos ne -1 do begin
            strput, fmt, string(sc,format='(I02)'), pos
            pos = strpos(fmt,'ss')
        endwhile
        pos = stregex(fmt, '\.f+',length=len, /fold_case)
        while pos ne -1 do begin
            tmp = 1d-3*(ms+1d-3*(mc+1d-3*pc))
            tmp2 = '(F'+strtrim(string(len+2),2)+'.'+strtrim(string(len),2)+')'
            strput, fmt, strmid(string(tmp,format=tmp2),1,len), pos
            pos = stregex(fmt, '\.f+',length=len, /fold_case)
        endwhile
        ; milli sec.
        pos = strpos(fmt,'lll')
        while pos ne -1 do begin
            strput, fmt, string(ms,format='(I03)'), pos
            pos = strpos(fmt,'lll')
        endwhile
        ; micro sec.
        pos = strpos(fmt,'ccc')
        while pos ne -1 do begin
            strput, fmt, string(mc,format='(I03)'), pos
            pos = strpos(fmt,'ccc')
        endwhile
        ; pico sec.
        pos = strpos(fmt,'ppp')
        while pos ne -1 do begin
            strput, fmt, string(pc,format='(I03)'), pos
            pos = strpos(fmt,'ppp')
        endwhile
        
        ; finally.
        etstrs[i] = fmt
    endfor
    return, etstrs
end

print, 'Test for Epoch ...'
cdf_epoch, et, 2014, 05, 19, 21, 06, 07, 123, /compute_epoch
print, sfmepoch(et)
print, sfmepoch(et,'yyYY-MM-DD/doy/hh:mm:ss.lll.ccc.ppp')
print, sfmepoch(et,'yyYY-MM-DD/doy/hh:mm:ss.fffffffffffff')
print, sfmepoch(et,'yyYY-MM-DD/doy/hh:mm:ss.ff')
print, time_string(sfmepoch(et,'unix'))

print, 'Test for Epoch16 ...'
cdf_epoch16, et, 2014, 05, 19, 21, 06, 07, 123, 456, 789, /compute_epoch
print, sfmepoch(et,/epoch16)
print, sfmepoch(et,/epoch16, 'yyYY-MM-DD/doy/hh:mm:ss.lll.ccc.ppp')
print, sfmepoch(et,/epoch16, 'yyYY-MM-DD/doy/hh:mm:ss.fffffffffffff')
print, sfmepoch(et,/epoch16, 'yyYY-MM-DD/doy/hh:mm:ss.ff')
print, time_string(sfmepoch(et,'unix'))

print, 'Test for arrays ...'
ets = [et,et]
print, time_string(sfmepoch(ets,'unix'),tformat='yy-MM-DD/hh:mm:ss.fffffffff')
end