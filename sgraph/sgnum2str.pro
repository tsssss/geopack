;+
; Type: function.
; Purpose: Convert number to string, no white space on both sides.
; Parameters:
; 	num, in, int/long/float/double, required. The given number.
; Keywords:
; Return: string. String of the number.
; Notes: Do not deal with complex(6) and dcomplex(9), maybe later.
;   Any settings, includes 3 flavors: natural, maximum, and required. For 
;   example, we can set how many decimals for 1.234. The natural form is 1.234,
;   which has 3 decimals. If you are not happy with it, you can say, maximum
;   decimal to be 2, then you get 1.23 for 1.234, 1.24 for 1.235 (round), 1.2
;   for 1.2 (it's natural form doesn't vialent the maximum decimal of 3).
;   Lastly, you can force certain number of decimals, say 2. Then you get
;   1.20 for 1.2, 1.23 for 1.234, 1.24 for 1.235.
; Dependence: none.
; History:
; 	2015-10-15, Sheng Tian, create based on snum2str.
;-

; nsigdigits  ; # of significant digit.
; shortest
; fixwidth

function sgnum2str, num0, ndec = ndec0, mdec = mdec0, $
    nsgn = nsgn0, msgn = msgn0, sci = sci, norm = norm
    
	compile_opt idl2
	
	if n_params() lt 1 then return, ''
	
	; num0 can be pointer.
	num = num0
	if size(num,/type) eq 10 then repeat num = *num until size(num,/type) ne 10
		
	tcode = size(num,/type)           ; type code.
	case tcode of
	   0: return, ''   ; undefined.
	   ; 1:byte, 2:int, 3:long, 4:float, 5:double
	   6: message, 'do not deal with complex ...'
	   7: return, ''   ; string.
	   8: message, 'do not deal with structure ...'
	   9: message, 'do not deal with complex ...'
	   10: ; impossible b/c we've treated pointer.
	   11: message, 'do not deal with object ...'
	   ; 12: uint, 13:ulong, 14:lon64, 15:ulon64.
	   else: ; do nothing.
	endcase

	isint = tcode ne 4 and tcode ne 5    ; not float or double
	iszero = num eq 0                   ; 0 needs special care.
	numsgn = (num ge 0)*2-1             ; sign, -1 for negativ, 1 for positive.
	num = double(num)                   ; convert to double.
	numabs = abs(num)                   ; absolute value.
    numexp = floor(alog10(numabs))      ; exponent.
    explen = floor(alog10(abs(numexp)))+1
    if numexp lt 1 then explen+= 1      ; length of exponent, includes sign.
    numdgt = num/10d^numexp             ; 1.345.

    if iszero then begin
        numexp = 0
        numdgt = 0d
    endif

    ; for integer, we can do
    ; 1, print as it is.
    ; 2, add more significant digits after decimal point.

    ; for scientific notation, we can do
    ; 1, "sthx10!UEXP!N", "x10!UEXP!N" part is 3+0.5*n, n is length of the exp.
    ; 2, "stheEXP", "eEXP" part is 1+n, n is length of the exp.
    
	; for float, we can do
	; 1, round to integer.
	; 2, scientific notation showing n significant digits.
	
    ; default settings.
	ndec = 1        ; 1 digit after decimal point.
	if isint then ndec = 0 ; integers do not have decimals.
	maxndec = 3     ; maximum of 3 digits after decimal point.
	maxlen = 5          ; maximum length of the final string, do not include sign.
    if num gt 1 then begin
        maxlen = 4+explen                       ; 3.4e9.
        if explen gt 1 then maxlen = 2+explen   ; 5e23.
    endif else begin
        maxlen = 4+explen                       ; 1.2e-5.
        if explen gt 2 then maxlen = 4+explen   ; 9.1e-31.
    endelse
	
	; what you get from IDL.
    str0 = strtrim(string(num0),2)
    ; remove trailing 0 in decimal part, and unnecessary decimal point.
    pos = strpos(str0, '.')
    if pos ne -1 then begin
        while strmid(str0,strlen(str0)-1) eq '0' do $
            str0 = strmid(str0,0,strlen(str0)-1)
        if strmid(str0,strlen(str0)-1) eq '.' then $
            str0 = strmid(str0,0,strlen(str0)-1)
    endif
    if n_elements(ndec0) eq 0 and n_elements(nsgn0) eq 0 $
    and n_elements(mdec0) eq 0 and n_elements(msgn0) eq 0 then return, str0
    str0nsgn = strlen(strjoin(strsplit(str0,'-.',/extract),''))
    
    ; when ndec is set.
    if n_elements(ndec0) ne 0 then ndec = ndec0
    pos = strpos(str0, '.')
    str0declen = (pos eq -1)? 0: strlen(str0)-pos-1
    if ndec ge str0declen then begin    ; add 0s to # of decimal digits.
        str1 = (pos eq -1)? str0+'.': str0  ; add decimal point if needed.
        if ndec eq 0 then str1 = str0   ; 0 decimal case.
        for i = 0, ndec-str0declen-1 do str1 = str1+'0'
    endif else begin                        ; round the number.
        numdec = num mod 1
        base = 10d^ndec
        num1 = num-numdec+round(numdec*base)/base
        str1 = strtrim(string(num1),2)
        pos = strpos(str1, '.')
        if pos ne -1 then begin
            if ndec eq 0 then str1 = strmid(str1,0,pos) $   ; no decimal point.
            else str1 = strmid(str1,0,pos+ndec+1)
        endif
    endelse
    if n_elements(ndec0) ne 0 then return, str1

    ; when nsgn is set.
    nsgn = 2
    if n_elements(nsgn0) ne 0 then nsgn = nsgn0
    if n_elements(msgn0) ne 0 then nsgn = str0nsgn<msgn0
    base1 = 10d^numexp          ; convert number into (-2,2).
    base2 = 10d^(numexp-nsgn+1) ; convert number into xxx, # of digit is nsgn.
    num2 = round(num/base2)*base2   ; around to meet # of significant digits.
    ; scientific notation.
    str21 = strtrim(string(num2/base1),2)
    pos = strpos(str21, '.')
    str2declen = strlen(str21)-pos-1
    if nsgn-1 ge str2declen then begin      ; add 0s to # of significant digits.
        if nsgn ne 1 then begin
            for i = 0, nsgn-str2declen-2 do str21 = str21+'0'
        endif
    endif else begin
        if nsgn eq 1 then str21 = strmid(str21,0,pos) $
        else str21 = strmid(str21,0,pos+nsgn)
    endelse
    strexp = (numexp ne 0)? 'e'+strtrim(string(numexp),2): ''
    str21 = str21+strexp
    
    ; normal notation.
    str22 = strtrim(string(num2),2)
    pos = strpos(str22, '.')
    if pos ne -1 then begin     ; remove trailing 0s and decimal point.
        while strmid(str22,strlen(str22)-1) eq '0' do $
            str22 = strmid(str22,0,strlen(str22)-1)
        if strmid(str22,strlen(str22)-1) eq '.' then $
            str22 = strmid(str22,0,strlen(str22)-1)
    endif
    str22ndigit = strlen(strjoin(strsplit(str22,'-.',/extract),''))
    if nsgn gt str22ndigit then begin     ; add 0s to meet # of sgn digits.
        if strpos(str22, '.') eq -1 then str22 = str22+'.'
        for i = 0, nsgn-str22ndigit-1 do str22 = str22+'0'
    endif
    
    str2 = str22
    if strlen(str2)-strlen(str21) gt 1 then str2 = str21  ; prefer normal notation.
    
    if keyword_set(sci) then begin
        str2 = str21
        if strpos(str2, 'e') eq -1 then str2 = str2+'e1'
    endif
    
    if keyword_set(norm) then str2 = str22 
    
    if n_elements(nsgn0)+n_elements(msgn0) ne 0 then return, str2
        
end

print, sgnum2str(0.001)
print, sgnum2str(120000, ndec = 1)
print, sgnum2str(120000, nsgn = 2)
print, sgnum2str(-120000)
print, sgnum2str(-12000)
print, sgnum2str(-1.23)
print, sgnum2str(1.73)
print, sgnum2str(-120000)
print, sgnum2str(120000)
end
