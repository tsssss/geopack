;+
; Type: function.
; Purpose: Convert CDF type to IDL type.
; Parameters:
;   cdftype, in, string/strarr[n], req. CDF type(s).
; Keywords: none.
; Return: int/intarr[n]. IDL data type code.
; Notes: none.
; Dependence: none.
; History:
;   2014-02-12, Sheng Tian, create.
;-
function scdf2idltype, types
    on_error, 2

    cdftypes = 'CDF_'+['XXX','BYTE','UINT1','INT1','CHAR','UCHAR',$
        'INT2','UINT2','INT4','UINT4','REAL4','FLOAT','DOUBLE','REAL8',$
        'EPOCH','EPOCH16','LONG_EPOCH']
    idltypes = [0,1,1,1,1,1,2,12,3,13,4,4,5,5,5,9,9]

    types = strupcase(types)
    ntype = n_elements(types)
    codes = intarr(ntype)
    for i = 0, ntype-1 do begin
        idx = where(cdftypes eq types[i],cnt)
        codes[i] = cnt? idltypes[idx[0]]: 0
    endfor
    if ntype eq 1 then codes = codes[0]
    return, codes
end
