;+
; Type: function.
; Purpose: Calculate Poynting flux from E and B field.
; Parameters:
;   de, in, dblarr[n,d]/dblarr[n,m,d], req. dE field in mV/m.
;       It is [n] when it is 1-d dE, in [n,2] for 2-d dE,
;       in [n,3] if dE is 3-d field.
;       If de is spectrogram, then it is in [n,m] or [n,m,d],
;       where m is # of scales.
;   db, in, dblarr[n,d]/dblarr[n,m,d], req. dB field in nT. See above.
;   coeff, in, double, opt. Set coeff to fix units.
;       Default is 1/(400*pi), gives Poynting flux in ergs/cm^2/s.
; Keywords:
;   spec, in, boolean, opt. Set to indicate de0 and db0 are spectrograms.
; Return: dblarr[n]/dblarr[n,3]. Poynting flux in ergs/cm^2/s.
;       It's in [n] if db and de are in [n], or in [n,2];
;       It's in [n,3] if db and de are in [n,3].
; Notes:
;   * For 2-d electric field, check the sign!!!
;   * For spectrogram, only 1-d and 2-d are supported.
; Dependence: none.
; History:
;   2010-07-27, Sheng Tian, create.
;   2012-10-01, Sheng Tian, rename.
;   2013-04-17, Sheng Tian, revise.
;-

function spoynt, de0, db0, c, spec = spec
    compile_opt idl2

    de = de0
    db = db0

    ; default constant.
    if n_elements(c) eq 0 then c = 1D/(400D*!dpi)

    sz = size(de)
    if sz[0] gt 3 then message, 'field dim > 2 ...'

    if keyword_set(spec) then begin
        ; 1-d field spectrogram in [n,m].
        if sz[0] eq 2 then return, (de*db)*c
        if sz[0] ne 3 then message, 'spec wrong dim ...'
        if sz[3] eq 2 then return, $
            (de[*,*,0]*db[*,*,1]-de[*,*,1]*db[*,*,0])*c
    endif
        
    ; 1-d field.
    if sz[0] eq 1 then return, (de*db)*c

    ; 2-d field.
    if sz[2] eq 2 then return, (de[*,0]*db[*,1]-de[*,1]*db[*,0])*c
    
    ; 3-d field.
    if sz[2] eq 3 then return, scross(de,db)*c
end
