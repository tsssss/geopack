;+
; Type:
; 	function.
; 
; Name:
; 	sang.
; 
; Purpose:
; 	Return angle between two vectors.
; 
; Parameters:
; 	vec1, in, type = dblarr[n,m], required.
;		Vector 1, m-dimensional, has n records.
; 
; 	vec2, in, type = dblarr[n,m], required.
;		Vector 2, m-dimensional, has n records.
;
; Keywords:
; 	degree = degree, in, type = boolean, optional.
;		Set it to return angle in degree.
;		Default is in radian.
; 
; Return:
; 	return, out, type = dblarr[n].
;		Angle between vec1 and vec2.
; 
; Example:
; 	angle = sang(vec1, vec2, /degree).
; 
; Notes:
; 	* if n = 1 or 0, then return is in double.
; 
; Dependence:
; 	none.
; 
; Author:
; 	Sheng Tian.
; 
; History:
; 	2013-04-18, Sheng Tian, create.
;-

function sang, vec1, vec2, degree = deg

	compile_opt idl2

	sz = size(vec1)
	if sz[0] gt 2 then message, 'vector dim > 2 ...'

	if sz[0] eq 1 then $
		ang = acos(total(vec1*vec2)/sqrt(total(vec1^2)*total(vec2^2))) $
    else $
	    ang = acos(total(vec1*vec2,2)/sqrt(total(vec1^2,2)*total(vec2^2,2)))

	if keyword_set(deg) then return, (1d/!dtor)*ang else return, ang

end
