;+
; Type: function.
; Purpose: Returns flag on if given tag exists.
; Parameters:
;   tag, in, string, req. The tag to be checked.
;   struct, in, struct, req. The structure to be checked.
; Keywords: none.
; Return: boolean. 1: tag exists, 0:tag does not exist.
; Notes: Tag is case-insensitive.
; Dependence: none.
; History:
;   2012-06-29, Sheng Tian, create.
;-

function stagexist, tag, struct
  compile_opt idl2 & on_error, 2
  
  if size(struct, /type) ne 8 then return, 0
  
  thetag = strupcase(tag)
  tags = tag_names(struct)
  
  if where(tags eq thetag) ne -1 then return, 1 $
  else return, 0
   
end
