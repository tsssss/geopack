;+
; Type: procedure.
; 
; Purpose: Print IDL !path to console or file.
; 
; Parameters:
;   path, in, string, optional. Path to be printed. Default is IDL !path.
; 
; Keywords:
;   filename, in, string, optional. File to which path are preinted. Omit it 
;       will print paths to console.
; 
; Notes: none.
; 
; Dependence: none.
; 
; History:
;   2012-05-17, Sheng Tian, create.
;-

pro slprpath, path, filename = file

    console = -1
    if n_elements(path) eq 0 then path = !path
    if n_elements(file) ne 0 then openw, lun, file, /get_lun else lun = console

    ; separator for !path, OS dependent.
    sep = path_sep(/search_path)
    
    paths = strsplit(path, sep, /extract)
    npath = n_elements(paths)
    printf, lun, 'number of paths: ', npath
    
    for ii = 0, npath-1 do printf, lun, paths[ii]
    
    if lun ne console then free_lun, lun
    
end
