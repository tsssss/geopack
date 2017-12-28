;+
; Type: function:
; Purpose: Make geometic progression from x0 to x1, x1 / x0 = dx^n.
;   * mode = 'x0', a = x0, b = dx, c = n
;   * mode = 'x1', a = x1, b = dx, c = n.
;   * mode = 'dx', a = x0, b = x1, c = dx.
;   * mode = 'n',  a = x0, b = x1, c = n.
; Parameters:
;   a, in, double, req. Lower limit of the array.
;   b, in, double, req. Upper limit of the array.
;   c, in, long/double, req. Number of records ('n' mode) or step ('dx' mode).
;   mode, in, string, opt. Use mode to specify various usages. 'dx' by default.
; Keywords: none.
; Return: double. Generated array.
; Notes: none.
; Dependence: none.
; History:
;   2012-09-19, Sheng Tian, create.
;   2012-10-29, Sheng Tian, add mode keyword.
;-

function smkgmtrc, a, b, c, mode
    if n_elements(mode) eq 0 then mode = 'dx'
    case mode of
        'x0' : begin
            arr = double(a)
            for i = 0, c-2 do arr = [arr, arr[i]*b]
        end
        'x1' : begin
            arr = double(a) & q = 1D/b
            for i = 0, c-2 do arr = [arr[0]*q, arr]
        end
        'dx' : begin
            ns = alog(double(b)/a)/alog(double(c)) & arr = a
            for i = 0, ns-1 do arr = [arr, arr[i]*c]
        end
        'n'  : begin
            q = (double(b)/a)^(1D/(c-1)) & arr = a
            for i = 0, c-2 do arr = [arr, arr[i]*q]
        end
    endcase
    return, arr
end
