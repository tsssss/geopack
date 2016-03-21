;+
; Type: function.
; Purpose: Breaks down time range into given time difference.
; Parameters:
; 	ets, in, dblarr[2], required. Time range in epoch.
; 	dt, in, double, optional. Time difference. Default in days.
; Keywords: none.
; Return: dblarr[2]/dblarr[2,n]. Each two pairs are [start, end] epochs.
; Notes: none.
; Dependence: none.
; History:
; 	2013-04-04, Sheng Tian, create.
;-
function sbreaktr, ets, dt
	compile_opt idl2

    if n_elements(dt) eq 0 then dt = 86400000d     ; 1 day.
;    dt0 = 1          ; to make time range do not equal to upper limit.
	et0 = ets[0]-(ets[0] mod dt)
	et1 = ets[1]-(ets[1] mod dt)
    netr = floor((et1-et0)/dt)
    if (ets[1] mod dt) eq 0 then netr-= 1
	if netr lt 0 then message, 'start time > end time ...'

	; only 1 time range.
	if netr eq 0 then return, [ets[0],ets[1]]
	; more 1 time range.
	tr0s = et0+findgen(netr+1)*dt
	tr1s = tr0s+dt
	tr0s[0] = ets[0]
	tr1s[netr] = ets[1]

	return, transpose([[tr0s],[tr1s]])
end
etr = stoepoch(['2008-05-01/13:33:26','2008-05-01/15:00'])
etr = sbreaktr(etr,3600000d)
print, sfmepoch(etr[*])
end
