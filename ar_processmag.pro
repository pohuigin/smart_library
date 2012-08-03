;Runs:
;	Remove non-finites
;	Does noise threshold
;	Does Cosine correction
function ar_processmag, map, noisethresh=noisethresh, $
	nocos=nocos, nonoise=nonoise, nofinite=nofinite, $
	mdi=mdi

dat=map.data

;Clean NaNs
if not keyword_set(nofinite) then $
	if (where(finite(dat) ne 1))[0] ne -1 then dat[where(finite(dat) ne 1)]=0.

;Threshold noisy values
if not keyword_set(nonoise) then begin
	if n_elements(noisethresh) lt 1 then begin
		params=ar_loadparam()
		if keyword_set(mdi) then noisethresh=params.mdi_noisethresh $
			else begin & print,'% AR_PROCESSMAG: Must input NOISETHRESH.' & return,'' & endelse
	endif

	dat[where(abs(dat) lt noisethresh)]=0.
endif

;Do magnetic field cosine correction
if not keyword_set(nocos) then begin
	coscor=ar_cosmap(map)

	coscorlim=ar_coscorlim(map)
	
	coscor=coscor < coscorlim
	dat=dat*coscor
endif

map.data=dat

outmap=map
return, outmap

end
