;Routine to process a LOS magnetogram for detecting ARs and determining magnetic properties.
;Runs:
;	Remove non-finites
;	Zero off-limb pixels 
;	Does Cosine correction
;	Do median filter
; 	Rotate solar north = up
function ar_processmag, inmap, limbmask=limbmask, $
	nocos=nocos, nofilter=nofilter, nofinite=nofinite, noofflimb=noofflimb, norotate=norotate

map=inmap
dat=map.data

;stop

;Clean NaNs
if not keyword_set(nofinite) then $
	if (where(finite(dat) ne 1))[0] ne -1 then dat[where(finite(dat) ne 1)]=0.

;Get the cosine map and off-limb pixel map using WCS
cosmap=ar_cosmap(map, rrdeg=rrdeg, offlimb=offlimb,/edge)
limbmask=offlimb

;zero off-limb pixels
;zero from 80 degrees to LOS
if not keyword_set(noofflimb) then begin
	wofflimb=where(rrdeg/2./!pi*360. ge 80.)
	if wofflimb[0] ne -1 then dat[wofflimb]=0.
	dat=dat*offlimb
endif

;Median filter noisy values
;Do a 3x3 median filter
if not keyword_set(nofilter) then begin
	dat=filter_image(dat,/MEDIAN)
endif

;Do magnetic field cosine correction
;Limit correction to having 1 pixel at edge of the Sun
if not keyword_set(nocos) then begin
	coscorlim=ar_coscorlim(map)
	cosmap=cosmap < coscorlim
	dat=dat*cosmap
endif

map.data=dat

;Rotate solar north = up
if not keyword_set(norotate) then $
	map=rot_map(map,(-map.roll_angle))

return, map

end
