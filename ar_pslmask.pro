;Create a PSL mask, given some LOS B input data
;
;Example: pslmask=ar_psl(datasm,radius=smoothhwhm,thresh=thresh)
;
;RADIUS = dilation radius in pixels
;THRESH = magnetic threshold in G
;STATUS = 	0 - initialised value
;			1 - a PSL detection was found
;			2 - no positive or negative pixels found above threshold
;			3 - no overlap between dilated positive and negative masks

function ar_pslmask, data, radius=inrad, thresh=inthresh, status=status
status=0

if n_elements(inrad) eq 1 then rad=inrad else rad=5.

if n_elements(inthresh) eq 1 then thresh=inthresh else thresh=100.

imgsz=size(data,/dim)

blank=fltarr(imgsz[0],imgsz[1])

;make a negative and positive mask

wpos=where(data gt thresh)

wneg=where(data lt -thresh)

if wpos[0] eq -1 or wneg[0] eq -1 then begin
	status=2
	return,blank
endif

posmsk=blank
posmsk[wpos]=1.

negmsk=blank
negmsk[wneg]=1.

;dilate the polarity masks

posmskgrw=ar_grow(posmsk, radius=rad)

negmskgrw=ar_grow(negmsk, radius=rad)

;Determine the overlap of the two masks

ovrmsk=posmskgrw+negmskgrw

wover=where(ovrmsk eq 2)

if wover[0] eq -1 then begin
	status=3
	return,blank
endif

outmsk=blank
outmsk[wover]=1.

status=1

return, outmsk

end