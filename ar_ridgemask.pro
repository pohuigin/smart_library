;Pull out a ridge skeleton and insert into a mask to create a 1 pixel-wide trace for 
;determining the PSL mask, using a watershed transform
;
;Example: ridgemask=ar_ridge(abs(datagrad))
;
;STATUS=	0 - initialised value
;			1 - no PSL/watershed pixels found
;			2 - PSL pixels found

function ar_ridgemask,data,status=status, thresh=inthresh
status=0

if n_elements(inthresh) eq 1 then thresh=inthresh else thresh=20

imgsz=size(data,/dim)

datwshed = WATERSHED(-(data>thresh), connect=8) ; [, CONNECTIVITY={4 | 8} ] [, /LONG] [, NREGIONS=variable] ) 

blank=fltarr(imgsz[0],imgsz[1])

msk=blank

wwshed=where(datwshed eq 0)

if wwshed[0] eq -1 then begin
	status=1
	return,blank
endif

msk[wwshed]=1

status=2

outmsk=msk

return,outmsk

end