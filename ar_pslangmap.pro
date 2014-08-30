;Determine the angle of the PSL at each point
;Use 3 or 5 point smoothing??

;All psl pixels should be = 1

function ar_pslangmap, pslmask, params=inparams

if data_type(inparams) eq 8 then params=inparams $
	else params=ar_loadparam() ;get the default SMART parameter list

radpslang=params.radpslang

imgsz=size(pslmask,/dim)

shiftval=imgsz/2

wpsl=where(pslmask ne 0)

kernelimg=dist(imgsz[0],imgsz[1]) ;,shiftval[0],shiftval[1])
wkgt=where(kernelimg gt radpslang)
wklt=where(kernelimg le radpslang)
kernelimg[wkgt]=0
kernelimg[wklt]=1


for i=0,wpsl-1 do begin
		
	thisx=wpsl mod imgsz[0]
	thisy=wpsl/imgsz[0]
	
	thiskern=shift(kernelimg,thisx,thisy)
	
	
	;shift kernel to the location of the pixel
	
	thismask=pslmask*thiskern
	
	wthispslpts=where(thismask eq 1)
	
	
	;might need to find a more smoothed out PSL detection
	
	
	
	;perhaps use core psl??
	
	
	
	
	;fit a line through the psl pixels that it finds
	
	thisang=(atan(linfit(wthispslpts mod imgsz[0],wthispslpts/imgsz[0])))[0]/!dtor
	
	
	
	
	
	
stop	
	
	
	
	
	
	
	
	
endfor







































return, pslangmap

end