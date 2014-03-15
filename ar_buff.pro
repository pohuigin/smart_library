;buffer an image or map by 100 px for doing 
;widthbuff=width of the buffer going around the image

function ar_buff, indat, widthbuff=inbuff, valbuff=inval, inverse=inverse

if data_type(indat) eq 8 then begin
	image=indat
	domap=1
endif else begin
	image=indat
	domap=0
endelse

if n_elements(inbuff) eq 1 then buff=inbuff else buff=100

if n_elements(inval) eq 1 then val=inval else val=0.

imgsz=size(image,/dim)

if keyword_set(inverse) then begin
	cutout=image[buff:imgsz[0]-buff-1,buff:imgsz[1]-buff-1]

	if domap then begin
		retdat=indat
		add_prop,retdat,data=cutout,/repl
	endif else begin
		retdat=cutout
	endelse

endif else begin
	blank=fltarr(imgsz[0]+2.*buff,imgsz[1]+2.*buff)+val
	
	blank[buff:buff+imgsz[0]-1,buff:buff+imgsz[1]-1]=image
	image=blank

	if domap then begin
		retdat=indat
		add_prop,retdat,data=image,/repl
	endif else begin
		retdat=image
	endelse

endelse




return,retdat

end