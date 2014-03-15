;Take the gradient in the horizontal plane of the LOS magnetic field.
;
;Example: datagrad=ar_losgrad(datasm)
;

function ar_losgrad, data

;Buffer the image to reduce edge effects
imgsz=size(data,/dim)

;xinterp=[findgen(imgsz[0]+1),fltarr(imgsz[1]+1)+imgsz[0],findgen(imgsz[0]+1)-1,fltarr(imgsz[1]+1)-1]
;yinterp=[fltarr(imgsz[0]+1)+imgsz[1],findgen(imgsz[1]+1)-1,fltarr(imgsz[0]+1)-1,findgen(imgsz[1]+1)]

dataint=fltarr(imgsz[0]+10,imgsz[1]+10)-1d6
dataint[5:imgsz[0]+4,5:imgsz[1]+4]=data

fill_missing,dataint,-1d6,1,/extrap

fill_missing,dataint,-1d6,2,/extrap

;winterp=where(dataint eq -1d6)
;xinterp=winterp mod imgsz[0]+10
;yinterp=winterp/(imgsz[0]+10)

;interps=interpolate(data,xinterp,yinterp)

;dataint[xinterp+5,yinterp+5]=interps

xgrad=deriv(dataint)

ygrad=rot(deriv(rot(dataint,-90)),90)

gradmag=sqrt(xgrad^2.+ygrad^2.)

return, gradmag[5:imgsz[0]+4,5:imgsz[1]+4]

end