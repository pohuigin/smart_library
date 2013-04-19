;SMART_GROW 
;Returns a dilated mask. If a NL mask is provided and GAUSSIAN is set, 
;then the result will be a Shrijver R-mask.
;
;Provide RADIUS or FWHM in pixels. FWHM is actually half width at half max!!!
;If GAUS is set, then the convolution stucture will fall off as a gaussian.
;If RADIUS is set then the FWHM of the gaussian will be half that.
;If FWHM is set, then RADIUS will be twice that.

function ar_grow, arr, radius=radius, gaus=gaus, fwhm=fwhm, _extra=_extra, kernal=outkernal

arr0=arr

if n_elements(radius) lt 1 then radius0=5. else radius0=radius

if n_elements(fwhm) gt 0 then begin
	radius0=fwhm;*2.
endif else fwhm=radius0;/2.
gsig=fwhm/(SQRT(2.*ALOG(2.)))

if keyword_set(gaus) then imgsz=[2,4.*radius0,4.*radius0] else imgsz=[2,3.*radius0,3.*radius0]

;make sure the kernal has an odd number of elements
if imgsz[1] mod 2 eq 0 then imgsz[1]=imgsz[1]+1
if imgsz[2] mod 2 eq 0 then imgsz[2]=imgsz[2]+1

struc=fltarr(imgsz[1],imgsz[2])

;Generate coordinate maps.
xcoord=rot(rebin(transpose(findgen(imgsz[1])),imgsz[1],imgsz[2]),90)
ycoord=rot(xcoord,-90)
rcoord=sqrt((xcoord-imgsz[1]/2.)^2.+(ycoord-imgsz[2]/2.)^2)

;The kernal apparently needs to be shifted left and down by 1 pixel
rcoord=shift(rcoord,[-1,-1])

struc[where(rcoord le radius0)]=1.

outkernal=struc

;stop

if keyword_set(gaus) then begin
	;gparams[0] = maximum value (factor) of Gaussian,
	;gparams[1] = mean value (center) of Gaussian,
	;gparams[2] = standard deviation (sigma) of Gaussian.
	
	gparams=[1.,0.,gsig]
	
;	dummy=gaussian([1,2,3], gparams)
	gstruc=gaussian(rcoord, gparams)
	gstruc=reform(gstruc,imgsz[1],imgsz[2])
	
	;;Normalize GSTRUC so that the "central" or max value of GROWNARR is 1.
	;gstruc=gstruc/total(radius0*2.)

	;Normalize GSTRUC so that the volume is 1
	gstruc=gstruc/total(gstruc)
	outkernal=gstruc

	grownarr=CONVOL( arr0, gstruc, _extra=_extra); < 1.
	
	;;grownarr=dilate(arr0,gstruc)

endif else grownarr=dilate(arr0,struc)

return, grownarr

end