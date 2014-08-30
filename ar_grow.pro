;SMART_GROW 
;Returns a dilated mask. If a NL mask is provided and GAUSSIAN is set, 
;then the result will be a Shrijver R-mask.
;
;Provide RADIUS or FWHM in pixels. FWHM is actually half width at half max!!!
;If GAUS is set, then the convolution stucture will fall off as a gaussian.
;If RADIUS is set then the FWHM of the gaussian will be half that.
;If FWHM is set, then RADIUS will be twice that.
;
;Notes:
;	1. For nice circular kernel binary kernel, width of kernel mask will be 2*radius+1, with a 1px boundary of 0 around the outside.
;	2. Setting radius to 1 will result in a 3x3 structuring element for binary kernels, with a total array size of 5x5

function ar_grow, arr, radius=radius, gaus=gaus, fwhm=fwhm, _extra=_extra, $
	kernal=outkernal, inkernal=inkernal

arr0=arr

if n_elements(radius) lt 1 then radius0=5. else radius0=radius

if n_elements(fwhm) gt 0 then begin
	radius0=fwhm;*2.
endif else fwhm=radius0;/2.
gsig=fwhm/(SQRT(2.*ALOG(2.)))

if radius0 eq 1 then imgsz=[2,6,6] else $
	if keyword_set(gaus) then imgsz=[2,4.*radius0,4.*radius0] else imgsz=[2,4.*radius0,4.*radius0] ;want it to be even. ; imgsz=[2,3.*radius0,3.*radius0]

;make sure the kernal has an odd number of elements
;if imgsz[1] mod 2 eq 0 then imgsz[1]=imgsz[1]+1
;if imgsz[2] mod 2 eq 0 then imgsz[2]=imgsz[2]+1

struc=fltarr(imgsz[1],imgsz[2])

;Generate coordinate maps.
;xcoord=rot(rebin(transpose(findgen(imgsz[1])),imgsz[1],imgsz[2]),90)
;ycoord=rot(xcoord,-90)
;rcoord=sqrt((xcoord-imgsz[1]/2.)^2.+(ycoord-imgsz[2]/2.)^2)

xyrcoord,imgsz,xcoord,ycoord,rcoord
struc[where(rcoord le radius0)]=1.

;Crop to the edges of kernel with 1px boundary
wxbound=minmax(where(total(struc,2)))
wybound=minmax(where(total(struc,1)))
struc=struc[wxbound[0]-1:wxbound[1]+1,wybound[0]-1:wybound[1]+1]

;The kernal apparently needs to be shifted left and down by 1 pixel
;rcoord=shift(rcoord,[-1,-1])



;Check for input kernal
if n_elements(inkernal) gt 0 then begin
	struc=inkernal
endif

wnogood=where(finite(struc) ne 1)
if wnogood[0] ne -1 then struc[wnogood]=0

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

;Check for input kernal
	if n_elements(inkernal) gt 0 then begin
		gstruc=inkernal
	endif

	wnogood=where(finite(gstruc) ne 1)
	if wnogood[0] ne -1 then gstruc[wnogood]=0

	outkernal=gstruc
	
	if min(size(arr0,/dim)) gt min(size(gstruc,/dim)) then grownarr=CONVOL( arr0, gstruc, _extra=_extra) $
		else begin
			print,'% AR_GROW: KERNEL is too big, compared to IMAGE!'
			return, arr
		endelse


	
	
	;;grownarr=dilate(arr0,gstruc)

endif else begin
	
	if min(size(arr0,/dim)) gt min(size(struc,/dim)) then grownarr=dilate(arr0,struc) $
		else begin
			print,'% AR_GROW: KERNEL is too big, compared to IMAGE!'
			return, arr
		endelse

endelse

return, grownarr

end