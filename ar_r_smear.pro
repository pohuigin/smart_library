;+
;  NAME:
;    ar_r_smear
;  PURPOSE:
;    convolve an image with [default] a gaussian profile of FWHM width
;    n and boxwidth 4n, or alternatively with a specified kernel
;
;  CALLING SEQUENCE:
;    image=ar_r_smear(image,n,kernel=kernel)
;
;  INPUTS:
;    image  image to be processed
;    n      fwhm value of the gaussian smearing that is applied
;
;  OUTPUTS:
;           returns convolved image
;
;  OPTIONAL INPUT KEYWORDS:
;    kernal        convolution kernel, 2d real array
;   [edge_truncate as in IDL routine convol - obsolete 2006/04/17, see below]
;
;  OPTIONAL OUTPUT KEYWORDS:
;
;  METHOD:
;    uses standard IDL routine convol
;
;  MODIFICATION HISTORY:
;	C.J. Shcrijver - 11-Feb-2014 - written
;	P.A. Higgins - 12-Feb-2014 - modified using M. Bobra's suggestion (changed 
;		kernal width to 4*n+1 rather than 4n) and standardised code to fit
;		within the SMART_LIBRARY repository: 
;		http://github.com/pohuigin/smart_library/
;
;-

function ar_r_smear,image, kernel=kernel, szkernel=inszkernel,edge_truncate=edge_truncate,nan=nan, param=param, fparam=fparam

;kernel=inkernel

if n_elements(param) eq 0 then param=ar_loadparam(fparam=fparam)

if n_elements(inszkernel) eq 0 then szkernel=param.r_kernsz else szkernel=inszkernel
n=szkernel

if n eq 0 then return,image

;if n mod 2 eq 0 then n=n+1

sigma=n/(2.*sqrt(2.*alog(2.)))

if n_elements(kernel) eq 0 then begin

	kernel=fltarr(4*fix(n)+1.)

	for i=0,4*fix(n) do kernel[i]=exp(-(i-float(2*fix(n)-.5))^2/(2*sigma^2))

	kernel=kernel#kernel

endif else n=(size(kernel))[1]/2

kernel=kernel/total(kernel)

; there is a non-repeating bug in convol on IDL5.6 on solserv, which
; creates an image artefact at the top of the fov for smoothing of
; synoptic maps. To avoid that, an explicit buffer zone is formed and
; later removed

sx=(size(image))[1]
sy=(size(image))[2]
padim=fltarr(sx+2*6*fix(n),sy+2*6*fix(n))
padim[6*fix(n),6*fix(n)]=image

;smeared=convol(float(padim),kernel,/edge_truncate,/nan)
smeared=convolve(float(padim),kernel)

out=smeared[6*fix(n):6*fix(n)+sx-1,6*fix(n):6*fix(n)+sy-1]

; restore the original edge of the image - now obsolete!
;if not(keyword_set(edge_truncate)) then begin
;  out=float(image)
;  out(2*n,2*n)=smeared(2*n:(size(image))(1)-2*n,2*n:(size(image))(2)-2*n)
;endif else out=smeared


return,out
end
