;Determine many PSL properties of an AR 
;
;INPUTS:
;	inmap = data map of magnetogram
;	inmask = indexed image mask of detected features
;
;OPT. INPUTS:
;	refimg = an image to be projected in the same manner as the data and mask
;
;KEYWORDS:
;	param = The SMART2 parameter structure 
;	fparam = Filename of desired parameter structure to use (PARAM takes precedence)
;	doproj = set to do a stereographic deprojection when determining PSL props.
;	projscl = DEFUNCT choose the factor to increase the projection image dimensions by, as compared to the original image
;	dobppxscl = set to use the bipole separation length to determine the deprojected pixel scaling (= bp_sep/bp_sep_proj)
;	outproj = ouput the STG projected magnetogram
;	outscl = output the 'true' pixel scaling for each projected pixel using the gradient of Rdeg map
;	outbpscl = output the conversion between great-circle and STG projected bipole separation length to use as conversion between scaling 
;	outmaskproj = STG projected mask
;	outrefproj = STG projected version of input REFIMG
;	projlimbxy = x,y positions of the limb in STG projected space
;	outpslmask = mask of PSL in STG space
;	outgradpsl = gradient image of PSL in STG space
;
;NOTES:
;	1. If /DOPROJ is NOT set, then all outputs labeled 'STG projected' will just be in LOS (HCP) space, corresponding to the input map.

function ar_pslprop, inmap, inmask, refimg, param=param, fparam=fparam, $
	doproj=indoproj, projscl=inprojscl, dobppxscl=dobppxscl, $
	outproj=projmag, outscl=projpxscl, outbpscl=projpxscl_bpsep, outmaskproj=projmask, outrefproj=projref, projlimbxy=projlimbxy, $
	outpslmask=pslmaskt, outgradpsl=gradpsl, projmaxscale=projmaxscale

magmap=inmap
mask=inmask
;arstr=inarstr

;MAKE A PSL PROP STRUCTURE IN THE STRUCTURE PARAM FILE






;Determine the PSL Length, number of strong gradient segments, the PSL curvature


;Determine R


;Determine WLSG (total gradient along PSL)




imgsz=size(magmap.data,/dim)

param=ar_loadparam(fparam=fparam)


if n_elements(indoproj) ne 1 then doproj=param.doproject else doproj=indoproj

if n_elements(inprojscl) ne 1 then projscl=param.projscl else projscl=inprojscl


DATE_OBS= magmap.time

xscale=magmap.dx

;Initialise PSL property structure
arpslstr={arid:0,psllength:0.,pslsglength:0.,pslcurvature:0d,rvalue:0.,wlsg:0.,bipolesep_mm:0.,bipolesep_px:0.,bipolesep_proj:0.}

nar=max(mask)

if nar lt 1 then return,arpslstr

;resize PSL property structure
arpslstr=replicate(arpslstr,nar)

for i=1,nar do begin

thispslstr=arpslstr[i-1]

	thismask=mask
	wni=where(thismask ne i)
	if wni[0] ne -1 then thismask[wni]=0
	wi=where(thismask eq i)
	if wi[0] ne -1 then thismask[wi]=1
	maskmap=magmap
	maskmap.data=thismask
	
	thisdat=magmap
	thisdat.data=magmap.data*thismask

;Check whether there is actually an AR there
if (where(finite(thisdat.data)))[0] eq -1 or (where(thisdat.data ne 0))[0] eq -1 then begin
	print,'% AR_PSLPROP: No AR appears to be present in the data!'
	projmask=-1 & projmag=-1 & projpxscl=-1 & projpxscl_bpsep=-1 & projref=-1 & pslmaskt=-1 & gradpsl=-1
	return, arpslstr
endif

;Take a sub-map around the AR
	sub_map,maskmap,submask,xran=minmax(where(thismask eq 1) mod imgsz[0])+[-1,1],yran=minmax(where(thismask eq 1)/imgsz[0])+[-1,1],/pixel,/noplot
	sub_map,thisdat,submag,ref=submask,/noplot

	map2wcs,submask,wcsmask & add_prop,submask,wcs=wcsmask,/repl
	map2wcs,submag,wcsmag & add_prop,submag,wcs=wcsmag,/repl
	
	if n_elements(refimg) gt 0 then begin
		refmap=maskmap
		refmap.data=refimg
		sub_map,refmap,subref,ref=submask,/noplot
		map2wcs,refmap,wcsref & add_prop,refmap,wcs=wcsref,/repl
	endif
	
	cutoutsz=size(thismask,/dim)
	
;Determine the bipole separation properties	
	bipsepstr=ar_bipolesep(submag)

	if doproj then begin

		projmag=map_hpc2stg(submag, refimg, projlonlatcent=[0.,0.], projscl=projscl, mask=submask.data,projmask=projmask,/doprojscl, projpxscl=projpxscl, outrefproj=projref, projlimbxy=projlimbxy, status=projstatus, projmaxscale=projmaxscale)
		projmask=round(projmask)
		
		if projstatus eq -1 then begin
			doproj=0 

			projpxscl=fltarr(cutoutsz[0],cutoutsz[1])+1.
			projmag=submag.data
			rim=projmag
			projmask=thismask
			
			bipsepstrproj=bipsepstr
			projpxscl_bpsep=1.


		endif else begin

			rim=projmag
			bipsepstrproj=ar_bipolesep(projmag)

;find conversion factor from HC-projected bipole sep. dist. to STG-projected bipole sep dist.
;Use to get a pixel scale conversion factor.
			projpxscl_bpsep=bipsepstr.gcdist_px/bipsepstrproj.pxsep
	
		endelse
	
		
	endif else begin
		projpxscl=fltarr(cutoutsz[0],cutoutsz[1])+1.
		projmag=submag.data
		rim=projmag
		projmask=submask.data
		
		bipsepstrproj=bipsepstr
		projpxscl_bpsep=1.

	endelse


	projsz=size(projmag,/dim)

;Choose whether to use the Rdeg gradient or the bipole separation conversion to determine the projected pixel scaling
	if keyword_set(dobppxscl) then $
		projmmscl=ar_pxscale(submag,/mmppx)*projpxscl_bpsep $ ; $
		else projmmscl=ar_pxscale(submag,/mmppx)*projpxscl ;

	kernpsl=[[0,0,1,1,1,0,0],[0,1,1,1,1,1,0],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1],[1,1,1,1,1,1,1],[0,1,1,1,1,1,0],[0,0,1,1,1,0,0]]
	kernsz=size(kernpsl,/dim)
	
;Resize the kernel based on the scale conversion
	if min(kernsz/projpxscl_bpsep) lt 1 or finite(min(kernsz/projpxscl_bpsep)) ne 1 then kernpsl=fix(round(congrid(kernpsl,kernsz[0],kernsz[1]))) $
		else kernpsl=fix(round(congrid(kernpsl,kernsz[0]/projpxscl_bpsep,kernsz[1]/projpxscl_bpsep)))

	projmagg=ar_grow(projmag,rad=1,/gaus)
	psz=size(projmagg,/dim) 
	nmask=fltarr(psz[0],psz[1])
	pmask=nmask
	if (where(projmagg lt (-param.mdi_noisethresh*2)))[0] ne -1 then nmask[where(projmagg lt (-param.mdi_noisethresh*2))]=1
	if (where(projmagg gt param.mdi_noisethresh*2))[0] ne -1 then pmask[where(projmagg gt param.mdi_noisethresh*2)]=1
	pmaskg=ar_grow(pmask,mrad=1,inkernal=kernpsl)
	nmaskg=ar_grow(nmask,mrad=1,inkernal=kernpsl)
	wpsl=where(pmaskg+nmaskg eq 2)
	pslmask=fltarr(psz[0],psz[1])
	if wpsl[0] ne -1 then pslmask[wpsl]=1

	gradpsl=pslmask*ar_losgrad(projmagg)*projmmscl/ar_pxscale(magmap,/mmppx) ;
	
	pslmaskthresh=pslmask
	if (where(gradpsl lt param.psl_grad))[0] ne -1 then pslmaskthresh[where(gradpsl lt param.psl_grad)]=0

	pslmaskt=thin(pslmask)>0<1
	pslmaskt_thresh=thin(pslmaskthresh)>0<1

;Find the largest segment of PSL and indicate terminals
	pslmaskt_skel=thin(ar_largest_blob(pslmask,gradpsl),/neighbor)
	
;Find largest possible distance between 2 terminal points 
;	wterm=where(pslmaskt_long eq 2)
;	w2xy,wterm,psz[0],xterm,yterm

;!!!NEED to run separately for each PSL fragment?! then can sum all to get total PSL length
;Also, make a structure of all the fragment infos
;Might give better stable PSL lengths as AR rotates.

;Determine the longest PSLs Skeleton length and curvature
;	skelstr=ar_skeleton(pslmaskt_skel,max=2) ;set max eq 2 so that pixels linked by corners wont break the chain
;	if data_type(skelstr) eq 8 then begin
;		wskellong=(where(skelstr.LENGTH eq max(skelstr.LENGTH)))[0]
;		pslcurvature=(skelstr.curvature)[wskellong]
;	endif else pslcurvature=0.
;!!!!TEMP
pslcurvature=0.

;	plot_mag, projmag*projmask
;	contour,pslmaskt,level=0.5,/over



	meanmmscl=mean(projmmscl)

	if (where(pslmaskt eq 1))[0] eq -1 then begin
		print,'NO PSL'
		psllength=0.
	endif else begin
		psllendat=pslmaskt*projmmscl
		if (where(finite(psllendat) ne 1))[0] ne -1 then psllendat[where(finite(psllendat) ne 1)]=0
	
		psllength=total(psllendat)
	endelse

	if (where(pslmaskt_thresh eq 1))[0] eq -1 then begin
		print,'NO STRONG PSL'
		psllengtht=0.
	endif else begin
		psllendatt=pslmaskt_thresh*projmmscl
		if (where(finite(psllendatt) ne 1))[0] ne -1 then psllendatt[where(finite(psllendatt) ne 1)]=0
	
		psllengtht=total(psllendatt)
	endelse

;DETERMINE R
;The R kernel is NOT being rescaled depending on position
;Appears to give better results without that. 
;Because its a gausian convolution??

; compute pos and neg polarity maps, with product defining polarity inversion line:
	p1p=(smooth(float(rim gt  150),3) gt 0)
	p1n=(smooth(float(rim lt -150),3) gt 0)
	pmap=ar_r_smear(float(p1p*p1n), szkernel=10, param=param)>0.  ; 15Mm wide ridge, gaussian weighted

	rmap=(pmap*abs(rim))>1
	rmasked=rmap*projmask
	wbadr=where(finite(rmasked) ne 1)
	if wbadr[0] ne -1 then rmasked[wbadr]=0

	if (where(finite(rmasked) ne 1))[0] ne -1 then rmasked[where(finite(rmasked) ne 1)]=0

	thisr=total(rmasked)

;	print,'R=',thisr
;	print,'log(R)=',alog10(thisr)


;DETERMINE SUMMED GRADIENT (WLsg)

	wlsgdat=gradpsl*pslmask ;thresh
	if (where(finite(wlsgdat) ne 1))[0] ne -1 then wlsgdat[where(finite(wlsgdat) ne 1)]=0

	thiswlsg=total(wlsgdat)


;DETERMINE STRONG GRADIENT PSL LENGTH






;!!!!!!!! FIX THIS. 
;Really need to get the PSL skeleton and compare the skeleton length with the great circle length...
;	wpslthresh=where(pslmaskthresh)
;	xpslthresh=minmax(float(wpslthresh mod projsz[0]))
;	ypslthresh=minmax(float(wpslthresh/projsz[0]))
;	pslcurvature=psllength/sqrt((xpslthresh[1]-xpslthresh[0])^2.+(ypslthresh[1]-ypslthresh[0])^2.)/meanmmscl





;FILL STRUCTURE

;	thispslstr.datafile=arstr.datafile

	thispslstr.arid=i

	thispslstr.psllength=psllength
	
	thispslstr.pslsglength=psllengtht

	thispslstr.pslcurvature=pslcurvature
	
	thispslstr.rvalue=thisr

	thispslstr.wlsg=thiswlsg
	
	thispslstr.bipolesep_mm=bipsepstr.gcdist_mm
	thispslstr.bipolesep_px=bipsepstr.gcdist_px
	thispslstr.bipolesep_proj=bipsepstrproj.pxsep

arpslstr[i-1]=thispslstr

endfor

return, arpslstr

end