;Input an euv map and smart mask
;cosine correction for magnetic field values should already be done
;The tot. area, pos. area, neg. area, and
;	total, signed, fractional signed, negative, and positive flux
;	are determined
;
;Mask should correspond to AIA map
;
;Properties:
;		...intense: the max(), mean(), or total() of the intensity values in the detection for a given EUV image
;		loopmedianshear: the median angle of overlying loop-orientation along a PSL segments. 90 degrees means NO SHEAR, 0 degrees means MAXIMAL SHEAR 
;		loopnpsl: the number of psl segments in the detection
;
;STATUS = output status array keyword (one status for each AR)
;		0: initialised value
;		7: Everything went swimmingly, and there should be valid magnetic properties for each AR
;		1: No ARs were present in the mask! Should only happen if a blank array was read in

function ar_euvprop, map=inmap, mask=inmask, pslmask=pslmask, params=inparams, fparam=fparam, status=status, datafile=indatafile, outpslangle=pslangmap, outloopangmap=loopangmap, outplotstr=outplotstr
status=0

if data_type(inparams) eq 8 then params=inparams $
	else params=ar_loadparam(fparam=fparam) ;get the default SMART parameter list

mask=inmask
map=inmap

pxmmsq=ar_pxscale(map,/mmsqr)
pxcmsq=ar_pxscale(map,/cmsqr)

if n_elements(indatafile) ne 0 then begin
	blankstr={datafile:indatafile[0], arid:0, maxintense:0d, minintense:0d, meanintense:0d, medianintense:0d, totalintense:0., loopmedianshear:0., loopnpsl:0, loopmedianshearlrg:0., looptotpsl:0., looptotpsllrg:0.}
endif else blankstr={arid:0, maxintense:0d, minintense:0d, meanintense:0d, medianintense:0d, totalintense:0., loopmedianshear:0., loopnpsl:0, loopmedianshearlrg:0., looptotpsl:0., looptotpsllrg:0.}

if data_type(mask) ne 8 then begin
   maskstr=map & maskstr.data=mask & mask=maskstr
endif
nmask=max(mask.data)

imgsz=size(mask.data,/dim)

;Check that there are ARs present in the mask
if nmask eq 0 then begin 
	status=1
	return,blankstr
endif

strarr=replicate(blankstr,nmask)

for i=1,nmask do begin

	strarr[i-1].arid=i

;Zero pixels outside of detection boundary
   thismask=mask.data
   wzero=where(mask.data ne i)
   if wzero[0] ne -1 then $
      thismask[wzero]=0
   
   thisdat=map.data
   if wzero[0] ne -1 then $
      thisdat[wzero]=0

   thisabs=abs(thisdat)

;Where are values within the detection boundary
   wval=where(thismask eq i)
   if wval[0] eq -1 then continue
   
   thismask[wval]=1

   nothresh=0 & nopos=0 & noneg=0 & noposbnd=0 & nonegbnd=0

;stop

   strarr[i-1].maxintense=max(thisdat)

   strarr[i-1].minintense=min(thisdat)

   strarr[i-1].meanintense=mean(thisdat)

   strarr[i-1].medianintense=median(thisdat)

   strarr[i-1].totalintense=total(thisdat)


if n_elements(pslmask) ne 0 then begin

	if max(pslmask) eq 1 then begin


;stop

;Crop around the AR mask to save computation time
cutbuffer=10
xrancrop=minmax(where(thismask eq 1) mod imgsz[0])+[-cutbuffer,cutbuffer]
yrancrop=minmax(where(thismask eq 1)/imgsz[0])+[-cutbuffer,cutbuffer]

thismaskcr=thismask[xrancrop[0]:xrancrop[1],yrancrop[0]:yrancrop[1]]
mapdatcr=(map.data)[xrancrop[0]:xrancrop[1],yrancrop[0]:yrancrop[1]]
pslmaskcr=pslmask[xrancrop[0]:xrancrop[1],yrancrop[0]:yrancrop[1]]
thisdatcr=thisdat[xrancrop[0]:xrancrop[1],yrancrop[0]:yrancrop[1]]

	   thispslmask=(pslmaskcr*thismaskcr)>0.<1.
	   grownpsl=thisdatcr*ar_grow(thispslmask,rad=10,/gaus)
	
	   loopangmap=directional_filter(mapdatcr, $
	   		hwid=params.dirhwid, r0=params.dirr0, l0=params.dirl0, b0=params.dirb0, $
	   		binfac=1, arrowbin=1, outmag=loopmagnitude, outplotstr=outplotstr)

;	   loopangarr=loopangmap[where(thispslmask eq 1)]/!dtor/2.
;	   pslangmap=ar_pslangmap(thispslmask)



pslangmap=directional_filter(thispslmask, $
       hwid=params.pslhwid, r0=params.dirr0, l0=params.dirl0, b0=params.dirb0, $
       binfac=1, arrowbin=1, outmag=pslmagnitude)

angdiffmap=abs(loopangmap-pslangmap)/!dtor/2.
angdiffarr=angdiffmap[where(thispslmask eq 1)]

angweights=loopmagnitude[where(thispslmask eq 1)] < 1.

;Take the longest PSL
thispslmasklrg=ar_largest_blob(thispslmask,/all_neighbors)
angdiffarrlrg=angdiffmap[where(thispslmasklrg eq 1)]

angweightslrg=loopmagnitude[where(thispslmasklrg eq 1)] < 1.

strarr[i-1].loopmedianshear=total(angdiffarr*angweights)/total(angweights)

strarr[i-1].loopmedianshearlrg=total(angdiffarrlrg*angweightslrg)/total(angweightslrg)

strarr[i-1].loopnpsl=max(label_region(thispslmask,/all))

strarr[i-1].looptotpsl=total(thispslmask)

strarr[i-1].looptotpsllrg=total(thispslmasklrg)


;plot_hist,angdiffmap[where(thispslmask eq 1)]

;JUST USE DIRECTION FINDING ON PSL
;NEED TO dECREASE LOBES SO<EHOW
;THEN SUBTRACT ANGLE AND TAKE PEAK(???) OF HISTOGRAM AS THE SHEAR ANGLE
;90=no shear
;loop over J and watch the shear angle

	endif

endif





   strarr[i-1].arid=i

endfor


outstr=strarr





return,outstr

end
