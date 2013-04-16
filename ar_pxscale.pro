;Calculate the area of an magnetogram pixel at disk-centre on the solar surface.

function ar_pxscale, map, mmsqr=mmsqr, cmsqr=cmsqr, mmppx=mmppx, cmppx=cmppx, asecppx=asecppx, rsunasec=rsunasec, mdi=mdi

if n_elements(map) lt 1 then begin
	params=ar_loadparam()
	if keyword_set(mdi) then $
		map={dx:params.mdi_dx,dy:params.mdi_dy,rsun:params.mdi_rsun} $
		else begin & print,'% AR_MDIPXAREA: Must input MAP structure.' & return,'' & endelse
endif

;Calculate the area of a pixel in Mm^2
rsunmm=wcs_rsun(units='Mm')
;rsunmm=6.955d2 ;Mm
mmperarcsec=rsunmm/map.rsun ;Mm/arcsec
pixarea=((map.dx*mmperarcsec)*(map.dy*mmperarcsec)) ;Mm^2

if keyword_Set(cmsqr) then pixarea=pixarea*1D16

;Length of a side of a pixel.
arcppx=map.dx
arcpsun=map.rsun
mmpsun=rsunmm;6.955*1D8/1D6
retmmppx=(arcppx/arcpsun)*mmpsun ;in Mm/px

if keyword_set(mmppx) then return, retmmppx

if keyword_set(cmppx) then return, retmmppx*1d16

if keyword_set(asecppx) then return, map.dx

if keyword_set(rsunasec) then return, map.rsun

return, pixarea

end