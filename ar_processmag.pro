;Routine to process a LOS magnetogram for detecting ARs and determining magnetic properties.
;Runs:
;	Remove non-finites
;	Zero off-limb pixels 
;	Does Cosine correction
;	Do median filter
; 	Rotate solar north = up
function ar_processmag, inmap, limbmask=limbmask, cosmap=cosmap, nocosmicray=nocosmicray, $
	nocosine=nocosine, nofilter=nofilter, nofinite=nofinite, noofflimb=noofflimb, norotate=norotate, fparam=fparam, params=inparams

map=inmap
dat=map.data
imgsz=size(dat,/dim)

if data_type(inparams) eq 8 then param=inparams else param=ar_loadparam(fparam=fparam)


indextag=strlowcase(tag_names(map))

;Search for cosmic rays using hard threshold. Remove if gt 3sig detection than neighbooring pixels
if param.docosmicray and not keyword_set(nocosmicray) then begin	

	wcosmic=where(dat gt param.cosmicthresh)
	ncosmic=n_elements(wcosmic)
	print,'Cosmic Ray Candidates Found: ',ncosmic
	n=0
	if wcosmic[0] ne -1 then begin
		for i=0l,ncosmic-1l do begin
			wcx=wcosmic mod imgsz[0]
			wcy=wcosmic/imgsz[0]
			wneighboorx=wcx[0]+long([-1,0,1,1,1,0,-1,-1])
			wneighboory=wcy[0]+long([1,1,1,0,-1,-1,-1,0])
			if dat[wcosmic[i]] ge (3*stddev(dat[wneighboorx,wneighboory])+mean(dat[wneighboorx,wneighboory])) $
				then dat[wcosmic[i]]=mean(dat[wneighboorx,wneighboory])
		endfor
	endif
endif
	

;Clean NaNs
if not keyword_set(nofinite) then begin
	wnan=where(finite(dat) ne 1)
	if wnan[0] ne -1 then begin
		dat[wnan]=-9999.
		fill_missing, dat, -9999.,1
	endif
endif

;Get the cosine map and off-limb pixel map using WCS
cosmap=ar_cosmap(map, rrdeg=rrdeg, offlimb=offlimb,/edge)
limbmask=offlimb

;zero off-limb pixels
;zero from 80 degrees to LOS
if not keyword_set(noofflimb) then begin
	wofflimb=where(rrdeg/2./!pi*360. ge 80.)
	if wofflimb[0] ne -1 then dat[wofflimb]=0.
	dat=dat*offlimb
endif

;Median filter noisy values
;Do a 3x3 median filter
if param.domedianfilt and not keyword_set(nofilter) then begin
	medianw=param.medfiltwidth
	dat=filter_image(dat,median=medianw)
endif

;Do magnetic field cosine correction
;Limit correction to having 1 pixel at edge of the Sun
if not keyword_set(nocosine) then begin
	coscorlim=ar_coscorlim(map)
	cosmap=cosmap < coscorlim
	dat=dat*cosmap
endif

map.data=dat

;Rotate solar north = up
if not keyword_set(norotate) then begin
	map=rot_map(map,(-map.roll_angle))
	maptag=strlowcase(tag_names(map))
	if (where(indextag eq 'crota'))[0] ne -1 then map.index.crota=0.
	if (where(indextag eq 'crota2'))[0] ne -1 then map.index.crota2=0.
endif

return, map

end
