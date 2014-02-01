;provide detections of more limited extent.
;
;STATUS =  -1 - the initialised value in the routine
;			0 - skiped file because no detections were found by normal SMART run
;			1 - skipped file because no PSL detections were found
;			2 - no ridge/PSL mask overlay pixels found
;			3 - no strong pixels found
;			4 - no strong blobby pixels found
;			5 - no stong-psl mask overlap found (using region grow)
;			6 - no final core detections (failed multi-polarity test)
;			7 - Core detections were found (lucky number seven!)
;
;NOTES:
;	1. 	For each AR, id values go like 100,200,300
;	2.	All positions with SMART mask pixels have 1. added to the mask value allowing one to subtract 1 and remove all SMART masks
;	3. 	If the SMART mask pixels are touching a core detection the 2. is added, giving a value of 2. or 3.
;	4. 	Thus to isolate only the core detections and have them increment by 1 (instead of 100), one would do: core_mask = ceil( ( ( combined_mask - 3. ) > 0 ) / 100. )
;

function ar_detect_core, inmap,smartmask=insmartmask, nosmart=nosmart, $
	doprocess=doprocess, $
	params=inparams, fparam=fparam, doplot=doplot, $
	mapproc=mapproc, pslmaskmap=outpslmaskmap, status=status, $
	cosmap=cosmap,limbmask=limbmask,_extra=_extra

if keyword_set(doplot) then doplot=1 else doplot=0

if keyword_set(nosmart) then nosmart=1 else nosmart=0

;!!! need to cross calibrate with HMI to find the correct parameters to use

;!!! need to run yafta on the core detections to track them over time
;-keep track of which ARs are connected to each given AR at each point in time
;-keep track of when AR cores merge, split, etc...
;-should the tracking be done on the cores or the cores grown with the smart masks??


maporig=inmap
if not nosmart then smartmask=insmartmask

status=-1

if data_type(inparams) eq 8 then params=inparams $
	else params=ar_loadparam(fparam=fparam) ;get the default SMART parameter list
	
cmpmm=params.cmpmm ;cm per Mm

szorig=size(maporig.data,/dim)

;DO PROCESSING ON MAGNETOGRAM------------------------------------------------->
if not params.domedianfilt then nofilter=1

;Make a blank array for making masks, etc.
blank=fltarr(szorig[0],szorig[1])

if keyword_set(doprocess) then begin

;process/calibrate magnetogram
	mapproc=ar_processmag(maporig, _extra=_extra, cosmap=cosmap,limbmask=limbmask, params=params, fparam=fparam)

endif else mapproc=maporig

;Make a blank map structure
coremapmask=mapproc
coremapmask.data=blank
outpslmaskmap=coremapmask

;Initialise mask map
mapmsk=mapproc

if not nosmart then begin
	if max(smartmask) eq 0 then begin
		;No detections were found by normal SMART run
		status=0
		return,coremapmask
	endif

;Make a MAP structure from the SMART data.
	mapmsk.data=smartmask
;!!!TEMP -> turns out it was because I accidentally did bytscl() instead of byte()
;!!!For some reason	the MAPMSK doesn't have integer values 1,2,3...
;Add a note about this to the dataset I have online in the blog post...
	mapmsk.data=ar_bytescl2index(mapmsk.data)

endif else mapmsk.data=blank

rsgrad=params.smoothphys ;radius of SG in Mm
smoothhwhm=rsgrad*params.cmpmm/ar_pxscale(mapproc,/cmppx) ;the smoothing gaussian kernal HWHM

;Smooth the data (used for finding the PSL and PSL mask)
datasm=ar_grow(mapproc.data, /gaus, fwhm=smoothhwhm)

;Take abs(horiz_gradient)?
;	datagrad=ar_losgrad(datasm)

;Pull out the ridge of the gradient map
;	ridgemask=ar_ridgemask(datagrad)
ridgemask=ar_ridgemask(abs(datasm),thresh=params.smooththresh)

;Detect PSL by dilating positive and negative blobs

pslmask=ar_pslmask(datasm,radius=smoothhwhm,thresh=params.SMOOTHTHRESH,status=pslstatus)

if pslstatus ne 1 then begin
	;No PSLs were found
	status=1
	return,coremapmask
endif


;Overlay ridge detection and PSL to get better PSL tracing

wpsl=where(ridgemask+pslmask eq 2)

if wpsl[0] eq -1 then begin
	;no overlay PSL/ridge pixels found
	status=2
	return,coremapmask
endif

psltrace=blank
psltrace[wpsl]=1.

;Dilate the PSL trace

pslblobmask=ar_grow(psltrace,rad=smoothhwhm) ;/2.)

;Then combine PSL trace with normal AR detect masks?
;Could region grow normal mask with dilated PSL trace to group SMART detections that are fractured due to wide polarity separation

;Make strongfield masks
;Maybe take >200G pixels and dilate by SG radius?

wstrong=where(abs(mapproc.data) gt params.strongthresh)
;	wstrong=where(abs(datasm) gt params.smooththresh)

if wstrong[0] eq -1 then begin
	;no strong pixels found
	status=3
	return,coremapmask
endif

strongmask=blank
strongmask[wstrong]=1.

datastrongsm=ar_grow(abs(mapproc.data*strongmask), /gaus, fwhm=smoothhwhm)

;Determine the strong BLOB pixels

wstrongblob=where(abs(datastrongsm) gt params.smooththresh) ;!!!IS THIS THRESHOLD TOO HIGH??

if wstrongblob[0] eq -1 then begin
	;no strong blobby pixels found
	status=4
	return,coremapmask
endif

;	strongblobmask=ar_grow(strongmask,rad=smoothhwhm/2.)
;	strongblobmask=strongmask
strongblobmask=blank
strongblobmask[wstrongblob]=1.

;Trim the PSL blob array to where it overlaps with the strong blob array
pslblobmasktrim=strongblobmask*pslblobmask

wpsltrim=where(pslblobmasktrim eq 1)
if wpsltrim[0] eq -1 then begin
	;no overlap was found between strong blobs and the blobby PSL mask
	status=8
	return,coremapmask
endif	

;Also for actual core detection, overlay smoothed, strong-field thresholded detection with dilated PSL
;Use region grow with the PSL map as a base and the strong field map as a kernal

wgrow=region_grow(strongblobmask,wpsltrim,thresh=[0.5,1.5])

if wgrow[0] eq -1 then begin
	;no PSL blob mask -> strong blob mask overlap pixels found
	status=5
	return,coremapmask
endif	

arcoremask=pslblobmasktrim
arcoremask[wgrow]=1.

if doplot then begin
	erase
	!p.color=0
	!p.background=255
	!p.multi=[0,2,1]
	loadct,0,/sil
	plot_image,magscl(mapproc.data)
	setcolors,/sys,/sil
	contour,mapproc.data,level=[-200,200],c_color=[!green,!red],/over
	contour,datasm,level=[-20,20],c_color=[!green,!red],/over
	;contour,strongblobmask,level=0.5,c_color=!blue,/over
	contour,arcoremask,level=0.5,c_color=!cyan,c_thick=2,/over

endif

;Filter to only take multi polar detections -> if >95% of pixels (above strong field threshold) are of one polarity then disregard the detection.

arcoremaskid=label_region(arcoremask)

nid=max(arcoremaskid)
polfracarr=fltarr(nid)
for j=1,nid do begin
	wthiscore=where(arcoremaskid eq j)
	thiscoremask=blank
	thiscoremask[wthiscore]=1.
	wpos=where(strongblobmask*thiscoremask*mapproc.data gt params.strongthresh)
	wneg=where(strongblobmask*thiscoremask*mapproc.data lt -params.strongthresh)
	if wpos[0] eq -1 then npos=0 else npos=n_elements(wpos)
	if wneg[0] eq -1 then nneg=0 else nneg=n_elements(wneg)
	if npos eq 0 or nneg eq 0 then polfrac=1. $
		else polfrac=abs(npos-nneg)/float(npos+nneg)
	polfracarr[j-1]=polfrac
	if polfrac ge params.polethresh then arcoremaskid[wthiscore]=0.
endfor


;Final detection mask is original detections + core detections. 1=original, 2=core, 3=overlap
;need to write a routine to separate the two maps
;For each AR, id values go like 100,200,300
;All positions with SMART mask pixels have 1. added to the mask value allowing one to subtract 1 and remove all SMART masks
;If the SMART mask pixels are touching a core detection the 2. is added, giving a value of 2. or 3.
;Thus to isolate only the core detections and have them increment by 1 (instead of 100), one would do: core_mask = ceil( ( ( combined_mask - 3. ) > 0 ) / 100. )

;reset all values to 1
arcoremaskmpole=arcoremaskid < 1.

if max(arcoremaskmpole) eq 0 then begin
	;no final core detections (failed multi-polarity test)
	status=6
	return,coremapmask
endif	

;Re-index the final core mask
arcoremaskmpoleid=label_region(arcoremaskmpole)

;Make combined core and normal AR detection using region grow
;identifies regions that are attached to AR cores
smartmask=mapmsk.data < 1
arcoresmartcomb=arcoremaskmpole
wcoresmartgrow=region_grow(smartmask,where(arcoremaskmpole eq 1.),thresh=[0.5,1.5])
if wcoresmartgrow[0] ne -1 then arcoresmartcomb[wcoresmartgrow]=1.

;change the values to uniquely identify the ARs and blobs so that only 1 output mask needs to be saved
arcoremaskfinal=smartmask+arcoresmartcomb*2.+arcoremaskmpoleid*100.


if doplot then begin

	contour,arcoresmartcomb,level=0.5, c_color=!white,/over
	contour,arcoremaskmpoleid,level=0.5, c_color=!black,c_thick=2,/over,c_lines=2

	plot_image, arcoremaskfinal < 10.

	stop
endif


outmaskmap=coremapmask
outmaskmap.data=arcoremaskfinal

;Trim the thinned PSL mask
outpslmaskmap.data=strongblobmask*psltrace

if max(arcoremaskmpoleid) gt 0 then status=7

return, outmaskmap

end
