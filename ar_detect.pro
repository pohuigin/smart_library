;Routine to make a detection mask (map structure) of ARs with mask values ordered from largest to smallest
;Assumes input mask has been:
;	-rotated to solar north up
;	-magnetic field values cosine corrected
;	-offlimb pixels zeroed
;DOPROCESS = set to do the above processing
;MAPPROC = Pull out the processed map
;REBIN1k = Do the detections on a magnetogram rebinned to 1kx1k

function ar_detect, inmap, doprocess=doprocess, mapproc=mapproc, rebin4k21k=rebin4k21k, reduce=reducemap, $
	params=inparams, doplot=doplot
map=inmap

if keyword_set(params) then params=inparams $
	else params=ar_loadparam() ;get the default SMART parameter list
	
cmpmm=params.cmpmm ;cm per Mm

szorig=size(map.data,/dim)

;DO PROCESSING ON MAGNETOGRAM------------------------------------------------->
if keyword_set(doprocess) then begin
	map=ar_processmag(map, _extra=_extra)
;		nocos=nocos, nofilter=nofilter, nofinite=nofinite, noofflimb=noofflimb, norotate=norotate
	mapproc=map
endif

;Optionally set the ammount to reduce the size of the array by or
;perform the detections on a 4kx4k map reduced to 1kx1k map
if keyword_set(rebin4k21k) then reducemap=[4,4]
map=map_rebin(map,reduce=reducemap)

;stop

;Gaussian smooth the magnetogram---------------------------------------------->
;Determine smoothing radius
rsgrad=params.smoothphys ;radius of SG in Mm
smoothhwhm=rsgrad*cmpmm/ar_pxscale(map,/cmppx) ;the smoothing gaussian kernal HWHM

;Smooth the data
datasm=ar_grow(map.data, /gaus, fwhm=smoothhwhm)

;Make a mask of detections--------------------------------------------------->
smthresh=params.smooththresh

sz=size(map.data,/dim)
mask=fltarr(sz[0],sz[1])
mask[where(abs(datasm) gt smthresh)]=1.

;stop

;Segment the non-smoothed magnetogram to grab near by fragments and connect adjacent blobs 
magthresh=params.magthresh
growrad=smoothhwhm/2.

fragmask=fltarr(sz[0],sz[1])
fragmask[where(abs(map.data) ge magthresh)]=1
smfragmask=ar_grow(fragmask,rad=growrad)
poismask=where(mask eq 1)

;Region grow the smooth detections
wgrow=region_grow(smfragmask,poismask,thresh=[0.5,1.5])
grmask=mask
grmask[wgrow]=1

;Mask the offlimb pixels
dum=ar_cosmap(map, offlimb=limbmask, /edge)
grmask=grmask*limbmask

;stop

if keyword_set(doplot) then begin
	loadct,0
	plot_mag,map.data
	setcolors,/sys
	contour,mask,level=0.5,color=!blue,/over
	contour,grmask,level=0.5,color=!red,/over
stop
endif

;Resize the mask to the full resolution
maskfull=congrid(grmask,szorig[0],szorig[1])
maskfull[where(maskfull lt 0.5)]=0
maskfull[where(maskfull ge 0.5)]=1

;stop

;Separate the detections by assigning numbers
maskfull=label_region(maskfull)

;stop

;Order the detections by number of pixels
nar=max(maskfull)
arnpix=histogram(maskfull,bin=1,loc=arind)

maskorder=fltarr(szorig[0],szorig[1])

;stop

for i=0,nar-1 do begin
	rank=reverse(arind[sort(arnpix)])
	
	maskorder[where(maskfull eq rank[i])]=i+1.

endfor

;stop

maskmap=mapproc
maskmap.data=maskorder

if keyword_set(doplot) then begin
	loadct,39
	plot_map,maskmap,/limb
stop
endif

return,maskmap

end