;Routine to make a detection mask of ARs with mask values ordered from largest to smallest
;Assumes input mask has been:
;	-rotated to solar north up
;	-magnetic field values cosine corrected
;	-offlimb pixels zeroed
;DOPROCESS = set to do the above processing
;MAPPROC = Pull out the processed map
;REBIN1k = Do the detections on a magnetogram rebinned to 1kx1k

function,ar_detect, inmap, doprocess=doprocess, mapproc=mapproc, rebin1k=rebin1k
map=inmap

params=ar_loadparam() ;get the default SMART parameter list
cmpmm=params.cmpmm ;cm per Mm

;DO PROCESSING ON MAGNETOGRAM------------------------------------------------->
if keyword_set(doprocess) then begin
	map=ar_processmag, map, _extra=_extra
;		nocos=nocos, nofilter=nofilter, nofinite=nofinite, noofflimb=noofflimb, norotate=norotate
	mapproc=map
endif

data=map.data

if keyword_set(rebin1k) then $
	map1k=ar_rebin(map,/rebin1k)

;Gaussian smooth the magnetogram---------------------------------------------->
;Determine smoothing radius
rsgrad=params.smoothphys ;radius of SG in Mm
smoothhwhm=rsgrad*cmpmm/ar_pxscale(map1k,/cmppx) ;the smoothing gaussian kernal HWHM

;Smooth the data
datasm=ar_grow(data, /gaus, fwhm=smoothhwhm, kernal=kernal) ;, radius=radius

;Make a mask of detections--------------------------------------------------->
thresh=params.hmi_noisethresh














end