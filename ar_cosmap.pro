;Generate a map of the solar disk that is 1 at disk center and goes radially outward as the cos(angle to LOS) 
;(= 2 at 60 degrees from LOS)
;optionally output:	rrdeg = gives degrees from disk center
;					wcs = wcs structure from input map file
;					offlimb = map of 1=on-disk and 0=off-disk
;					edgefudge = take off an extra half percent from the disk to get rid of limb effects
function ar_cosmap, map, rrdeg=rrdeg, wcs=wcs, offlimb=offlimb, edgefudge=edgefudge, outcoord=outcoord

if n_elements(edgefudge) eq 1 then begin
   if edgefudge eq 1 then fudge=0.999
   if edgefudge ne 1 then fudge=edgefudge
endif else fudge=1.

wcs=map.wcs ;fitshead2wcs(map.index)
coord=wcs_get_coord(wcs)
xx=reform(coord[0,*,*])
yy=reform(coord[1,*,*])
rr=(xx^(2.)+yy^(2.))^(0.5)
coscor=rr
rrdeg=asin(coscor/map.rsun)
coscor=1./cos(rrdeg)
wgt=where(rr gt map.rsun*fudge)
if wgt[0] ne -1 then coscor[wgt]=1.

offlimb=rr
wgtrr=where(rr ge map.rsun*fudge)
if wgtrr[0] ne -1 then offlimb[wgtrr]=0
wltrr=where(rr lt map.rsun*fudge)
if wltrr[0] ne -1 then offlimb[wltrr]=1

;stop
outcoord=coord

cosmap=coscor
return, cosmap

end
