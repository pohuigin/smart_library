;Generate a map of the solar disk that is 1 at disk center and goes radially outward as the cos(angle to LOS) 
;(= 2 at 60 degrees from LOS)
;optionally output:	rrdeg = gives degrees from disk center
;					wcs = wcs structure from input map file
;					offlimb = map of 1=on-disk and 0=off-disk
;					edgefudge = take off an extra half percent from the disk to get rid of limb effects
function ar_cosmap, map, rrdeg=rrdeg, wcs=wcs, offlimb=offlimb, edgefudge=edgefudge

if keyword_set(edgefudge) then fudge=0.999 else fudge=1.

wcs=fitshead2wcs(map)
coord=wcs_get_coord(wcs)
xx=reform(coord[0,*,*])
yy=reform(coord[1,*,*])
rr=(xx^(2.)+yy^(2.))^(0.5)
coscor=rr
rrdeg=asin(coscor/map.rsun)
coscor=1./cos(rrdeg)
coscor[where(rr gt map.rsun*fudge)]=1.

offlimb=rr
offlimb[where(rr ge map.rsun*fudge)]=0
offlimb[where(rr lt map.rsun*fudge)]=1

;stop

cosmap=coscor
return, cosmap

end