;Generate a map of the solar disk that is 1 at disk center and goes radially outward as the cos(angle to LOS) 
;(= 2 at 60 degrees from LOS)
function ar_cosmap, map

wcs=fitshead2wcs(map)
coord=wcs_get_coord(wcs)
xx=reform(coord[0,*,*])
yy=reform(coord[1,*,*])
rr=(xx^(2.)+yy^(2.))^(0.5)
coscor=rr
coscor=1./cos(asin(coscor/map.rsun))
coscor[where(rr gt map.rsun)]=1.

cosmap=coscor
return, cosmap

end