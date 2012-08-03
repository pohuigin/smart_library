;The maximum cosine correction for a magnetogram
;this is the maximum factor of pixel area covered by a single pixel at the solar limb as compared with at disk centre
function ar_coscorlim, map, thetalim=thetalim

thetalim=asin(1.-map.dx/map.rsun)

corlim=1./cos(thetalim)

outcor=corlim
return, outcor

end