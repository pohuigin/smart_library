;Read a magnetogram fits file into a map with a complete header
function ar_readmag, fmag

mreadfits,fmag,ind,dat

;Make a map with the full FITS header in addition to the map-specific keywords
mindex2map, ind, dat, map,/nest

;Gets around a bug in DROT_MAP
;TEMP??
map=add_tag(map,map.time,'rtime',/no_copy)

outmap=map

return,outmap
end