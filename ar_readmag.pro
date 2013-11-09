;Read a magnetogram fits file into a map with a complete header
function ar_readmag, fmag, inindex, mread=mread

;Check to see if a full header is bing read in to replace the dummy header in the FITS file.
if data_type(inindex) eq 8 then begin
	index=inindex
	doindex=1
endif else doindex=0

;Cheack if fits file has a header

;READ_SDO should work for both HMI and MDI data...
if keyword_set(mread) then mreadfits,fmag,ind,dat $
	else read_sdo,fmag,ind,dat

;Make a map with the full FITS header in addition to the map-specific keywords
if doindex then mindex2map, index, dat, map,/nest $
	else mindex2map, ind, dat, map,/nest

;Gets around a bug in DROT_MAP
;TEMP??
map=add_tag(map,map.time,'rtime',/no_copy)

outmap=map

return,outmap
end