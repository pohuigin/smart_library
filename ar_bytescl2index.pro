;I accidentally used bytscl() instead of byte() when making the SMART masks
;Run this routine on the mask after reading from the fits file to make the AR mask values go from 0=bg, ar1=1,ar2=2, ... arN=N.

;mskdat=mapmsk.data
;mskdat=mskdat[sort(mskdat)]
;mskvals=mskdat[uniq(mskdat)]
;mapmsk.data=round((mapmsk.data-min(mapmsk.data))/(mskvals[1]-mskvals[0]))
function ar_bytescl2index, inmask

mask=inmask

mskdat=mask[sort(mask)]

mskvals=mskdat[uniq(mskdat)]

if n_elements(mskvals) eq 1 then return,fltarr((size(mask,/dim))[0],(size(mask,/dim))[1])

outmask=round((mask-min(mask))/(mskvals[1]-mskvals[0]))


return,outmask

end