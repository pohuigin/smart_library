;Input the output mask from AR_DETECT_CORE
;
;Convert the 3-layered core mask to:
;	-an indexed core mask going 0=bg,1,2,3,4...
;	-a SMART mask of 0's and 1's
;	-a SMART mask with only the blobs that are touching cores

function ar_core2mask, inmask,smartmask=outsmartmask,coresmartmask=outsmartconn
mask=inmask

core_mask = ceil( ( ( mask - 3. ) > 0 ) / 100. )


smartmask=(mask-core_mask*100.) < 1

smartconn=(mask-core_mask*100.) > 1 < 2


outsmartmask=ceil(smartmask)

outsmartconn=ceil(smartconn)

return, core_mask

end