;Input MASK with 1=feature, 0=quiet
;Returns MASK with all features zeroed except for the largest

function ar_largest_blob, inmask, indata, flux=flux,nozero=nozero,narr=narr, nosep=nosep

mask=inmask
if keyword_set(flux) gt 0 then data=indata

if keyword_set(nosep) then masksep=mask else masksep=LABEL_REGION(mask)

ncont=max(masksep)
narr=fltarr(ncont)
i=1
while i le ncont do begin
	if n_elements(flux) gt 0 then narr[i-1]=total(abs(data[where(masksep eq i)])) $
		else narr[i-1]=n_elements(where(masksep eq i))
	i++
endwhile

if keyword_set(nozero) then return,''

wnbest=(where(narr eq max(narr)))[0]+1.
wbig=(where(masksep eq wnbest))
w0=(where(masksep ne wnbest))
mask[w0]=0

return,mask

end