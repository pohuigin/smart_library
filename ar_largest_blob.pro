;Input MASK with 1=feature, 0=quiet
;Returns MASK with all features zeroed except for the largest
;set flux to take the flux weighted largest blob

function ar_largest_blob, inmask, indata, flux=flux,nozero=nozero,narr=narr, nosep=nosep, _extra=_extra

mask=inmask
if keyword_set(flux) gt 0 then data=indata

if keyword_set(nosep) then masksep=mask else masksep=LABEL_REGION(mask, _extra=_extra)

ncont=max(masksep)

if ncont lt 1 then return,mask

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