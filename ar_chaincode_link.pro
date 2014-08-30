;Make a x,y chain code, given a thinned mask

function ar_chaincode_link,chainmask,wstart,xmask,ymask, maxjump=inmaxjump, subsamp=subsamp, subchain=subchain

if n_elements(inmaxjump) ne 0 then maxjump=inmaxjump else maxjump=0

wthisstart=wstart
wbest=-1
maskother=chainmask
maskother[wstart]=0

chainarrx=xmask[wthisstart]
chainarry=ymask[wthisstart]

while wbest ne wstart do begin
	
	xthis=xmask[wthisstart]
	ythis=ymask[wthisstart]

	wother=where(maskother eq 1)

;Check to see if we've run out of links
	if wother[0] eq -1 then break

	xother=xmask[wother]
	yother=ymask[wother]

;Find the best match to the 
	dother=sqrt((xthis-xother)^2.+(ythis-yother)^2.)
	
	if maxjump ne 0 then $
		if min(dother) gt maxjump then break
	
	wwother=(where(dother eq min(dother)))[0]
	wbest=wother[wwother]

;Append the new link to the chain
	chainarrx=[chainarrx,xmask[wbest]]
	chainarry=[chainarry,ymask[wbest]]

;update the chain mask and zero the new reference pixel	
;	maskother=chainmask
	maskother[wbest]=0

;update the reference pixel	
	wthisstart=wbest

;loadct,0,/sile
;setcolors,/sys,/sil
;plot_image,chainmask
;plots,chainarrx,chainarry,ps=4,color=!cyan
;plots,xthis,ythis,ps=5,color=!red

;stop

endwhile

nchain=n_elements(chainarrx)

;Subsample the chain if requested
if n_elements(subsamp) ne 0 then begin
	
;require that subchain be at least 10 points
	if nchain ge 10.*subsamp then begin
		indarr=indgen(nchain)
		wgood=where(indarr mod subsamp eq 0)

		chainarrxs=[chainarrx[wgood],chainarrx[0]]
		chainarrys=[chainarry[wgood],chainarry[0]]
		subchain=[transpose(chainarrxs),transpose(chainarrys)]

	endif else begin
		subchain=[transpose([chainarrx,chainarrx[0]]),transpose([chainarry,chainarry[0]])]
	endelse
endif

chainarrx=[chainarrx,chainarrx[0]]
chainarry=[chainarry,chainarry[0]]



return,[transpose(chainarrx),transpose(chainarry)]

end


;----------------------------------------------------------------------------->