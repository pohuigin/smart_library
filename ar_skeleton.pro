;Find the longest path along the (possibly) branched skeleton
;Input a mask with 2's at the terminals, and 3s in the intermediate points
;Allow jumps of max jump
;Set to xmask, ymask to lon, lat arrays if /spherical is set.

function ar_skeleton,inchainmask,inxmask,inymask, maxlump=inmaxjump,spherical=spherical, wstart=inwstart, status=status
chainmask=inchainmask
status=0

if n_elements(inmaxjump) ne 0 then maxjump=inmaxjump else maxjump=1

if keyword_set(spherical) then spherical=1 else spherical=0

imgsz=size(chainmask,/dim)

if n_elements(inxmask) lt 1 or n_elements(inymask) lt 1 then begin
	spherical=0
	xyrcoord,[0,imgsz],xmask,ymask

endif else begin
	xmask=inxmask
	ymask=inymask
endelse

if n_elements(inwstart) lt 1 then wstart=where(chainmask eq 2)
if wstart[0] eq -1 then begin
	print,'% AR_SKELETON: NO START POINTS FOUND! Returning.'
	status=-1
	return,''
endif

chainmask=round(inchainmask)>0<1

nstart=n_elements(wstart)

chainlenarr=fltarr(nstart)
chainstendist=fltarr(nstart)
chaincurvature=fltarr(nstart)

for i=0,nstart-1 do begin

	wthisstart=wstart[i]
	wbest=-1
	maskother=chainmask
	maskother[wthisstart]=0
	
	chainarrx=xmask[wthisstart]
	chainarry=ymask[wthisstart]
	
	while wbest ne wstart[i] do begin
		
		xthis=xmask[wthisstart]
		ythis=ymask[wthisstart]
	
		wother=where(maskother eq 1)
		nother=n_elements(wother)
	
;Check to see if we've run out of links
		if wother[0] eq -1 then begin
			break
		endif
	
		xother=xmask[wother]
		yother=ymask[wother]
	
;Find the best match to the current reference pixel
		if spherical then begin
			if nother gt 1 then dother=gc_dist([transpose(fltarr(nother)+xthis[0]),transpose(fltarr(nother)+ythis[0])],[transpose(xother),transpose(yother)]) $
				else dother=gc_dist([(xthis),(ythis)],[(xother),(yother)])
		endif else dother=sqrt((xthis[0]-xother)^2.+(ythis[0]-yother)^2.)
		
		
		if maxjump ne 0 then $
			if min(dother) gt maxjump then break

;update the total chain length array			
		chainlenarr[i]=chainlenarr[i]+min(dother)
		
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

;print,chainlenarr[i]
	
	endwhile

;only when making loops, which we are not, here	
;	chainarrx=[chainarrx,chainarrx[0]]
;	chainarry=[chainarry,chainarry[0]]
	
	if n_elements(chainarrx) gt 1 then status=execute('chain'+strtrim(fix(i),2)+'=[transpose(chainarrx),transpose(chainarry)]') $
		else status=execute('chain'+strtrim(fix(i),2)+'=transpose([(chainarrx),(chainarry)])')

;	status=execute('chstruct = CREATE_STRUCT(''ch''+strtrim(fix(i),2), thischain'+strtrim(fix(i),2)+', chstruct)') 


	if spherical then chainstendist[i]=gc_dist([chainarrx[0],chainarry[0]],[(reverse(chainarrx))[0],(reverse(chainarry))[0]]) $
		else chainstendist[i]=sqrt( ( chainarrx[0]-(reverse(chainarrx))[0] )^2.+( chainarry[0]-(reverse(chainarry))[0] )^2. )

endfor


chaincurvature=chainlenarr/chainstendist

status=execute('chstruct={'+strjoin('ch'+strtrim(indgen(nstart),2)+':chain'+strtrim(indgen(nstart),2),',')+'}')

chstruct = CREATE_STRUCT('length',chainlenarr,chstruct,'st_en_dist',chainstendist,'curvature',chaincurvature)



;help,chstruct,/str
;stop

return,chstruct

end


;----------------------------------------------------------------------------->