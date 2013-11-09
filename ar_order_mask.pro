;Order the detections by number of pixels
;input a mask that has been indexed (ar1 =1, ar2=2, bg=0) using, for example, LABEL_REGION()

function ar_order_mask, inmask

mask=inmask
szorig=size(mask,/dim)

nar=max(mask)
arnpix=histogram(mask,bin=1,loc=arind)

maskorder=fltarr(szorig[0],szorig[1])

;stop

for i=0,nar-1 do begin
	rank=reverse(arind[sort(arnpix)])
	
	maskorder[where(mask eq rank[i])]=i

endfor

return,maskorder

end