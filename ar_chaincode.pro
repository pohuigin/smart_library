;Out put a chain code for each AR in the mask.
;Uses HEK specifications
;STATUS:
;		0 = initialised value
;		1 = No ARs found in mask

;----------------------------------------------------------------------------->

function ar_chaincode,map,mask,arstruc=arstruc,params=params,dxdy=indxdy,xcyc=inxcyc, $
	hekstrout=hekstrout, hekstructhc=hekstructhc, hekstructpx=hekstructpx, status=status, aridsout=arids
status=0

;Check for data types and pull out origin and scaling
if data_type(map) ne 8 then begin
	mask=map
	dohc=0
endif else begin
	if n_elements(mask) eq 0 then mask=map.data
	dxdy=[map.dx,map.dy]
	xcyc=[map.xc,map.yc]
	dohc=1
endelse

;If origin and scale manually entered, then override map
if n_elements(indxdy) eq 2 and n_elements(inxcyc) eq 2 then begin
	dxdy=indxdy
	xcyc=inxcyc
	dohc=1	
endif

;obtain list of AR IDs from input info structure
if n_elements(arstruc) gt 0 then begin
	arids=arstruc.arid
	nar=n_elements(arids)
	if nar eq 1 and arids eq 0 then begin
		status=1
		return,''
	endif
endif else begin
	nar=max(mask)
	if nar eq 0 then begin
		status=1
		return,''
	endif
	arids=findgen(nar)+1
endelse


imgsz=size(mask,/dim)
blank=fltarr(imgsz[0],imgsz[1])

xyrcoord,[0,imgsz],xx,yy


;Initialise structure to load in chain codes

hekstructpx=ar_struct_init(structid='ar_chaincodehek')
hekstructpx=replicate(hekstructpx,nar)

hekstructhc=hekstructpx


;Loop over ARs
for i=0,nar-1 do begin

	thisarmask=blank
	warpx=where(fix(mask) eq fix(arids[i]))
	if warpx[0] eq -1 then continue
	thisarmask[warpx]=1.

;plot_image,thisarmask
;stop
	
	kernal=[[0,0,0,0,0],[0,1,1,1,0],[0,1,1,1,0],[0,1,1,1,0],[0,0,0,0,0]]
	gthisarmask=ar_grow(thisarmask,mrad=1,inkernal=kernal)
	
	chainmask=thin((gthisarmask-thisarmask)>0<1,/prune)<1

	wchain=where(chainmask eq 1)

	xchain=xx[wchain]
	ychain=yy[wchain]
	rchain=sqrt(xchain^2.+ychain^2.)
	
	wmin=wchain[(where(rchain eq min(rchain)))[0]]

;get the px chain code coordinates
	pxthischain=ar_chaincode_link(chainmask,wmin,xx*chainmask,yy*chainmask,maxjump=5)
	xthischain=pxthischain[0,*]
	ythischain=pxthischain[1,*]

if dohc then begin
;convert to HC
	hcxthischain=(xthischain-imgsz[0]/2.)*dxdy[0]-xcyc[0]
	hcythischain=(ythischain-imgsz[0]/2.)*dxdy[1]-xcyc[1]
	hcthischain=[transpose(reform(hcxthischain)),transpose(reform(hcythischain))]
endif

;Convert the chains to strings
	pxchainstring=strjoin(strtrim(fix(round(reform(pxthischain,n_elements(pxthischain)))),2),',')
	if dohc then hcchainstring=strjoin(strtrim(fix(round(reform(hcthischain,n_elements(hcthischain)))),2),',')


	hekstructpx[i].bound_ccnsteps= n_elements(xthischain)
	hekstructpx[i].bound_ccstartc1= fix(round(xthischain[0]))
	hekstructpx[i].bound_ccstartc2= fix(round(ythischain[0]))
	hekstructpx[i].bound_chaincode= pxchainstring
	hekstructpx[i].chaincodetype= 'ordered list'
	
	if dohc then begin
		hekstructhc[i].bound_ccnsteps= n_elements(hcxthischain)
		hekstructhc[i].bound_ccstartc1= fix(round(hcxthischain[0]))
		hekstructhc[i].bound_ccstartc2= fix(round(hcythischain[0]))
		hekstructhc[i].bound_chaincode= hcchainstring
		hekstructhc[i].chaincodetype= 'ordered list'
	endif

endfor




return,1

end