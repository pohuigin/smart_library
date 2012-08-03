;Crop around an isolated active region (everything else has already been zeroed)
;IMG = image to be cropped
;MASK = mask corresponding to feature detections
;ID = the detection ID to plot (0->99)
;STANDARD = Set to crop to the standard size (300x300)
;BOXDIM = Set a custom box size (XxY)
;BOXXYC = Set a custom box center (X,Y) 
;PLOT = Create a plot of the cropped image
;MMPLOT = Plot the axes in Mm. (Uses the MDI conversion Px->Mm)
;PSPLOT = Write the plot to a PS file
;OUTFILE = File to write the PS file to

function ar_crop, inimgin, inmaskin, id, arstruct=arstruct, zerother=zerother, standard=standard, $
	boxdim=bdim, boxxyc=boxxyc, plot=plot, mmplot=mmplot, asecplot=asecplot, $
	contourplot=contourplot, psplot=psplot, outfile=outfile, outmask=outmask, zoutmask=outmask0, $ ;, limbplot=limbplot
	minplot=minplot,maxplot=maxplot,_extra=_extra, outxxyy=outxxyy

imgin=inimgin & maskin=inmaskin

if n_elements(bdim) eq 2 then boxdim=bdim

if keyword_set(standard) then boxdim=[300,300]

if n_elements(id) ne 1 then id=1

if n_elements(arstruct) gt 0 then artitle=arstruct[id-1].smid+' '+arstruct[id-1].time else artitle=string(id,format='(I02)')

;if keyword_set(limbplot) then restore,smart_paths(/resmap,/no_calib)+'mdi_limbmask_map.sav'

;mmppx=smart_mdipxarea(/mmppx) ;mmsqr=mmsqr, cmsqr=cmsqr, mmppx=mmppx, cmppx=cmppx)
;asecppx=smart_mdipxarea(/asecppx)

img=imgin
mask=maskin

;Crop the data.
imgsz=size(img)

wzero=where(mask ne id)
mask[wzero]=0
if keyword_set(zerother) then img[wzero]=0
wn0x=where(total(mask,2) ne 0)
wn0y=where(total(mask,1) ne 0)
wx=[min(wn0x), max(wn0x)]
wy=[min(wn0y), max(wn0y)]
wxc=(wx[1]-wx[0])/2.+wx[0]
wyc=(wy[1]-wy[0])/2.+wy[0]

if n_elements(boxxyc) eq 2 then begin & wxc=boxxyc[0] & wyc=boxxyc[1] & endif

if n_elements(boxdim) eq 2 then begin
	bbx=boxdim[0]
	bby=boxdim[1]
	ximgrng=[wxc-bbx/2., wxc+bbx/2.]
	if ximgrng[0] lt 0 then ximgrng=[0,bbx]
	if ximgrng[1] gt (imgsz[1]-1) then ximgrng=[(imgsz[1]-1)-bbx[0],(imgsz[1]-1)]
	yimgrng=[wyc-bby/2., wyc+bby/2.]
	if yimgrng[0] lt 0 then yimgrng=[0,bbx]
	if yimgrng[1] gt (imgsz[2]-1) then yimgrng=[(imgsz[2]-1)-bbx,(imgsz[2]-1)]
	img=img[ximgrng[0]:ximgrng[1]-1,yimgrng[0]:yimgrng[1]-1]
	mask0=mask[ximgrng[0]:ximgrng[1]-1,yimgrng[0]:yimgrng[1]-1]
	mask=maskin[ximgrng[0]:ximgrng[1]-1,yimgrng[0]:yimgrng[1]-1]
	if keyword_set(limbplot) then limbmask=limbmask[ximgrng[0]:ximgrng[1]-1,yimgrng[0]:yimgrng[1]-1]

	xrng=[ximgrng[0],ximgrng[1]-1] & yrng=[yimgrng[0],yimgrng[1]-1]
endif else begin
	xrng=[wx[0]-50. > 0.,wx[1]+49. < 1023.] & yrng=[wy[0]-50. > 0.,wy[1]+49. < 1023.]

	img=img[xrng[0]:xrng[1], yrng[0]:yrng[1]]
	mask0=mask[xrng[0]:xrng[1], yrng[0]:yrng[1]] 
	mask=maskin[xrng[0]:xrng[1], yrng[0]:yrng[1]] 
;	if keyword_set(limbplot) then limbmask=limbmask[xrng[0]:xrng[1], yrng[0]:yrng[1]]

endelse

outmask0=mask0
outmask=mask
outxxyy=[xrng[0],yrng[0],xrng[1],yrng[1]]

crsz=size(img)

if keyword_set(plot) then begin
	loadct,0,/silent
	plotscale=1
	plotorigin=[(xrng[0])*plotscale,(yrng[0])*plotscale]
	xtitle='[pixels]' & ytitle='[pixels]'
	if keyword_set(mmplot) then begin
		plotscale=mmppx
		plotorigin=[-crsz[1]*plotscale/2.,-crsz[2]*plotscale/2.]
		xtitle='[Mm]' & ytitle='[Mm]'
	endif
	if keyword_set(asecplot) then begin
		plotscale=asecppx
		plotorigin=[(xrng[0]-imgsz[1]/2.)*plotscale,(yrng[0]-imgsz[2]/2.)*plotscale]
		xtitle='[arcseconds]' & ytitle='[arcseconds]'
	endif
	if not keyword_set(minplot) then minplot=-1200 & if not keyword_set(maxplot) then maxplot=1200
	plot_image,img,title=artitle,_extra=_extra, scale=plotscale, origin=plotorigin,min=minplot,max=maxplot, xtitle=xtitle,ytitle=ytitle
	if keyword_set(contourplot) then begin
		contx=findgen(crsz[1])*plotscale+plotorigin[0]
		conty=findgen(crsz[2])*plotscale+plotorigin[1]
		setcolors,/sys,/silent
		contour,mask,contx,conty,level=.5,c_color=!blue,/over
		contour,mask0,contx,conty,level=id-.5,c_color=!red,/over,c_thick=2
;		if keyword_set(limbplot) then contour,limbmask,contx,conty,level=.5,c_color=!black,/over
	endif
	loadct,0,/silent
endif

img_crop=img

return, img_crop

end