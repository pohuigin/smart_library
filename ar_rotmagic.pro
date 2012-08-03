;INDD - preprocessed data array 
;(defunct) OUTPLOT is filename for output PNG file
;
;Error Codes:
;	0 - Initialised value (no predicted error occurred)
;	1 - Data array is all 0's
;	2 - Detection is 100% unipolar
;	3 - Positive blob HG lon crossing detection unsuccessful 
;	4 - Negative blob HG lon crossing detection unsuccessful 
;	5 - No PSL between two main blobs detected
;	6 - No pixels in gradient map above the threshold
;	7 - No pixels in gradient map above the strong threshold

function ar_rotmagic,indd, rsundeg=rsundeg, map=inmap, argstr=argstr, genparams=genparams, $
	docrop=docrop, noderot=noderot, $
	_extra=_extra,debug=debug, silent=insilent, error=err ; arid=arid, plot=plot,outplot=outplot, outps=outps, ps=ps
err=0
silent=keyword_set(insilent)

dd=indd
map=inmap
map.data=dd
dx=map.dx
dy=map.dy
maptime=map.time

;if n_elements(arid) lt 1 then arid='01'

;if keyword_set(ps) then $
;	if n_elements(outps) lt 1 then outps='./'

if n_elements(argstr) gt 0 then begin
	gradthreshw=argstr.gradthreshw
	gradthreshs=argstr.gradthreshs
	smoothpole=argstr.smoothpole
	threshpole=argstr.threshpole

endif else begin
	if not silent then print,'% AR ROTMAGIC: No parameter structure input. Applying default values.'

	;Threshold for PSL detection
	gradthreshw=50.
	gradthreshs=150.
	
	;Threshold for polarity blob counting
	smoothpole=5. ;smoothing radius/FWHM for AR prior to polarity counting
	threshpole=.1 ;percent of total area a polarity has to be to be counted
endelse

;strrot={thetanl:10000d, thetasnl:10000d, thetabcl:10000d, lbcl:0d, npole:0d, polarity:0d, hglonlatpos:[0d,0d], hglonlatneg:[0d,0d]}
rotstr=ar_blanknar(/rotstr)

;for i=0,nfiles-1 do begin

imgsz=size(dd)
blankarr=fltarr(imgsz[1],imgsz[2])

;Do coordinate stuff
dxdy=[dx,dy]
xyrcoord,imgsz,xx,yy,rr
hglatlon_map, hglonmap,hglatmap, dxdy, imgsz, time=maptime, /mdi
;;hglonmap[where(hglonmap eq hglonmap[0,0])]=genparams.nan
;;hglatmap[where(hglatmap eq hglatmap[0,0])]=genparams.nan
;hgmissx=where(hglonmap eq hglonmap[0,0]) mod imgsz[1]
;hgmissy=where(hglonmap eq hglonmap[0,0])/imgsz[1]
;hglonmissval=interpolate(hglonmap,hgmissx,hgmissy)
;hglonmap[hgmissx,hgmissy]=hglonmissval
;FILL_MISSING, hglonmap, hglonmap[0,0], 2
;FILL_MISSING, hglonmap, hglonmap[0,0], 1
;FILL_MISSING, hglatmap, hglatmap[0,0], 2
;FILL_MISSING, hglatmap, hglatmap[0,0], 1

;Rotate data to disk center
mm=blankarr
wddn0=where(dd ne 0)
if wddn0[0] eq -1 then begin
	err=1
	goto,getout
endif
mm[wddn0]=1.

centall=[total(mm*abs(dd)*xx)/total(mm*abs(dd)),total(mm*abs(dd)*yy)/total(mm*abs(dd))]
hgcent=[hglonmap[centall[0],centall[1]],hglatmap[centall[0],centall[1]]]
map2=drot_map(map,-hgcent[0],/degrees,/rigid)

dd=map2.data

;stop

;Create polarity masks
if (where(dd ne 0))[0] eq -1 or (where(dd gt 0))[0] eq -1 or (where(dd lt 0))[0] eq -1 then begin
	err=2
	goto,getout
endif

mm=blankarr
mm[where(dd ne 0)]=1.
mp=blankarr
mp[where(dd gt 0)]=1
mn=blankarr
mn[where(dd lt 0)]=1

if keyword_set(docrop) then begin
	;Crop all arrays to the feature detection boundary
	dd=ar_crop(dd, mm, 1)
	;xx=smart_crop_ar(xx, mm, 1)
	;yy=smart_crop_ar(yy, mm, 1)
	;rr=smart_crop_ar(rr, mm, 1)
	hglonmap=ar_crop(hglonmap, mm, 1)
	hglatmap=ar_crop(hglatmap, mm, 1)
	blankarr=ar_crop(blankarr, mm, 1)
	mp=ar_crop(mp, mm, 1)
	mn=ar_crop(mn, mm, 1)
	mm=ar_crop(mm, mm, 1)
	
	;After crop, recreate coordinate arrays
	imgsz=size(dd)
	xyrcoord,imgsz,xx,yy,rr
endif



;Determine overall pos/neg centroids
rotstr.xypos=[total(mp*dd*xx)/total(mp*dd),total(mp*dd*yy)/total(mp*dd)]
rotstr.xyneg=[total(mn*dd*xx)/total(mn*dd),total(mn*dd*yy)/total(mn*dd)]

;Calculate HG centroid of positive and negative parts of detection
rotstr.hglonlatpos=[total(mp*hglonmap)/total(mp),total(mp*hglatmap)/total(mp)]
rotstr.hglonlatneg=[total(mn*hglonmap)/total(mn),total(mn*hglatmap)/total(mn)]

;test centroids TEMP!!!
;plot_image,mp-mn
;contour,hglonmap,level=[-10,0,10],color=0,/over
;contour,hglatmap,level=[-10,0,10],color=0,/over
;setcolors,/sys,/sil
;contour,hglonmap,level=-2.78170,color=!red,/over 
;contour,hglatmap,level=1.96685,color=!red,/over
;contour,hglonmap,level=2.47168,color=!green,/over  
;contour,hglatmap,level=1.12631,color=!green,/over 
;window_capture,file='largest_blobs'

;Calculate number of polarity blobs
smdd=ar_grow(dd,/gaus,rad=smoothpole)
smddp=blankarr
smddn=blankarr
if (where(smdd gt 0))[0] ne -1 then smddp[where(smdd gt 0)]=1.
if (where(smdd lt 0))[0] ne -1 then smddn[where(smdd lt 0)]=1.
smddp=LABEL_REGION(smddp)
smddn=LABEL_REGION(smddn)
if (where(smddn gt 0))[0] ne -1 then smddn[where(smddn gt 0)]=smddn[where(smddn gt 0)]+max(smddp)
smddcont=smddp+smddn
dummy=ar_largest_blob(smddcont,/nozero,narr=contnarr,/nosep)
wsm=where(contnarr ge threshpole*total(contnarr))
npole=n_elements(wsm)
rotstr.npole=npole

mp=ar_largest_blob(mp,dd,/flux)
mn=ar_largest_blob(mn,dd,/flux)

;Calculate HG centroid of largest positive and negative blobs
rotstr.hglonlatlgpos=[total(mp*hglonmap)/total(mp),total(mp*hglatmap)/total(mp)]
rotstr.hglonlatlgneg=[total(mn*hglonmap)/total(mn),total(mn*hglatmap)/total(mn)]

;test centroids TEMP!!!
;plot_image,mp-mn
;contour,hglonmap,level=[-10,0,10],color=0,/over
;contour,hglatmap,level=[-10,0,10],color=0,/over
;setcolors,/sys,/sil
;contour,hglonmap,level=3.2899742,color=!green,/over
;contour,hglatmap,level=1.5591141,color=!green,/over
;contour,hglonmap,level=-5.0236497,color=!red,/over 
;contour,hglatmap,level=3.8390343,color=!red,/over 
;window_capture,file='largest_blobs'

centall=[total(mm*abs(dd)*xx)/total(mm*abs(dd)),total(mm*abs(dd)*yy)/total(mm*abs(dd))]
centpos=[total(mp*dd*xx)/total(mp*dd),total(mp*dd*yy)/total(mp*dd)]
centneg=[total(mn*dd*xx)/total(mn*dd),total(mn*dd*yy)/total(mn*dd)]
;distpos=shift(rr,centpos[0]-imgsz[1]/2.,centpos[1]-imgsz[2]/2.)
;distneg=shift(rr,centneg[0]-imgsz[1]/2.,centneg[1]-imgsz[2]/2.)
rotstr.xylgpos=centpos
rotstr.xylgneg=centneg

;test centroids
;plot_image,mp-mn
;plots,[187.464,110.648,151.213],[114.990,133.560,118.015],ps=4,color=0
;plots,[187.464,110.648,151.213],[114.990,133.560,118.015],ps=2,color=255
;plots,[187.464,110.648,151.213],[114.990,133.560,118.015],ps=4,color=0

;Calculate bipole connection line length
lbcl=sqrt((centpos[0]-centneg[1])^2.+(centpos[1]-centneg[1])^2)*ar_pxscale(map,/mmppx)
rotstr.lbcl=lbcl

;hc_lb=[(((centall[0]-100) > 0)-imgsz[1]/2.)*dxdy[0],(((centall[1]-100) > 0)-imgsz[2]/2.)*dxdy[1]]
;hc_rt=[(((centall[0]+100) < 1023)-imgsz[1]/2.)*dxdy[0],(((centall[1]+100 < 1023))-imgsz[2]/2.)*dxdy[1]]

;hg_lb=conv_a2h(hc_lb,maptime)
;hg_rt=conv_a2h(hc_rt,maptime)

;Create contour mask of latitude line through centroid
hglatmask=blankarr
hglonposmask=blankarr
hglonnegmask=blankarr

;stop

contour,hglatmap,findgen(imgsz[1]),findgen(imgsz[2]),level=hglatmap[centall[0],centall[1]],path_info=path_info,/path_data_coords,path_xy=hglatcontxy,/over
hglatmask[hglatcontxy[0,*],hglatcontxy[1,*]]=1.
hglatmask=hglatmask*mm

contour,hglonmap,findgen(imgsz[1]),findgen(imgsz[2]),level=hglonmap[centpos[0],centpos[1]],path_info=path_info,/path_data_coords,path_xy=hglonposcontxy,/over
hglonposmask[hglonposcontxy[0,*],hglonposcontxy[1,*]]=1.
hglonposmask=hglonposmask*mm
hglonposmask=smart_grow(hglonposmask,radius=5)+hglatmask
if (where(hglonposmask ne 2))[0] eq -1 or (where(hglonposmask eq 2))[0] eq -1 then begin
	err=3
	goto,getout
endif
hglonposmask[where(hglonposmask ne 2)]=0 & hglonposmask[where(hglonposmask eq 2)]=1.

contour,hglonmap,findgen(imgsz[1]),findgen(imgsz[2]),level=hglonmap[centneg[0],centneg[1]],path_info=path_info,/path_data_coords,path_xy=hglonnegcontxy,/over
hglonnegmask[hglonnegcontxy[0,*],hglonnegcontxy[1,*]]=1.
hglonnegmask=hglonnegmask*mm
hglonnegmask=smart_grow(hglonnegmask,radius=5)+hglatmask
if (where(hglonnegmask ne 2))[0] eq -1 or (where(hglonnegmask eq 2))[0] eq -1 then begin
	err=4
	goto,getout
endif
hglonnegmask[where(hglonnegmask ne 2)]=0 & hglonnegmask[where(hglonnegmask eq 2)]=1.

;stop

;Calculate HG centroid of largest positive and negative blobs (X,Y pos??)
centhgpos=[total(hglonposmask*xx)/total(hglonposmask),total(hglonposmask*yy)/total(hglonposmask)]
centhgneg=[total(hglonnegmask*xx)/total(hglonnegmask),total(hglonnegmask*yy)/total(hglonnegmask)]

;Calculate angle between HG centroid arc and bipole connecting line 

;test for angle direction (inclined or declined to equator)
if centpos[0] gt centneg[0] then begin
	if centpos[1] gt centneg[1] then angsign=1. else angsign=-1.
endif else begin
	if centneg[1] gt centpos[1] then angsign=1. else angsign=-1.
endelse

alpha=angsign*vangle([(centpos[0]-centneg[0]),(centpos[1]-centneg[1]),0],[(centhgpos[0]-centhgneg[0]),(centhgpos[1]-centhgneg[1]),0])/!dtor

rotstr.thetabcl=alpha

;Find main PSL
pslmp=smart_grow(mp,rad=4)
pslmn=smart_grow(mn,rad=4)
pslmm=blankarr
if (where((pslmp+pslmn) eq 2.))[0] eq -1 then begin
	err=5
	goto,getout
endif
pslmm[where((pslmp+pslmn) eq 2.)]=1
pslmm=m_thin(pslmm)

;Create B field weak gradient mask
gradb=sqrt(deriv(dd)^2.+rot(deriv(rot(dd,-90)),90)^2.)/ar_pxscale(map, /mmppx)
gradbmaskw=blankarr
if (where(gradb ge gradthreshw))[0] eq -1 then begin
	err=6
	goto,getout
endif
gradbmaskw[where(gradb ge gradthreshw)]=1.
pslmmw=pslmm*gradbmaskw

;Find line through PSL and angle to HG arc
;weak psl
wpsl=where(pslmmw eq 1)
pslx=wpsl mod imgsz[1]
psly=wpsl/imgsz[1]
pslbm=linfit(pslx,psly)
pslx12=[min(pslx),max(pslx)]
psly12=pslbm[0]+pslbm[1]*pslx12
wpslalpha=vangle([pslx12[0]-pslx12[1], psly12[0]-psly12[1], 0],[(centhgpos[0]-centhgneg[0]),(centhgpos[1]-centhgneg[1]),0])/!dtor
rotstr.thetapsl=wpslalpha

rotstr.xy1psl=[pslx12[0],psly12[0]]
rotstr.xy2psl=[pslx12[1],psly12[1]]

;Create B field strong gradient mask
gradbmasks=blankarr
if (where(gradb ge gradthreshs))[0] eq -1 then begin
	err=7
	goto,skipstrong
endif
gradbmasks[where(gradb ge gradthreshs)]=1.
pslmms=pslmm*gradbmasks

;strong psl
wpsl=where(pslmms eq 1)
spslx=wpsl mod imgsz[1]
spsly=wpsl/imgsz[1]
spslbm=linfit(spslx,spsly)
spslx12=[min(spslx),max(spslx)]
spsly12=spslbm[0]+spslbm[1]*spslx12
spslalpha=vangle([spslx12[0]-spslx12[1], spsly12[0]-spsly12[1], 0],[(centhgpos[0]-centhgneg[0]),(centhgpos[1]-centhgneg[1]),0])/!dtor
rotstr.thetaspsl=spslalpha

rotstr.xy1spsl=[spslx12[0],spsly12[0]]
rotstr.xy2spsl=[spslx12[1],spsly12[1]]

skipstrong:

getout:

if keyword_set(debug) then begin
;if keyword_set(ps) then setplotenv,/ps,file=outplot
	loadct,0
	plot_image,dd,_extra=_extra

	if err eq 1 or err eq 2 or err eq 3 or err eq 4 then goto, skipplot
	xyouts,0.05,0.05,'theta_bcl='+strtrim(alpha,2),/norm
	setcolors,/sys
	oplot,[centhgpos[0],centhgneg[0]],[centhgpos[1],centhgneg[1]],ps=-4,color=!red
	oplot,[centpos[0],centneg[0]],[centpos[1],centneg[1]],ps=-1,color=!red
	oplot,hglatcontxy[0,*],hglatcontxy[1,*],lines=1
	oplot,hglonposcontxy[0,*],hglonposcontxy[1,*],lines=1
	oplot,hglonnegcontxy[0,*],hglonnegcontxy[1,*],lines=1

	if err eq 5 then goto, skipplot
	plots,pslx,psly,color=!cyan,ps=1
	oplot,pslx12,psly12,color=!blue,lines=0
;	contour,pslmm,level=.5,/over,color=!white
	xyouts,0.33,0.05,'theta_psl='+strtrim(wpslalpha,2),/norm

	if err eq 6 or err eq 7 then goto, skipplot
	plots,spslx,spsly,color=!orange,ps=1
	oplot,spslx12,spsly12,color=!yellow,lines=0
	xyouts,0.66,0.05,'theta_spsl='+strtrim(spslalpha,2),/norm

;	dumm=''
;	help,rotstr,/str
;	read,dumm
;if keyword_set(ps) then closeplotenv else $
;	window_capture,file=outplot,/png

print,'centpos',rotstr.xypos
print,'centneg',rotstr.xyneg
print,'lgcentpos',rotstr.xylgpos
print,'lgcentneg',rotstr.xylgneg
print,'psl1',rotstr.xy1psl
print,'psl2',rotstr.xy2psl
print,'spsl1',rotstr.xy1spsl
print,'spsl2',rotstr.xy2spsl

skipplot:
;stop

endif

;strrot={thetapsl:10000d, thetaspsl:10000d, thetabcl:10000d, lbcl:0d, npole:0d, polarity:0d, $ ;angle of weak neutral line, angle of strong neutral line, angle of bipole connecting line, length of bipole connecting line, number of polarity blobs, ???
;	hglonlatpos:[0d,0d], hglonlatneg:[0d,0d], hglonlatlgpos:[0d,0d], hglonlatlgneg:[0d,0d], $ ;overall pos/neg centroids, largest pos/neg blob centroids
;	xypos:[0d,0d],xyneg:[0d,0d],xylgpos:[0d,0d],xylgneg:[0d,0d], $ ;x/y pixel locations of overall and largest pos/neg blob centroids
;	xy1psl:[0d,0d],xy2psl:[0d,0d],xy1spsl:[0d,0d],xy2spsl:[0d,0d]}

;Derotate all coordinate properties so that they work with the original masks

;round coordinates so the decimals don't just get truncated
rotstr.xypos=round(rotstr.xypos) & rotstr.xyneg=round(rotstr.xyneg) & rotstr.xylgpos=round(rotstr.xylgpos) & rotstr.xylgneg=round(rotstr.xylgneg)
rotstr.xy1psl=round(rotstr.xy1psl) & rotstr.xy2psl=round(rotstr.xy2psl) & rotstr.xy1spsl=round(rotstr.xy1spsl) & rotstr.xy2spsl=round(rotstr.xy2spsl)

;stop

if not keyword_set(noderot) then begin
	if err ne 1 or err ne 2 then begin
		;de-rotate centroids to original centering of AR detection

		;DVD's method:
		;get_map_coord,map,xmap,ymap
		map_x = map
		map_y = map
		map_x.data = xx ;xmap
		map_y.data = yy ;ymap
		map2_x=drot_map(map_x,-hgcent[0],/degrees,/rigid)
		map2_y=drot_map(map_y,-hgcent[0],/degrees,/rigid)
		;plot_image,dd
		;cursor,x,y,/up
		;map2_x.data[x,y]
		;map2_y.data[x,y]

		rotstr.hglonlatpos=rotstr.hglonlatpos+[hgcent[0],0]
		rotstr.hglonlatneg=rotstr.hglonlatneg+[hgcent[0],0]
		rotstr.hglonlatlgpos=rotstr.hglonlatlgpos+[hgcent[0],0]
		rotstr.hglonlatlgneg=rotstr.hglonlatlgneg+[hgcent[0],0]

;		rotstr.xypos=hel2arcmin((rotstr.hglonlatpos)[1],(rotstr.hglonlatpos)[0],date=map.time)*60./[map.dx,map.dy]+[map.xc,map.yc]+[map.naxis1,map.naxis2]/2.
;		rotstr.xyneg=hel2arcmin((rotstr.hglonlatneg)[1],(rotstr.hglonlatneg)[0],date=map.time)*60./[map.dx,map.dy]+[map.xc,map.yc]+[map.naxis1,map.naxis2]/2.
;		rotstr.xylgpos=hel2arcmin((rotstr.hglonlatlgpos)[1],(rotstr.hglonlatlgpos)[0],date=map.time)*60./[map.dx,map.dy]+[map.xc,map.yc]+[map.naxis1,map.naxis2]/2.
;		rotstr.xylgneg=hel2arcmin((rotstr.hglonlatlgneg)[1],(rotstr.hglonlatlgneg)[0],date=map.time)*60./[map.dx,map.dy]+[map.xc,map.yc]+[map.naxis1,map.naxis2]/2.
		rotstr.xypos=[(map2_x.data)[(rotstr.xypos)[0],(rotstr.xypos)[1]], (map2_y.data)[(rotstr.xypos)[0],(rotstr.xypos)[1]]]
		rotstr.xyneg=[(map2_x.data)[(rotstr.xyneg)[0],(rotstr.xyneg)[1]], (map2_y.data)[(rotstr.xyneg)[0],(rotstr.xyneg)[1]]]
		rotstr.xylgpos=[(map2_x.data)[(rotstr.xylgpos)[0],(rotstr.xylgpos)[1]], (map2_y.data)[(rotstr.xylgpos)[0],(rotstr.xylgpos)[1]]]
		rotstr.xylgneg=[(map2_x.data)[(rotstr.xylgneg)[0],(rotstr.xylgneg)[1]], (map2_y.data)[(rotstr.xylgneg)[0],(rotstr.xylgneg)[1]]]
		
		if err lt 5 then begin
			;de-rotate PSL end positions 

;			psl1hc=(rotstr.xy1psl-[map.xc,map.yc]-[map.naxis1,map.naxis2]/2.)*[map.dx,map.dy]/60.
;			psl1hg=reverse(arcmin2hel(psl1hc[0],psl1hc[1],date=map.time))+[hgcent[0],0]
;			rotstr.xy1psl=hel2arcmin(psl1hg[1],psl1hg[0],date=map.time)*60./[map.dx,map.dy]+[map.xc,map.yc]+[map.naxis1,map.naxis2]/2.
			rotstr.xy1psl=[(map2_x.data)[(rotstr.xy1psl)[0],(rotstr.xy1psl)[1]], (map2_y.data)[(rotstr.xy1psl)[0],(rotstr.xy1psl)[1]]]

;			psl2hc=(rotstr.xy2psl-[map.xc,map.yc]-[map.naxis1,map.naxis2]/2.)*[map.dx,map.dy]/60.
;			psl2hg=reverse(arcmin2hel(psl2hc[0],psl2hc[1],date=map.time))+[hgcent[0],0]
;			rotstr.xy2psl=hel2arcmin(psl2hg[1],psl2hg[0],date=map.time)*60./[map.dx,map.dy]+[map.xc,map.yc]+[map.naxis1,map.naxis2]/2.
			rotstr.xy2psl=[(map2_x.data)[(rotstr.xy2psl)[0],(rotstr.xy2psl)[1]], (map2_y.data)[(rotstr.xy2psl)[0],(rotstr.xy2psl)[1]]]

			if err lt 6 then begin
;				spsl1hc=(rotstr.xy1spsl-[map.xc,map.yc]-[map.naxis1,map.naxis2]/2.)*[map.dx,map.dy]/60.
;				spsl1hg=reverse(arcmin2hel(spsl1hc[0],spsl1hc[1],date=map.time))+[hgcent[0],0]
;				rotstr.xy1spsl=hel2arcmin(spsl1hg[1],spsl1hg[0],date=map.time)*60./[map.dx,map.dy]+[map.xc,map.yc]+[map.naxis1,map.naxis2]/2.
				rotstr.xy1spsl=[(map2_x.data)[(rotstr.xy1spsl)[0],(rotstr.xy1spsl)[1]], (map2_y.data)[(rotstr.xy1spsl)[0],(rotstr.xy1spsl)[1]]]

;				spsl2hc=(rotstr.xy2spsl-[map.xc,map.yc]-[map.naxis1,map.naxis2]/2.)*[map.dx,map.dy]/60.
;				spsl2hg=reverse(arcmin2hel(spsl2hc[0],spsl2hc[1],date=map.time))+[hgcent[0],0]
;				rotstr.xy2spsl=hel2arcmin(spsl2hg[1],spsl2hg[0],date=map.time)*60./[map.dx,map.dy]+[map.xc,map.yc]+[map.naxis1,map.naxis2]/2.
				rotstr.xy2spsl=[(map2_x.data)[(rotstr.xy2spsl)[0],(rotstr.xy2spsl)[1]], (map2_y.data)[(rotstr.xy2spsl)[0],(rotstr.xy2spsl)[1]]]

			endif
		endif
	endif
endif

return,rotstr

end

;------------------------------------------------------------------------------>