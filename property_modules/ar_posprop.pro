;Input a processed (pre-run mag cos-corr) data map and mask
;The tot. area, pos. area, neg. area, and
;	total, signed, fractional signed, negative, and positive flux
;	are determined
;The returned structure contains position info for the whole
;detection.
;OUTPOSSTR has position information for the positive pixels in the
;detection
;OUTNEGSTR has position '' for negative ''
;
;Status:
;	0 = initialised value
;	7 = All desired positions should have been calculated
;	1 = No positive pixels found; positive positions not calulated/output
;	2 = No negative pixels found; negative positions not calulated/output
;	3 = No positive or nagative pixels found; Shouldn't happen!! (expect crash in this case)
;	4 = No detection 'where' values found for detection mask; Shouldn't happen!!
;	5 = No mask detection pixels found (max(mask)); Shouldn't happen unless you input a blank map!
;------------------------------------------------------------------------------>

;Determine the AR positions, given mask indices
;Input/Output:
;   inoutstr = a blank structure that will be filled
;Input:
;   WAR = an array of the AR pixel indices
;   DATA = the masked abs(data) array (cosine-area and -mag corrected; in flux units)
;

pro ar_posprop_findpos, inoutstr, war, indata, map=map
data=indata
data=abs(data)

date=map.time
dx=map.dx
dy=map.dy
rsun=map.rsun
xc=map.xc
yc=map.yc

sz=size(data,/dim)

xyrcoord,[0,sz],xx,yy

;X/Y indices
  xwar=war mod sz[0]
  ywar=war/sz[0]

;Determine the bounding box center positions
   xcenbnd=mean(float(minmax(xwar))) ;in pixels from lower left corner of the image FOV
   ycenbnd=mean(float(minmax(ywar))) ;in pixels

   px2hc, xcenbnd, ycenbnd, hcxbnd, hcybnd, dx=dx, dy=dy, xc=xc, yc=yc, xs=sz[0], ys=sz[1]

   hc2hg, hcxbnd, hcybnd, hgxbnd, hgybnd, carxbnd, date=date, rsunarcsec=rsun

;Determine the area weighted centroids
   xcenarea=total(xwar*xx[war])/total(xx[war])
   ycenarea=total(ywar*yy[war])/total(yy[war])

   px2hc, xcenarea, ycenarea, hcxarea, hcyarea, dx=dx, dy=dy, xc=xc, yc=yc, xs=sz[0], ys=sz[1]

   hc2hg, hcxarea, hcyarea, hgxarea, hgyarea, carxarea, date=date, rsunarcsec=rsun

;Determine the Flux weighted centroids
   xcenflx=total(xwar*data[war])/total(data[war])
   ycenflx=total(ywar*data[war])/total(data[war])

   px2hc, xcenflx, ycenflx, hcxflx, hcyflx, dx=dx, dy=dy, xc=xc, yc=yc, xs=sz[0], ys=sz[1]

   hc2hg, hcxflx, hcyflx, hgxflx, hgyflx, carxflx, date=date, rsunarcsec=rsun

;Fill structure
inoutstr.xcenbnd=xcenbnd
inoutstr.ycenbnd=ycenbnd
inoutstr.xcenflx=xcenflx
inoutstr.ycenflx=ycenflx
inoutstr.xcenarea=xcenarea
inoutstr.ycenarea=ycenarea
inoutstr.hcxbnd=hcxbnd
inoutstr.hcybnd=hcybnd
inoutstr.hcxflx=hcxflx
inoutstr.hcyflx=hcyflx
inoutstr.hcxarea=hcxarea
inoutstr.hcyarea=hcyarea
inoutstr.hglonbnd=hgxbnd
inoutstr.hglatbnd=hgybnd
inoutstr.hglonflx=hgxflx
inoutstr.hglatflx=hgyflx
inoutstr.hglonarea=hgxarea
inoutstr.hglatarea=hgyarea
inoutstr.carlonbnd=carxbnd
inoutstr.carlonflx=carxflx
inoutstr.carlonarea=carxarea

return

end


;------------------------------------------------------------------------------>


function ar_posprop, map=inmap, mask=inmask, cosmap=incosmap, params=inparams, $
                     outpos=outpos, outneg=outneg, nosigned=nosigned, status=status

status=0
mask=inmask
map=inmap

if data_type(inparams) eq 8 then params=inparams $
	else params=ar_loadparam() ;get the default SMART parameter list

if n_elements(incosmap) eq 0 then cosmap=ar_cosmap(map) $
	else cosmap=incosmap

pxmmsq=ar_pxscale(map,/mmsqr)
pxcmsq=ar_pxscale(map,/cmsqr)

blankstr={xcenbnd:0d, ycenbnd:0d, xcenflx:0d, ycenflx:0d, xcenarea:0d, ycenarea:0d, $
          hcxbnd:0d, hcybnd:0d, hcxflx:0d, hcyflx:0d, hcxarea:0d, hcyarea:0d, $
          hglonbnd:0d, hglatbnd:0d, hglonflx:0d, hglatflx:0d, hglonarea:0d, hglatarea:0d, $
          carlonbnd:0d,  carlonflx:0d,  carlonarea:0d}
arstr=blankstr
arposstr=blankstr
arnegstr=blankstr

if data_type(mask) ne 8 then begin
   maskstr=map & maskstr.data=mask & mask=maskstr
endif
nmask=max(mask.data)

if nmask eq 0 then begin
      status=5
      return,arstr
endif

;Make a status for each AR
status=fltarr(nmask)

strarr=replicate(blankstr,nmask)
strarrpos=strarr
strarrneg=strarr

for i=1,nmask do begin

;Zero pixels outside of detection boundary
   thismask=mask.data
   wzero=where(mask.data ne i)
   if wzero[0] ne -1 then $
      thismask[wzero]=0
   
   thisdat=map.data
   if wzero[0] ne -1 then $
      thisdat[wzero]=0

   thisabs=abs(thisdat)
   thisflx=thisabs*cosmap

;Where are values within the detection boundary
   wval=where(thismask eq i)
   if wval[0] eq -1 then begin
      status[i-1]=4
      continue
   endif

   thismask[wval]=1

   nothresh=0 & nopos=0 & noneg=0 & noposbnd=0 & nonegbnd=0

;Where are the signed values
   wneg=where(thisdat lt 0)
   wpos=where(thisdat gt 0)

;Fill the position structure for - whole detection
   ar_posprop_findpos, arstr, wval, thisflx, map=map

;The basic positions were calculated without incident
   status[i-1]=7

;Fill position structure
   strarr[i-1]=arstr

if not keyword_set(nosigned) then begin
   if wpos[0] eq -1 then arposstr=blankstr else $
;Fill the position structure for - Positive regions
      ar_posprop_findpos, arposstr, wpos, thisflx, map=map

   strarrpos[i-1]=arposstr

   if wneg[0] eq -1 then arnegstr=blankstr else $
;Fill the position structure for - Negative regions
      ar_posprop_findpos, arnegstr, wneg, thisflx, map=map

   strarrneg[i-1]=arnegstr

endif

	if wpos[0] eq -1 then status[i-1]=1
	if wneg[0] eq -1 then status[i-1]=2
	if wneg[0] eq -1 and wpos[0] eq -1 then status[i-1]=3
	
	strarr[i-1]=arstr

endfor

outstr=strarr
outpos=strarrpos
outneg=strarrneg

return,outstr

end
