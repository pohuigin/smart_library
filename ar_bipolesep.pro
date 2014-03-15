;Determine the flux-weighted bipole separation distance between the pos and neg centroids 
;in degrees if a map is input, or in Px if only an image is input

function ar_bipolesep,indat

if data_type(indat) eq 8 then begin
	image=indat.data
	domap=1
endif else begin
	image=indat
	domap=0
endelse

imgsz=size(image,/dim)

xyrcoord,[0,imgsz],xx,yy

pxpxloc=total(xx*(image>0))/total(image>0)

nxpxloc=total(xx*abs(image<0))/total(abs(image<0))

pypxloc=total(yy*(image>0))/total(image>0)

nypxloc=total(yy*abs(image<0))/total(abs(image<0))

pxsep=sqrt((pxpxloc-nxpxloc)^2.+(pypxloc-nypxloc)^2.)

sepstr={pxcen:pxpxloc, pycen:pypxloc, nxcen:nxpxloc, nycen:nypxloc, plon:0., plat:0., nlon:0., nlat:0., pxsep:pxsep, gcdist_deg:0., gcdist_mm:0., gcdist_px:0.}

if domap then begin
	px2hc, pxpxloc, pypxloc, phcxflx, phcyflx, dx=indat.dx, dy=indat.dy, xc=indat.xc, yc=indat.yc, xs=imgsz[0], ys=imgsz[1]
	hc2hg, phcxflx, phcyflx, phgxflx, phgyflx, date=indat.time, rsunarcsec=indat.rsun

	px2hc, nxpxloc, nypxloc, nhcxflx, nhcyflx, dx=indat.dx, dy=indat.dy, xc=indat.xc, yc=indat.yc, xs=imgsz[0], ys=imgsz[1]
	hc2hg, nhcxflx, nhcyflx, nhgxflx, nhgyflx, date=indat.time, rsunarcsec=indat.rsun

	sepstr.plon=phgxflx
	sepstr.plat=phgyflx
	sepstr.nlon=nhgxflx
	sepstr.nlat=nhgyflx
	
	sepstr.gcdist_deg=gc_dist([phgxflx,phgyflx],[nhgxflx,nhgyflx])
	
	sepstr.gcdist_mm=sepstr.gcdist_deg/360.*2.*!pi*wcs_rsun(unit='Mm')
	
	sepstr.gcdist_px=sepstr.gcdist_deg/360.*2.*!pi*(indat.rsun/indat.dx)
	
endif

return,sepstr

end
