pro run_ar_extract_template

;Set up parameters---------------------------------->

;Read in the AR detection parameter file
fparam='./ar_param.txt'
params=ar_loadparam(fparam=fparam)

;Input magnetogram
fmdi='./mag.fits'

;Output AR info CSV file
;fcsv=dpath+'smart_cutouts_metadata_alpha.txt'

;Run version of code (just to keep track of output CSV files)
;runvers='ALPHA.0'

;Cut-out size for AR zooms
;xycutsz=[600,600]
;determine from detection bounding box

;+/- Dynamic range for scaling the magnetograms 
magdisplay=1000


;Start detecting ARs-------------------------------->

;Read in a fits file (including WCS and full header)
thismap=ar_readmag(fmdi)
maporig=thismap

;filedate=anytim(file2time(fmdi),/vms)
fdate=time2file(thismap.time,/date)
fyyyy=strmid(fdate,0,4)

;Create AR mask (includes processing of MDI image -> read out into THISPROC)
thisarstr=ar_detect(thismap, /doprocess, mapproc=thisproc, params=params, /doplot, status=status, cosmap=cosmap, limbmask=limbmask)

;Overwrite original data with processed
thismap=thisproc

;Get number of ARs in image
nar=max(thisarstr.data)

;For each detection, pull out its X,Y position, make cut-out map (for plotting),
;determine mag properties

;Loop through each detection
for k=0l,nar-1l do begin

;Determine bounds of box surrounding AR detection
   imgsz=size(thisarstr.data,/dim)
   war=where(thisarstr.data eq k+1)
   wnotar=where(thisarstr.data ne k+1)

;Isolate the current AR of interest
   thisarmask=thisarstr.data
   thisarmask[wnotar]=0
   thisarmask[war]=1

;stop

;Determine bounding box in pixels
   ;x1,y1,x2,y2
   arbox=[min(war mod imgsz[0]), min(war/imgsz[0]), max(war mod imgsz[0]), max(war/imgsz[0])]
   ;get geometric centroid and width of box
   arxywidth=[arbox[2]-arbox[0], arbox[3]-arbox[1]]
   arxycent=[arbox[0]+arxywidth[0]/2., arbox[1]+arxywidth[1]/2.]

;Convert to image coordinates (arcsecs)
   hcarxywidth=arxywidth*thismap.dx
   hcarxycent=(arxycent-imgsz/2.)*thismap.dx+[thismap.xc,thismap.yc]


;make cutout image   
   narmap=map_cutout(thismap, xycen=hcarxycent, xwidth=hcarxywidth[0], yheight=hcarxywidth[1], auxdat1=thisarstr.data, outauxdat1=narmask, auxdat2=thisarmask, outauxdat2=narthismask, auxdat3=limbmask, outauxdat3=narlimb)
      
   narsz=size(narmap.data,/dim)
   
;Isolate the AR within the FD magnetogram
   thisarmap=thismap
   thisarmap.data=(thisarmap.data)*thisarmask
   
;Determines magnetic properties of the AR
;-> need to add x,y position, neg./pos. centroid in HG coord
   magpropstr=ar_magprop(map=thisarmap, mask=thisarmask, cosmap=cosmap, params=params)

;Extra property to determine=polarity imbalanced area coverage
;-> need to add this to AR_MAGPROP
   areafrac=(magpropstr.posareabnd-magpropstr.negareabnd)/magpropstr.areabnd
      

;Examine the output for each AR---------------------------------->

;Display the AR cutout
   loadct,0
   plot_image,magscl(narmap.data)
   setcolors,/sys
   contour,narlimb,level=0.5,/over,color=!white,c_thick=2
   contour,narmask,level=0.5,/over,color=!blue,c_thick=2
   contour,narthismask,level=0.5,/over,color=!red,c_thick=2

;Output AR data/detection/property structure
   help,thismap,/str
   help,thismap.index
   help,thismap.wcs
   help,magpropstr,/str

stop

endfor

stop

end
