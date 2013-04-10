;RUN AR ROTATION PROPERTIES ON GRID
;FMAG = magnetogram filename
;FSAV = SMART save file
;Error codes:
;	0 = no error
;	1 = missing user input (FMAG or FSAV)
;
;
;
;
pro run_grid_rotmagic, fmag=fmag, fsav=fsav, arpath=inarpath, arparam=inarparam, $
	err=err, debug=indebug, silent=insilent
err=0
debug=keyword_set(indebug)
silent=keyword_set(insilent)

;Define paths and files to run the code
if n_elements(inarpath) eq 1 then arpath=inarpath else arpath='~/science/procedures/smart_auxiliary/'
if n_elements(inarparam) eq 1 then arparam=inarparam else arparam='ar_param.txt'
ar_setup, ar_path=arpath, ar_param=arparam

;Load the parameter file- read in the parameters for running the code
defparams=ar_loadparam(silent=insilent)


;Read in data
if n_elements(fsav) ne 1 or n_elements(fmag) ne 1 then begin
	if not silent then print,'% RUN_GRID_ROTMAGIC: Must supply FMAG magnetogram filename, FSAV smart save file.'
	err=1 & return
endif
restore,fsav,/ver
map=ar_readmag(fmag) ;Make map with complete header (and rtime for drot_map bug)

;Form mask map
mmask=map
mmask.data=armask

;Run processing on magnetogram
map=ar_processmag(map, noisethresh=defparams.mdi_noisethresh) ;ar_argdefault(/mdi,/noisethresh))

;TEMP!!!
;window,xs=1000,ys=1000
if debug then loadct,0 & plot_image,armask & contour,armask,lev=0.5,c_col=255,/over
if debug then setcolors,/sys,/sil 
if debug then plots,arstruct.xpos,arstruct.ypos,color=!red,ps=4

;Set custom thresholds for running the code
;Thresholds for PSL detection, thresholds for polarity blob counting, smoothing radius/FWHM for AR prior to polarity counting
rotargstr={gradthreshw:50.,gradthreshs:150.,smoothpole:5.,threshpole:.1}

arid=string(findgen(max(armask))+1,form='(I02)')
for i=0,n_elements(arid)-1 do begin
	;isolate detection
	thismask=armask
	war=where(thismask ne float(arid[i]))
	thismask[war]=0
	thisdat=map.data*thismask
	;run rotation routine
	rotstr=ar_rotmagic(thisdat,map=map,argstr=rotargstr, genparams=defparams, debug=debug);,/noderot)

	if i eq 0 then rotarr=rotstr else rotarr=[rotarr,rotstr]
	if i eq 0 then roterrarr=err else roterrarr=[roterrarr,err]

endfor

if debug then print,'Detection error codes (ARID=err):'+strjoin(arid+'='+strtrim(roterrarr,2),', ')
;!temp!!!
if debug then plots,(rotarr.xypos)[0,*],(rotarr.xypos)[1,*],color=!cyan,ps=4

;Join (embed) the new rotation structure to the original and resave the file
rottag='rotstr' ;tag name for rotation structure in original one
arjoin=create_struct(arstruct[0],rottag,rotarr[0])
arjoin=replicate(arjoin,n_elements(arstruct))
;Fill new structure with the original data (upto but not including the last (new embedded one))
for i=0,n_elements(tag_names(arstruct))-1 do arjoin.(i)=arstruct.(i)
arjoin.(where(strlowcase(tag_names(arjoin)) eq strlowcase(rottag)))=rotarr ;fill the new embedded one with data

ARSTRUCT=arjoin
save,ARSTRUCT,EXTENTSTR,NOAASTR_DAILY,ARMASK,S_F,file=fsav+'.new'

;test output by plotting
;loadct,0 & plot_image,armask & contour,armask,lev=0.5,c_col=255,/over
;setcolors,/sys,/sil & oplot,arstruct.xpos,arstruct.ypos,color=!red,ps=4
;oplot,(rotarr.xypos)[0,*],(rotarr.xypos)[1,*],color=!yellow,ps=4

if debug then print,'DONE!!'
;stop

end
