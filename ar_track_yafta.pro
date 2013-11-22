;Run YAFTA for a range of times on a stack of SMART detection masks.


;----------------------------------------------------------------------------->
;initialise the state structure so that AR_TRACK_YAFTA can be continued in 
;subsequent runs for subsequent times.

;That won't work because it will be a different size each time!
;!! but maybe I should standardise the AR track struct so the tracking routines all use the same blank structures
;----------------------------------------------------------------------------->
;initialise the state structure so that AR_TRACK_YAFTA can be continued in 
;subsequent runs for subsequent times.
;
;function ar_track_yafta_state_init
;
;
;
;return,state
;
;end


;----------------------------------------------------------------------------->
;Run the tracking for a stack of masks
;
;Mandatory Keywords:
;	MAGSTACK = 
;	MASKSTACK = 
;	MDIMETA = 
;	SMARTMETA = 
;
;Optional Keywords:
;	PARAMS = 
;	OUTFILE = 
;
;In/Out Variable: STATE =  a structure that can be read out from and read back into 
;			AR_TRACK_YAFTA so the tracking can be continued over subsequent runs.
;NOTE: 		if a previous state is read-in then it must have been run on the 
;			same size masks/data that are read in also.
;
;Resturns: STRUCTURE =  The YAFTA tracking structure that corresponds to the input 
;			SMART AR detection meta info structure 
;

function ar_track_yafta,state,magstack=inmags,maskstack=insmmask, $
	mdimeta=inthissmart,smartmeta=inmdimetastrar, $
	params=inparams, fparam=fparam, $
	doplot=doplot, plotdir=plotdir, $
	outsingle=outsingle 	;out put a single map if there was only one on the day
								;so that it can be included with the next set and still be useful

if data_type(inparams) eq 8 then params=inparams $
	else params=ar_loadparam(fparam=fparam) ;get the default SMART parameter list

outsinglemap=''

help,state,/str

mags=inmags
smmask=insmmask
thissmart=inthissmart
mdimetastrar=inmdimetastrar

;Initialise a blank tracking structure
yaftastrblank=ar_struct_init(yaftaformat,structid='ar_track_yafta')

if data_type(state) eq 8 then begin
;This is for a continued run...

;Pull out the variables for the current YAFTA run

	mgram1 = state.mgram1
	ar_mask1 = state.ar_mask1
	ars1 = state.ars1
	ip1 = state.ip1 ;saves the last tracking 'step'

	if data_type(ars1) ne 8 then undefine,ars1
	
	tlastarsfound = state.tlastarsfound

;!!! Will need to add the following check into the loop 
;	 if delta time matching threshold is below 24 hrs

;Check to make sure the map was not from too long ago so that incorrect associations won't be made
	thisdeltatlastfound=thissmart[0].tim-tlastarsfound
	if thisdeltatlastfound ge params.tlastfoundthresh then begin
;Resets the tracking		
		undefine,ars1
	endif

	
	maxtrackid=state.maxtrackid
	
endif else begin

;This is for an initialised tracking run...
;Might need to make the state structure at the end...
;	state=ar_track_yafta_state_init(size())

;to guard against a crash when this variable is checked for, if there are no ARs found in an initial YAFTA run
	tlastarsfound=0.
	maxtrackid=0.
	ip1=0l
endelse



;Determine the size of input data

imgsz=size(mags,/dim)
nx=imgsz[0]
ny=imgsz[1]
if n_elements(imgsz) eq 2 then begin
	outsingle={mags:mags,smmask:smmask,thissmart:thissmart,mdimetastrar:mdimetastrar}
	print,'ONLY ONE FILE TO TRACK. ABORTING.'
	return,''
endif
nt=imgsz[2]



for i = 0,nt-1 do begin 

;Pull out the AR meta structure for this time
	wthisarmeta=where(mdimetastrar.datafile eq thissmart[i].datafile)
	if wthisarmeta[0] ne -1 then begin
		thisarmeta=mdimetastrar[wthisarmeta]
		nsmars=n_elements(wthisarmeta)
	endif else nsmars=0


	print, "Tracking step:",string(i+1) 

	mags_i = mags[*,*,i]
	ar_mask2 = smmask[*,*,i]

;FOR DEBUGGING!!!	
print,'SM MAX MASK = ',max(ar_mask2)
print,'SM NARS = ',n_elements(uniqpx(ar_mask2))-1
print,'NARS SM META = ',n_elements(thisarmeta)

	create_features, mags_i, ar_mask2, ars2, dx=params.yaftadx, min_size=params.yaftaminsize, peakthreshold=params.yaftapeakthresh


;If no features were found to track then keep track of how much time has elapsed since features were last found

if n_elements(ars2) eq 0 then begin

	if keyword_Set(doplot) then begin
		plot_image,magscl(mags_i),title='NO ARS FOUND!'
		contour,ar_mask2,level=0.5,color=0,/over
	endif

	thisdeltatlastfound=thissmart[i].tim-tlastarsfound

	if thisdeltatlastfound lt params.tlastfoundthresh then begin

;Pretends like the current map didn't happen and continues tracking
		continue
	endif else begin

;Resets the tracking		
		mgram1 = temporary(mags_i)
		ar_mask1 = temporary(ar_mask2)
		undefine,ars1

		continue
	endelse
endif else begin

;Save the time that ARs were last found
	tlastarsfound=thissmart[i].tim

	ip1=ip1+1l
	ars2.step = ip1
	
endelse

;Check if ARs existed previously and if they exist now
print,'line 131'

    n_ar1 = n_elements(ars1)

    n_ar2 = n_elements(ars2)

help,n_ar1,n_ar2


    possible_match = 0          ; default assumption

;If previous and current time have ARs, tracking will be attempted	
print,'line 138'

    if (n_ar1 gt 0) and (n_ar2 gt 0) then possible_match = 1


    if (possible_match eq 1) then begin

       step1 = mean(ars1.step)  ; it's clumsy to use "mean", but oh well

       step2 = mean(ars2.step)

       if (step2 - step1 ne 1) then possible_match = 0 ; not consecutive

;Update Max Label variable
		maxtrackid = maxtrackid > max(ars1.LABEL)

help,maxtrackid


    endif

;Attempt tracking
print,'line 148'

    if (possible_match eq 1) then begin

       match_features_v01, ars1, ars2, ar_mask1, ar_mask2, mgram1, mags_i, old_max_label = maxtrackid


       orig_ars2 = ars2

       orig_armask = ar_mask2

;Attempt to prevent fragmentation
print,'line 156'

       merge_fragments, ars2, ar_mask2

    endif                       ; ends if check for possible match

;Save the tracking YAFTA ID, ARID, data file name, and tim into a structure
	if (n_ar2 gt 0) then begin
print,'line 162'

if nsmars ne n_ar2 then print,'YAFTA NARS NE SM NARS!!!!'

		thisyaftastr=replicate(yaftastrblank,n_elements(thisarmeta))

		thisyaftastr.datafile=thisarmeta.datafile
		thisyaftastr.arid=thisarmeta.arid
		thisyaftastr.yaftaid=-1		
print,'line 170'

;Loop over ARS2 and find which tracked feature structure corresponds to which SMART detection
;AR detections with tracking (YAFTA) ID = -1 were not tracked/were disgarded by YAFTA for some reason

		for ys=0,n_elements(ars2)-1 do begin

;Check the value of the original mask at the first point of the tracked mask to determine what the original detection ID was.
			thissmid=(smmask[*,*,i])[(long(strsplit(ars2[ys].mask_str,/extract)))[0]]
			
			wsmid=where(thisyaftastr.arid eq thissmid)
			
			if wsmid[0] eq -1 then continue
			
			thisyaftastr[wsmid].yaftaid=ars2[ys].label
			thisyaftastr[wsmid].src=ars2[ys].src
			thisyaftastr[wsmid].trm=ars2[ys].trm
			thisyaftastr[wsmid].step=ars2[ys].step
		endfor

print,'line 185'
		if n_elements(thisyaftastrarr) ne 0 then thisyaftastrarr=[thisyaftastrarr,thisyaftastr] $
			else thisyaftastrarr=thisyaftastr

	endif	
    

;Concatenate current mask with array of previous masks
print,'line 192'


    if (n_elements(ar_masks) eq 0) then ar_masks = ar_mask2 $
	
    else ar_masks = [[[ar_masks]],[[ar_mask2]]]

;Over write previous arrays with current arrays for the next iteration
print,'line 198'


    mgram1 = temporary(mags_i)


    ar_mask1 = temporary(ar_mask2)


    if (n_ar2 gt 0) then begin
		ars1 = temporary(ars2)

;Concatenate current AR meta tracking info with array of previous ARs meta info
print,'line 208'

       if (n_elements(all_ars) eq 0) then  all_ars = ars1 $


       else all_ars = [all_ars, ars1]


;Plot the AR detections with the persistent tracking names
print,'line 215'

		if keyword_Set(doplot) then begin
			loadct,0,/sil
			display_yafta,magscl(mags[*,*,i]),/asp
			setcolors,/sys
			contour,smmask[*,*,i],level=0.5,color=0,thick=1,/over
			plot_edges, ar_mask1,thick=2
			plot_labels, ars1
		
			if n_elements(plotdir) eq 1 then window_capture,file=plotdir+thissmart[i].maskfile
		endif

	endif

endfor  ; ends loop over a given day's data

print,'line 230'


if n_elements(ars1) eq 0 then ars1=''

;Update the STATE variable and output it through the procedural call.

state={mgram1:mgram1,ar_mask1:ar_mask1,ars1:ars1,tlastarsfound:tlastarsfound,maxtrackid:maxtrackid,ip1:ip1}

help,state,/str

if n_elements(thisyaftastrarr) gt 0 then outstruct=thisyaftastrarr else outstruct=''

return,outstruct

end