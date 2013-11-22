;Read in a SMART2 meta structure (and optionally the corresponding AR meta structure) and filter the entries to get rid of files with missing/corrupted blocks.
;Use a reference meta file to take only the 'data file' entries corresponding to those in the reference file.

function smart2_filter_meta_ref, inmeta, armeta=inarmeta, outarmeta=armeta, refmeta=inrefmeta, $
	dotracking=dotracking, domagprop=domagprop
meta=inmeta
refmeta=inrefmeta


;#DATAFILE; LOCFILE; MASKFILE; DATE; TIM; TOT; MEAN; NAN; D00; V00; INTERVAL; ROLL_ANGLE; NAR; STATUS

;#DATAFILE; ARID; NARPX; ARPXPOSX; ARPXPOSY; ARHCPOSX; ARHCPOSY; ARPXWIDTHX; ARPXWIDTHY; ARHCWIDTHX; ARHCWIDTHY



match,meta.datafile,refmeta.datafile,wmetaref,wrefmeta

meta=meta[wmetaref]


goodfiles=meta.datafile


if n_elements(inarmeta) ne 0 then begin
	armeta=inarmeta

;pull out unique data file names from AR meta file
	uardatafile=uniq(armeta.datafile)

	match,armeta[uardatafile].datafile,goodfiles,armatch,strmatch

;Make a list of where positions for good AR entries	

	ardatafile=armeta.datafile
	ardatafilenum=long(strmid(ardatafile,13,4)+strmid(ardatafile,18,4))
	
	datafilegood=goodfiles
	datafilegoodnum=long(strmid(datafilegood,13,4)+strmid(datafilegood,18,4))

	wargood=-1l
	for i=0l,n_elements(goodfiles)-1l do $
		wargood=[wargood,where(ardatafilenum eq datafilegoodnum[i])]
	
;		thissubran=[uardatafile[armatch[i]]-20,uardatafile[armatch[i]]]		
;		wargood=[wargood,long(thissubran[0]+where(ardatafile[thissubran[0]:thissubran[1]] eq datafilegood[i])]
;	endfor
	
	wfiltar=wargood[where(wargood ne -1)]
	
	armeta=armeta[wfiltar]
	
;	stop
	
	
endif




return, meta

end