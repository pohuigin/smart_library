;Read in a SMART2 meta structure (and optionally the corresponding AR meta structure) and filter the entries to get rid of files with missing/corrupted blocks.

function smart2_filter_meta, inmeta, armeta=inarmeta, outarmeta=armeta, $
	dotracking=dotracking, domagprop=domagprop
meta=inmeta



;#DATAFILE; LOCFILE; MASKFILE; DATE; TIM; TOT; MEAN; NAN; D00; V00; INTERVAL; ROLL_ANGLE; NAR; STATUS

;#DATAFILE; ARID; NARPX; ARPXPOSX; ARPXPOSY; ARHCPOSX; ARHCPOSY; ARPXWIDTHX; ARPXWIDTHY; ARHCWIDTHX; ARHCWIDTHY


window,xs=1400,ys=1400
!p.multi=[0,2,4]

tim=(meta.tim-anytim('1-jan-1990 00:00'))/(3600.*24.*365.)+90.

;plot,tim,abs(meta.mean),ytit='MEAN',chars=4,/ylog

;plot,tim,abs(meta.tot),ytit='TOT',chars=4,/ylog

;plot,tim,meta.nan,ytit='NAN',chars=4,/ylog

;plot,tim,meta.D00,ytit='D00',chars=4

;plot,tim,meta.V00,ytit='V00',chars=4

;plot,tim,meta.INTERVAL,ytit='INTERVAL',chars=4

;plot,tim,meta.STATUS,ytit='STATUS',chars=4

;plot,tim,meta.D00,ytit='D00',chars=4,yran=[2d5,4d5]

wgood=where(meta.d00 gt 2.6d5 and meta.d00 lt 3.3d5)
wbetter=where(meta[wgood].v00 gt -1d4 and meta[wgood].v00 lt 0)

;Probably need to filter for INTERVAL = 300 when doing magnetic property determination, 
;but for tracking, we want as many masks as possible.

wbest=where(meta[wgood[wbetter]].interval eq max(meta.interval))

;plot,tim[wgood[wbetter]],meta[wgood[wbetter]].D00,ytit='D00',chars=4

if keyword_set(domagprop) then wfilt=wgood[wbetter[wbest]] $
	else wfilt=wgood[wbetter]

;stop

meta=meta[wfilt]



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