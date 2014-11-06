pro smart_eds, fmag=fmag, indmag=indmag, fparam=fparam, fhekparam=fhekparam, $
				write_fits=write_fits, $
				outroot=outroot, cadence=incadence, $ ; rebin1k=rebin1k, $
				events=events, error=error, imagerejected=imagerejected, $
				noevents=noevents, $
				debug=debug, status=status, detstr=detstr
;                       encode=encode

;----------------------------------------------------------------------------->
;+
;	PROCEDURE 
;			running AR detection algorithm on SDO HMI Mag LOS data
;
;	INPUT
;			fmag			list of mag files (45s or 720s, B LOS magnetograms)
;			fparam			file listing all detection parameters
;			fhekparam		file listing all parameters for making SMART VO events
;			indmag			(OPTIONAL) input index structure corresponding to FITS file
;							necessary if fits file only contains 'stub' header
;			status			(INPUT/OUTPUT) structure with info from previous runs. May be initialised and input as:
;							status={sszn:0l,imagerejected:0,noevents:0,detstatus:'',posstatus:'',magstatus:'',chainstatus:'', edserror:-1}
;	
;	OUTPUT
;			VO Events and (optionally) detection FITs files are written if >/= 1 detections are made.
;			detstr			(OPTIONAL) Structure array with detection properties for each AR detected
;			outfits			output FITs file name of detection mask
;
;	KEYWORDS
;			debug			Set if running the code manually (plots to Xwindow, etc.)
;			X (DEFUNCT) rebin1k			scaling factor (either Full 4k res. or reduced 1k res.)
;			X				0	Full resolution is used
;			X				1	the images are reduced to 1kx1k (SOHO resolution)
;			write_fits		write-out the byte (0-255) mask as a fits file
;			outroot			folder for writing output (fits and xml) files
;			cadence			time spacing between images in seconds (default is 6hr cadence)
;			events			output IDL structure of filled HEK events
;			error			output error code
;							-1	Initialised value (probably means crash early on)
;							0	SMART EDS ran without incident
;			imagerejected	-1	Initialised value
;							1	image was rejected
;							0	SMART2 was run on the code
;			noevents		-1	Initialised value
;							1	if no events were found
;							0	if events were found
;
;	NOTES
;			1. Keyword /noshell was used for READ_SDO to read-in the magnetogram 
;
;	MODIFICATIONS
;			2014-10-28	Written (P.A. Higgins)
;-
;----------------------------------------------------------------------------->

error=-1
imagerejected=-1
noevents=-1

;Check all inputs

if data_type(status) ne 8 or n_elements(status) ne 1 then message,'Input STATUS must be 1-element structure.'

if n_elements(cadence) eq 1 then cadence=incadence else cadence=6.*3600.

if data_type(fmag) ne 7 or n_elements(fmag) ne 1 then message,'Input FMAG must be 1-element string.'

if file_exist(fmag) ne 1 then message,'Magnetogram file specified by FMAG not found.'

if data_type(fparam) ne 7 or n_elements(fparam) ne 1 then message,'Keyword FPARAM must be 1-element string.'

if file_exist(fparam) ne 1 then message,'Parameter file specified by FPARAM not found.'

if data_type(fhekparam) ne 7 or n_elements(fhekparam) ne 1 then message,'Keyword FHEKPARAM must be 1-element string.'

if file_exist(fhekparam) ne 1 then message,'Parameter file specified by FHEKPARAM not found.'

;Load parameter files

params=ar_loadparam(fparam=fparam)

hekparam=ar_loadparam(fparam=fhekparam)

;Read in mag file

read_sdo,fmag,magind,magdat,/nosh

;Test mag file header and check for supplemental input header

if (where(tag_names(magind) eq 'date_obs'))[0] eq -1 then stubheader=1 else stubheader=0

if data_type(indmag) eq 8 and n_elements(indmag) eq 1 then suppheader=1 else suppheader=0

;If supplemental header is input then the FITS header is overwritten, whether or not it was a stub or full header.

if suppheader then magind=indmag

if stubheader and not suppheader then message,'FITS header is missing essential keywords. Must input full header structure (1-element structure variable) through INDMAG keyword to run SMART on this file.'
	
;Make a map structure from the index and data array

mindex2map,magind,magdat,magmap,/nest

magorig=magmap

;Process magnetogram and detect ARs

magproc=ar_processmag(magorig, cosmap=cosmap,limbmask=limbmask, params=params, /nofilter, /nocosmicray)

;Reduce to 1kx1k resolution
magprocrb=map_rebin(magproc,/rebin1k) ; else magprocrb=magproc



maskstrrb=ar_detect_core(magprocrb, /nosmart, params=params,doplot=debug, status=detstatus, cosmap=cosmap, limbmask=limbmask)
detmaskrb=ar_core2mask(maskstrrb.data)
maskstrrb.data=detmaskrb



;Expand mask to full res: 4k x 4k
maskstr=map_rebin(maskstrrb,xy=[4096,4096])
detmask=fix(round(maskstr.data))
maskstr.data=detmask

if keyword_set(write_fits) then begin

;Write fits file of full-res. mask as byte array

	outfits=outroot+'smart_mask_'+time2file(magind.date_obs)+'.fits'
	mwritefits, magind, byte(detmask), outfile=outfits
	spawn,'gzip -f '+outfits,/sh

endif

;Get number of ARs in image

nar=max(detmask)

if nar lt 1 then noevents=1 else noevents=0

;measure:
;	positions
;	mag
;	psl

if not noevents then begin



	posstr=ar_posprop(map=magproc, mask=detmask, cosmap=cosmap, params=params, status=posstatus)
    

    
	magpropstr=ar_magprop(map=magproc, mask=detmask, cosmap=cosmap, params=params, status=magstatus)
	areafrac=(magpropstr.posareabnd-magpropstr.negareabnd)/magpropstr.areabnd



	pslstr=ar_pslprop(magproc, detmask, param=params, doproj=0)
;	pslstr.psllength
;	pslstr.rvalue
;	pslstr.wlsg



;make chain code using the 

	dum=ar_chaincode(maskstrrb,detmaskrb,hekstructhc=chainstr,status=chainstatus, subsamp=hekparam.chainsubsamp, params=params)



;make cutout image of stereoscopic deprojection and mask fits????

;make HEK events


endif 


;Initialize data structure
strpropdet={sszn:0l,arid:0,datafile:'',date:'',hgpos:'',hcpos:'',pxpos:'', $
	PXSCL_HPC2STG:0.,DEG2DC:0.,NPSL:0, $
	AREA:0.,AREATHRESH:0.,BMAX:0.,FLUX:0.,FLUXFRAC:0., $
	BIPOLESEP:0d,AREAFRAC:0d,psllength:0d,PSLCURVATURE:0d,rvalue:0.,wlsg:0.,posstatus:-1,MAGSTATUS:-1,DETSTATUS:-1,chainstatus:-1,edserror:-1}
strpropdet=replicate(strpropdet,nar)

if not noevents then begin

;ARSTR
	strpropdet.SSZN=status.sszn+1l+lindgen(nar)          	;STRING    '0000001'
	strpropdet.DATAFILE=time2file(magind.date_obs)        	;STRING    'fd_M_96m_01d.1218.0000.fits'
	strpropdet.ARID=posstr.arid            					;INT              1
	strpropdet.DATE=anytim(magind.date_obs,/vms)            ;STRING    '2-May-1996 23:28:04.861'
	
	for i=0,nar-1 do begin
		strpropdet[i].HGPOS=strjoin(strtrim([string(posstr[i].HGLONFLX,form='(F15.5)'),string(posstr[i].HGLATFLX,form='(F15.5)')],2),',')           	;STRING    '-39.55054,-6.76069'
		strpropdet[i].HCPOS=strjoin(strtrim([string(posstr[i].HCXFLX,form='(F15.5)'),string(posstr[i].HCYFLX,form='(F15.5)')],2),',')      	;STRING    '-609.77055,-62.65432'
		strpropdet[i].PXPOS=strjoin(strtrim([string(posstr[i].XCENFLX,form='(F15.5)'),string(posstr[i].YCENFLX,form='(F15.5)')],2),',')    	;STRING    '203.71601,480.74146'
	endfor
	
	strpropdet.PXSCL_HPC2STG=-1   							;DOUBLE           0.0000000
	strpropdet.DEG2DC=-1          							;DOUBLE           40.026300
;	strpropdet.NPSL=pslstr.NARPSL          					;INT              1
	strpropdet.BMAX=abs(magpropstr.bmax) > abs(magpropstr.bmin)            						;FLOAT           1638.00
	strpropdet.AREA=magpropstr.AREABND            			;FLOAT           8797.80
	strpropdet.AREAFRAC=(magpropstr.POSAREA-abs(magpropstr.negAREA))/magpropstr.TOTAREA        	;DOUBLE          0.71740000
	strpropdet.AREATHRESH=magpropstr.totarea	      		;FLOAT           1514.40
	strpropdet.FLUX=magpropstr.totFLX            				;FLOAT       1.02050e+22
	strpropdet.FLUXFRAC=magpropstr.FRCFLX        			;DOUBLE          0.44610000
	strpropdet.BIPOLESEP=pslstr.bipolesep_mm       			;DOUBLE           30.269100
	strpropdet.PSLLENGTH=pslstr.PSLLENGTH       			;DOUBLE           24.172300
	strpropdet.PSLCURVATURE=pslstr.pslcurvature    			;DOUBLE           0.0000000
	strpropdet.RVALUE=pslstr.RVALUE          				;FLOAT           8879.90
	strpropdet.WLSG=pslstr.WLSG            					;FLOAT           9006.00
	strpropdet.POSSTATUS=posstatus       					;INT              7
	strpropdet.MAGSTATUS=magstatus       					;INT              0
	strpropdet.DETSTATUS=detstatus       					;INT              7
	strpropdet.chainstatus=chainstatus
	strpropdet.edserror=error       						;INT              1

endif

;Output detection structure

detstr=strpropdet


if not noevents then begin



	hekstr=ar_smart2hek(strpropdet, chain=chainstr, map=magproc, params=params, hekparam=hekparam)


	
	hekstr.required.EVENT_ENDTIME=anytim(anytim(magind.date_obs)+cadence,/ecs)
	
	for j=0,nar-1 do begin
	
		thisxml='smart_eds_voevent_'+strmid(time2file(magind.date_obs),0,13)+'_arid'+string(strpropdet[j].arid,form='(I02)')+'_prob'+strtrim(string(hekstr[j].optional.EVENT_PROBABILITY*10.,form='(I02)'),2)+'.xml'
			
		export_event,hekstr[j],/write,outfil=thisxml,outdir=outroot
	
	endfor

endif

;FUTURE addition: add deltaspot code

error=0




;Update the in/out status structure

status.sszn=max(strpropdet.SSZN)
status.imagerejected=imagerejected
status.noevents=noevents
status.posstatus=strjoin(posstatus,',')
status.magstatus=strjoin(magstatus,',')
status.detstatus=detstatus
status.chainstatus=strjoin(chainstatus,',')
status.edserror=error
































































stop

end