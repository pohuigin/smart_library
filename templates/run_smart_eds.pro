;wrapper for smart_eds.pro

pro run_smart_eds

ssw_jsoc_time2data,dum1,dum2,drmsmag,filesmag,/files_only,cadence='3600s',ds='hmi.M_720s_nrt',/jsoc2,lastn=6.

fmag=filesmag[0]
indmag=drmsmag[0]
fparam='./ar_param_hmi_eds.txt'
fhekparam='./ar_param_hek_eds.txt'
write_fits=1
;rebin1k=1
outroot='./'
cadence=3600.*6.
debug=0

status={sszn:0l,imagerejected:0,noevents:0,detstatus:0,posstatus:0,magstatus:0,chainstatus:0,edserror:-1, lastmag:fltarr(4096,4096)-1., lastmask:fltarr(4096,4096)-1., trackstatus:ar_struct_init(yaftaformat,structid='ar_track_yafta', fstruct='./ar_struct_param.txt'),trackmetasm:{},trackmetadat:{}}

timst=anytim(systim(/utc))

smart_eds, fmag=fmag, indmag=indmag, fparam=fparam, fhekparam=fhekparam, $
				write_fits=write_fits, $ ;rebin1k=rebin1k, $
				outroot=outroot, cadence=cadence, $
				events=events, error=error, imagerejected=imagerejected, $
				noevents=noevents, $
				debug=debug, status=status

timen=anytim(systim(/utc))

print,'Run time = ',timen-timst,'seconds

help,noevents
help,imagerejected
help,error
help,events


stop

end