;Output a blank structure and format key (string) corresponding to the input STRUCTID.

function ar_struct_init,format,structid=structid, silent=silent, fstruct=fstruct

paramstruct=ar_loadstruct(fstruct=fstruct, silent=silent)

wid=where(strlowcase(strtrim(paramstruct.id,2)) eq strlowcase(strtrim(structid,2)))

if wid[0] eq -1 then message,'% AR_STRUCT_INIT: Structure ID not found in '+ar_global(/fstruct)

thisparam=paramstruct[wid]

format=thisparam.format
formatarr=str_sep(format,',')
tags=str_sep(thisparam.tags,',')

;if n_elements(tags) gt 20 then begin
;	create_struct,strblank,'',tags[0:19],format2create_struct(strjoin(formatarr[0:19],','))
;	create_struct,strblank2,'',tags[20:*],format2create_struct(strjoin(formatarr[20:*],','))
;	strblank=merge_struct(strblank,strblank2)
;endif else begin
	create_struct,strblank,'',tags,format2create_struct(format)
;endelse

outstrblank=strblank

return,outstrblank

end