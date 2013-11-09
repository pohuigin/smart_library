;Output a blank structure and format key (string) corresponding to the input STRUCTID.

function ar_struct_init,format,structid=structid, silent=silent, fstruct=fstruct

paramstruct=ar_loadstruct(fstruct=fstruct, silent=silent)

wid=where(strlowcase(strtrim(paramstruct.id,2)) eq strlowcase(strtrim(structid,2)))

if wid[0] eq -1 then message,'% AR_STRUCT_INIT: Structure ID not found in '+ar_global(/fstruct)

thisparam=paramstruct[wid]

format=thisparam.format
tags=str_sep(thisparam.tags,',')

create_struct,strblank,'',tags,format2create_struct(format)

outstrblank=strblank

return,outstrblank

end