;----------------------------------------------------------------------------->

;Load the structure meta data file (e.g., ar_struct_param.txt).
;This routine is called in AR_STRUCT_INIT or the output is fed in.

function ar_loadstruct, fstruct=fstructin, silent=insilent
silent=keyword_set(insilent)

;Determine name of meta data file
if not keyword_set(fstructin) then fstruct=ar_global(/fstruct, silent=silent) else fstruct=fstructin

if file_exist(fstruct) ne 1 then begin
	if not silent then print,'% AR_LOADSTRUCT: The file FSTRUCT = '+fstruct+' does not exist.'
	return,''
endif

if not silent then print,'% AR_LOADSTRUCT: AR Parameter File is '+fstruct

;Read parameters from meta data file
readcol, fstruct, id, format, tags, meta, comment='#', format='A,A,A,A,A', delim=';',/silent
id=strtrim(id,2)
format=strtrim(format,2)
tags=strtrim(tags,2)
meta=strtrim(meta,2)

;Create empty structure
create_struct, paramstruct, '', ['id','format','tags','meta'], 'A,A,A,A'

paramstruct=replicate(paramstruct,n_elements(id))

;Fill the structure
paramstruct.id=id
paramstruct.format=format
paramstruct.tags=tags
paramstruct.meta=meta

return, paramstruct

end

;----------------------------------------------------------------------------->