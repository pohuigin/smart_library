;----------------------------------------------------------------------------->

;Load the parameter meta data file for the detection algorithm.
function ar_loadparam, fparam=fparamin, meta=outmeta, silent=insilent
silent=keyword_set(insilent)

;Determine name of meta data file
if not keyword_set(fparamin) then fparam=ar_global(/fparam, silent=silent) else fparam=fparamin

if file_exist(fparam) ne 1 then begin
	if not silent then print,'% AR_LOADPARAM: The file FPARAM = '+fparam+' does not exist.'
	return,''
endif

if not silent then print,'% AR_LOADPARAM: AR Parameter File is '+fparam

;Read parameters from meta data file
readcol, fparam, param, val, type, meta, comment='#', format='A,A,A,A', delim=';',/silent
param=strtrim(param,2)
val=strtrim(val,2)
type=strtrim(type,2)
meta=strtrim(meta,2)

;Make array of data types for each field in structure
dataspec=strjoin(type,',')

;Create empty structure
create_struct, paramstruct, '', param, dataspec

;Fill the structure
for i=0,n_elements(param)-1 do paramstruct.(i)=val[i]

;Output the description of each parameter
outmeta=[[param],[meta]]

return, paramstruct

end

;----------------------------------------------------------------------------->