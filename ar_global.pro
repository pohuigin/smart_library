;----------------------------------------------------------------------------->

function ar_global, fparam=fparam, proot=proot, fstruct=fstruct, silent=insilent
silent=keyword_set(insilent)

;test if setup has been run. if not, then run it.
ar_setup,status=status,/get

if not status then begin 
	if not silent then print,'% AR_GLOBAL: AR_SETUP has not been run. Running with default paths and config...'
	ar_setup 
	ar_setup,defvar=defvar,/get
	if not silent then print,'% AR_GLOBAL: '+strjoin(defvar.(0),' = ')
	if not silent then print,'% AR_GLOBAL: '+strjoin(defvar.(1),' = ')
	if not silent then print,'% AR_GLOBAL: '+strjoin(defvar.(2),' = ')
endif

retval=''

if keyword_set(fparam) then retval=!AR_PATH+!AR_PARAM ;'~/science/procedures/chole_detect/chole_param.txt'

;if n_elements(pdata) gt 0 then retval='~/science/procedures/chole_detect/data/'

if keyword_set(proot) then retval=!AR_PATH

if keyword_set(fstruct) then retval=!AR_PATH+!AR_STRUCT


return, retval

end

;----------------------------------------------------------------------------->