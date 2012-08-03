;A list of default arguments for various AR detection routines
;if no instrument is set, it assumes MDI
;needs the complete param structure as input (automatically reads it in if not supplied)
function ar_argdefault, inparams, instrument=ininst, $
	noisethresh=noisethresh, rotarg=rotarg

if n_elements(inparams) lt 1 then params = ar_loadparam() else params = inparams
	
if n_elements(ininst) lt 1 then inst='mdi' else inst=ininst

case inst of 
	'mdi' : 
	'hmi' :
	'gong' :
	'sot' :
endcase

;return the particular thresholds for a given instrument


return,retval