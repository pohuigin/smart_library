;WRITE ROUTINE TO TAKE DET STRUCT READ FROM DET ARR AND BUFFER IT OUT SO THAT IT HAS SAME NUMBER OF 
;ELEMENTS AS POSPROPARR/MAGPROPARR

function ar_detstr2arstr,detstrarr,datafilearr

nar=n_elements(datafilearr)

ardetstrarr=detstrarr[0]

ardetstrarr=replicate(ardetstrarr,nar)


datas=datafilearr[sort(datafilearr)]
datau=datas[uniq(datas)]
nu=n_elements(datau)

for i=0,nu-1 do begin

	wmatch1=where(datafilearr eq datau[i])

	wmatch2=(where(detstrarr.datafile eq datau[i]))[0]

;DANGEROUS!!!!!!!! 
	if wmatch1[0] eq -1 or wmatch2[0] eq -1 then begin
		print,'!!!!!!! ar_detstr2arstr - datafile not matched in structure; i='+strtrim(i,2)+' !!!!!!!!!!'
		continue
	endif
	
	ardetstrarr[wmatch1]=detstrarr[wmatch2]

endfor

return,ardetstrarr

end