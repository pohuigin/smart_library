;Take the gradient in the horizontal plane of the LOS magnetic field.
;
;Example: datagrad=ar_losgrad(datasm)
;

function ar_losgrad, data

xgrad=deriv(data)

ygrad=rotate(deriv(rotate(data,-90)),90)

gradmag=sqrt(xgrad^2.+ygrad^2.)

return, gradmag

end