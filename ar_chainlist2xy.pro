;Convert ordered chain list of x1,y1,x2,y2,...,xn,yn into a x and y position array

function ar_chainlist2xy, chain

chainarr=str_sep(chain,',')
nchain=n_elements(chainarr)
chainxy=float(reform(chainarr,2,nchain/2.))

return,chainxy

end