function ar_blanknar, flare=flare, arstr=arstr, sfflare=sfflare, nlstr=nlstr, extentstr=extentstr, blank=blank, issi=issi, rotstr=rotstr, isingstr=isingstr, manolisstr=manolisstr

;NOAA
noaastr={time:0L,day:0,noaa:0,location:intarr(2),longitude:0,area:0,st$macintosh:byte(intarr(3)),long_ext:intarr(2),num_spots:intarr(2),st$mag_type:byte(intarr(16)),spare:byte(intarr(9))}
str=noaastr
;srs_str={name:'', loc:'', mtwil:'', mcint:'', area:0, lonlen:0, nspots:0, mu:0D, cor_area:0D}

;Polarity separation line
strnl={lnl:0d, lsg:0d, gradmax:0d, gradmean:0d, gradmedian:0d, $
	rval:0d, wlsg:0d, r_star:0d, wlsg_star:0d, $
	thetanl:0d, thetabcl:0d, lbcl:0d, npole:0d, polarity:0d}
if keyword_set(nlstr) then str=strnl

;AR rotation properties -> 
strrot={thetapsl:10000d, thetaspsl:10000d, thetabcl:10000d, lbcl:0d, npole:0d, polarity:0d, $ ;angle of weak neutral line, angle of strong neutral line, angle of bipole connecting line, length of bipole connecting line, number of polarity blobs, ???
	hglonlatpos:[10000d,10000d], hglonlatneg:[10000d,10000d], hglonlatlgpos:[10000d,10000d], hglonlatlgneg:[10000d,10000d], $ ;overall pos/neg centroids, largest pos/neg blob centroids
	xypos:[0d,0d],xyneg:[0d,0d],xylgpos:[0d,0d],xylgneg:[0d,0d], $ ;x/y pixel locations of overall and largest pos/neg blob centroids
	xy1psl:[0d,0d],xy2psl:[0d,0d],xy1spsl:[0d,0d],xy2spsl:[0d,0d]}
if keyword_set(rotstr) then str=strrot

;AR Ising energy
strising={energy:0d,energy_fits:0d,energyppx:0d,energy_fitsppx:0d}
if keyword_set(isingstr) then str=strising

;Defines Feature's Extent (NSEW & ~Sun center to edge)
extstr={xylon:[0D,0D],xylat:[0D,0D],rdeglon:[0D,0D],rdeglat:[0D,0D],hglon:[0D,0D],hglat:[0D,0D],xymean10lon:[0D,0D],xymean10lat:[0D,0D],hglonwidth:0D,hglatwidth:0D}
if keyword_set(extentstr) then return, extstr

;SEC Event Listing
if keyword_Set(flare) then str={eventnum:0L, start_time:'', max_time:'', end_time:'', satellite:'', q:'', type:'', freq:'', fclass:'', fbase:0., flux:0., p1:0., region:'', hglat:1000, hglon:1000}

;SMART AR Structure
if keyword_Set(arstr) then str={smid:'', id:'', class:'', type:['','',''], time:'', $
	hglon:10000D, hglat:10000D, hclon:10000D, hclat:10000D, carlon:10000D, carlat:10000D, $
	xpos:0D, ypos:0D, xbary:0D, ybary:0D, meanval:0D, stddv:0D, kurt:0D, narpx:0D, $
	bflux:0D, bfluxpos:0D, bfluxneg:0D, bfluxemrg:0D, area:0D, bmin:0D, bmax:0D, $
	nlstr:strnl, extstr:extstr} ;r_star:0D, wlsg_star:0D, schrijver_r:0D, wlsg:0D, nl_sg:0D, nl_length:0D, 

;Sam Freeland Flares
if keyword_set(sfflare) then str={date_obs:'',ename:'',class:'',fstart:'',fstop:'',fpeak:'',xcen:0,ycen:0,helio:'',lfiles:'',recok:0,url_index:'',url_movie:''}

;ISSI
if keyword_set(issi) then str={Event_Type: 'ActiveRegion', SMID: '', ID: '', CLASS: '', TYPE: '', Event_CoordSys: 'UTC-HPC-TOPO', $
Event_CoordUnit: 'arcsec, arcsec', Event_EndTime: '', Event_StartTime: '', Event_Coord1: 0D, Event_Coord2: 0D, Event_C1Error: 0D, Event_C2Error: 0D, $
HCLON: 0D, HCLAT: 0D, HGLON: 0D, HGLAT: 0D, CARLON: 0D, CARLAT: 0D, XPOS: 0D, YPOS: 0D, $
FRM_Contact: 'pohuigingmail', FRM_DateRun: '', FRM_HumanFlag: 'F', FRM_Identifier: 'phiggins', FRM_Institute: 'TCD', FRM_Name: 'SMART', FRM_ParamSet: '', FRM_URL: 'SM/smart_disk', $
OBS_Observatory: 'SOHO', OBS_ChannelID: 'V band', OBS_Instrument: 'MDI', OBS_MeanWavel: 6768, OBS_WaveUnit: 'angstroms', $
Bound_CCNsteps: 0D, Bound_CCStartC1: 0D, Bound_CCStartC2: 0D, Bound_ChainCode: '', ChainCodeType: 'ordered list', $
BoundBox_C1LL: 0D, BoundBox_C2LL: 0D, BoundBox_C1UR: 0D, BoundBox_C2UR: 0D, $ 
BoundHg_C1LL: 0D, BoundHg_C2LL: 0D, BoundHg_C1UR: 0D, BoundHg_C2UR: 0D, $ 
X1X2: '', Y1Y2: '', RDEGLON12: '', RDEGLAT12: '', HGLON12: '', HGLAT12: '', X1X2MEAN10: '', Y1Y2MEAN10: '', HGLONW: 0D, HGLATW: 0D, $
MEANVAL: 0D, STDDV: 0D, KURT: 0D, Event_Npixels: 0D, BFLUX: 0D, BFLUXPOS: 0D, BFLUXNEG: 0D, BFLUXEMRG: 0D, Area_AtDiskCenter: 0D, Area_Unit: 'Mm^2', BMIN: 0D, BMAX: 0D, $
LNL: 0D, LSG: 0D, GRADMAX: 0D, GRADMEAN: 0D, GRADMEDIAN: 0D, RVAL: 0D, WLSG: 0D, R_STAR: 0D, WLSG_STAR: 0D, $
thetanl:0d, thetasnl:0d, thetabcl:0d, lbcl:0d, npole:0d, $
isinge:0d, isinge_ppx:0d}


return, str

end