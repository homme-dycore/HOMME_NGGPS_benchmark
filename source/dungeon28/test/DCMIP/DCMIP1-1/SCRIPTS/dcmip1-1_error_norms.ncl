load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   

begin

  ; parameters

  pi  = acos(-1.)
  g   = 9.80616d0       ; gravity (m/s^2)
  rad = pi/180.         ; radians per degree
  R   = 6.37122e6       ; radius of the earth

  ; open netcdf file

	f  = addfile("./movies/DCMIP1-11.nc","r")
 
  ; read grid information

	lat    = f->lat * rad
	lon    = f->lon * rad
	lev    = f->lev
  time   = f->time
  
	nlat   = getfilevardimsizes(f, "lat" )
	nlon   = getfilevardimsizes(f, "lon" )
	nlev   = getfilevardimsizes(f, "lev" )
  nt     = dimsizes(time)

  ; read initial and final states (time,lev,lat,long)

  q_i = f->Q(0,:,:,:)               ; initial state
	q_f = f->Q(nt-1,:,:,:)              ; final state
  dq  = q_f - q_i

  ; get vertical extent from geopotential
 
  h   = f->geo(0,:,0,0)/g              ; get level height
  dh  = new(nlev,typeof(h))
  base= 0.0*h(0)                        ; assume bottom is at z=0

  do i = 1,nlev
    dh(nlev-i) = 2.0*(h(nlev-i)-base)   ; get level width
    base       = base + dh(nlev-i)      ; move interface up
  end do

  ; get sample volumes

  dlat  = lat(1)-lat(0)
  dlon  = lon(1)-lon(0)
  dx = R*dlon
  dy = conform(q_i, R*cos(lat)*dlat, 1)
  dz = conform(q_i, dh, 0)
  dV = dx * dy * dz

  ; get error norms

  L1  =      sum(abs(dq)^1*dV) /      sum(abs(q_i - avg(q_i))^1*dV )
  L2  = sqrt(sum(abs(dq)^2*dV))/ sqrt(sum(abs(q_i - avg(q_i))^2*dV ))
  Linf=      max(abs(dq)  *dV) /      max(abs(q_i - avg(q_i))*dV   )
  q_max = max(q_f)
  q_min = min(q_f)
  
  print("DCMIP 1-1: L1="+L1+" L2= "+L2+" Linf= "+Linf+" q_max= "+q_max+" q_min= "+q_min)

end