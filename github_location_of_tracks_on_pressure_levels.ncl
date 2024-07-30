
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;********************************************
begin

f  = addfile("/public/home/wangji19/data/track_data/TRACK_ff_trs_pos_ERA5_GRID1_1950101_20230720_NH_CYC.nc","r")
; f  = addfile("/cygdrive/e/data/track_data/TRACK_ff_trs_pos_ERA5_GRID1_1950101_20230720_NH_CYC.nc","r")

tlon=f->traj_lon
tlat=f->traj_lat
vor=f->traj_rel_vor_pos
time=f->traj_time 
time1=time*6
time1@units="hours since 1950-01-01 00:00:0.0"
time1@long_name="Time of the trajectory"

time1@_FillValue=-2147483647
tlon@_FillValue=-9.96921e+36
tlat@_FillValue=-9.96921e+36
; TT1=tt%10000
; t=ispan(1950,2022,1)
delete(time)
delete(f)


 data=asciiread("/public/home/wangji19/data/ERA5_TRACK/grid1/1950-2022/MAMmg_1950-2022.txt",-1,"integer")
 npts=dimsizes(data)
print("Number of values: "+npts)

tttime1=cd_calendar(time1(:,0),-5)
; x=ind(tttime1(data,0).eq.2021)
x=ind(tttime1(data,0).ge.2011.and.tttime1(data,0).le.2022)
x=data(x);
  delete(tttime1)
  delete(data) 
  
asciiwrite("/public/home/wangji19/program3/vor/mg_id_2011_2022.txt",x)

; exit
dim=dimsizes(x)
 tlon_p=tlon(x,0:149)
tlat_p=tlat(x,0:149)
  Tt=cd_calendar(time1(x,0:149),-3)
  delete(tlon) 
  delete(tlat)
  delete(time1)
  
;;;-----------------------------------------------
 ; id =  115304;;120886
 ; xx=ind(id.eq.x)
 
 ; print(x(xx))
 ; Tt=cd_calendar(time1(i,:),-3)
; lat1=lat(i,:)
; lon1=lon(i,:)
; vor1=vor(i,:)
; numm=num(.not.ismissing(lat1))
 ;;;;;;;;;;
 	
fvor  = addfile("/public/home/wangji19/data/era5_daily/vor/vor.2011-2022.nc","r")
; fvor  = addfile("/public/home/wangji19/data/era5_daily/vor/2021.vor0.25.nc","r")

 tim   = fvor->time                          ; time:units = "hours since 1-1-1 00:00:0.0"                               
  Time2=cd_calendar(tim,-3)
delete(tim)
  p= short2flt(fvor->vo(:,:,:,:))
   Pp= p(:,:,::-1,:)
  delete(fvor)
  delete(p)
printVarSummary(Pp)
  
   dims = dimsizes(Pp)
      lat1d = Pp&latitude 
      lon1d = Pp&longitude 
	  level = Pp&level 
      misva = Pp@_FillValue

  path_lat = new((/dim,150,7/),float)
  path_lon = new((/dim,150,7/),float)
;;;------------900hPa-----------------
  do o=0, dim-1;;; 气旋id 
  ; print(o)
  ; ; nnn=num(.not.ismissing(tlat_p(o,:)))
  ; o=xx
  y1=ind((tlon_p(o,:).ge.75).and.(tlon_p(o,:).le.145))

  ; exit
  y2=ind((tlat_p(o,:).ge.25).and.(tlat_p(o,:).le.65))
  yy=venn2_intersection(y1,y2)
 
 ; print(tlon_p(o,0:50))
 ; print(y1)
 ; print(y2)
 ; print(yy)
 ; exit
  do i=0,dimsizes(yy)-1;;;;;
    idx1	 =	ind(Tt(o,yy(i)).eq.Time2(:));;;;
    londex   =  ind( toint(tlon_p(o,yy(i)) ).eq.lon1d(:)   )
; print(yy)

    fvor = new((/dims(2),dims(3)/),float)
	fvor@_FillValue = 32767
    fvor(:,:)  =   Pp(idx1,8,:,:)
   copy_VarMeta(Pp(0,0,:,:),fvor(:,:))
    fvor(:,:londex)=fvor@_FillValue;;;;;
 ; printVarSummary(fvor)
    ; do j=0,dims(2)-1
     clat = new(25,float)    ; Create arrays to hold great circle.
     clon = new(25,float)    ; As of NCL 6.5.0, can be float or double

    nggcog(tlat_p(o,yy(i)),tlon_p(o,yy(i)),5.0,clat,clon);;;
	min_lat = min(clat)
    min_lon = min(clon)
    max_lat = max(clat)	
    max_lon = max(clon)
 
    newTRMMyearly := fvor({min_lat:max_lat},{min_lon:max_lon})
    delete(fvor)
 ; printVarSummary(newTRMMyearly)
    lat2d  = conform(newTRMMyearly,newTRMMyearly&latitude,0)
    lon2d  = conform(newTRMMyearly,newTRMMyearly&longitude,1)
	in_circle  = gc_inout(lat2d,lon2d,clat,clon)
    newTRMMyearly = where(in_circle,newTRMMyearly,newTRMMyearly@_FillValue)
    x1D = ndtooned(newTRMMyearly)     ; only do this once
	xMax = max(newTRMMyearly);;;;;;;;;
	xind = ind(x1D.eq.xMax)
    inds = ind_resolve(xind(0),dimsizes(newTRMMyearly)) 
    ; inds = ind_resolve(maxind (x1D), dimsizes(newTRMMyearly))    ; convert 1D array back to 2D 
        
    ; print(tlon_p(o,i))
	; exit
	ilat = inds(0,0)        ; select the latitude index where the X array is at its' maximum  
    ilon = inds(0,1)        ; select the longitude index where the X array is at its' maximum

    lat_max = newTRMMyearly&latitude(ilat)   ; insert the latitude index into the lat coordinate variable
    lon_max = newTRMMyearly&longitude(ilon)   ; insert the longitude index into the lon coordinate variable
    
	path_lat(o,yy(i),5) = lat_max;;
	path_lon(o,yy(i),5) = lon_max;;
	 
   delete(newTRMMyearly)
   delete(idx1)
   delete(londex)
   delete(lat2d)
   delete(lon2d)
   delete(clat)
   delete(clon)
   delete(x1D)
   delete(in_circle)
   delete(xind)
  end do
    print(yy )
     delete(y1)
   delete(y2)
   delete(yy)
   
 end do

 ;;;------------1000hPa-----------------
 do o=0 ,dim-1;;
   y1=ind((tlon_p(o,:).ge.75).and.(tlon_p(o,:).le.145))
  y2=ind((tlat_p(o,:).ge.25).and.(tlat_p(o,:).le.65))
  yy=venn2_intersection(y1,y2)
 
  do i=0,dimsizes(yy)-1;;;;;
    idx1	 =	ind(Tt(o,yy(i)).eq.Time2(:));;;;	
    londex   =  ind( toint(path_lon(o,yy(i),5)).eq.lon1d(:)   )

    fvor = new((/dims(2),dims(3)/),float)
	fvor@_FillValue = 32767

    fvor(:,:)  =   Pp(idx1,9,:,:);;;
    copy_VarMeta(Pp(0,0,:,:),fvor(:,:))
    fvor(:,:londex)=fvor@_FillValue;;;;;
 printVarSummary(fvor)
    ; do j=0,dims(2)-1
     clat = new(25,float)    ; Create arrays to hold great circle.
     clon = new(25,float)    ; As of NCL 6.5.0, can be float or double

    nggcog(path_lat(o,yy(i),5) ,path_lon(o,yy(i),5) ,5.0,clat,clon);;;
	min_lat = min(clat)
    min_lon = min(clon)
    max_lat = max(clat)	
    max_lon = max(clon)
 
    newTRMMyearly := fvor({min_lat:max_lat},{min_lon:max_lon})
    delete(fvor)
 printVarSummary(newTRMMyearly)
    lat2d  = conform(newTRMMyearly,newTRMMyearly&latitude,0)
    lon2d  = conform(newTRMMyearly,newTRMMyearly&longitude,1)
	in_circle  = gc_inout(lat2d,lon2d,clat,clon)
    newTRMMyearly = where(in_circle,newTRMMyearly,newTRMMyearly@_FillValue)
    x1D = ndtooned(newTRMMyearly)     ; only do this once
	xMax = max(newTRMMyearly)
    xind = ind(x1D.eq.xMax)
    inds = ind_resolve(xind(0),dimsizes(newTRMMyearly)) 
    ilat = inds(0,0)        ; select the latitude index where the X array is at its' maximum  
    ilon = inds(0,1)        ; select the longitude index where the X array is at its' maximum
   
   lat_max = newTRMMyearly&latitude(ilat)   ; insert the latitude index into the lat coordinate variable
    lon_max = newTRMMyearly&longitude(ilon)   ; insert the longitude index into the lon coordinate variable
    
	path_lat(o,yy(i),6) = lat_max;;
	path_lon(o,yy(i),6) = lon_max;;
  
   delete(newTRMMyearly)
   delete(idx1)
   delete(londex)
   delete(lat2d)
   delete(lon2d)
   delete(clat)
   delete(clon)
   delete(x1D)
   delete(in_circle)
   delete(xind)
   
  end do	
     delete(y1)
   delete(y2)
   delete(yy)
   
 end do
 
 ;;;------------850hPa-----------------
 do o=0 ,dim-1;;
   y1=ind((tlon_p(o,:).ge.75).and.(tlon_p(o,:).le.145))
  y2=ind((tlat_p(o,:).ge.25).and.(tlat_p(o,:).le.65))
  yy=venn2_intersection(y1,y2)
 
  do i=0,dimsizes(yy)-1;;;;;
 
 	path_lat(o,yy(i),4) = tlat_p(o,yy(i));;
	path_lon(o,yy(i),4) = tlon_p(o,yy(i));;
	
   end do	
   delete(y1)
   delete(y2)
   delete(yy)
 end do
 ;;;-------------800hPa--------------------------   

 do o=0  ,dim-1;;
   y1=ind((tlon_p(o,:).ge.75).and.(tlon_p(o,:).le.145))
  y2=ind((tlat_p(o,:).ge.25).and.(tlat_p(o,:).le.65))
  yy=venn2_intersection(y1,y2)
 
  do i=0,dimsizes(yy)-1;;;;;
    idx1	 =	ind(Tt(o,yy(i)).eq.Time2(:));;;;	
    londex   =  ind( toint(tlon_p(o,yy(i)) ).eq.lon1d(:)   )

    fvor = new((/dims(2),dims(3)/),float)
	fvor@_FillValue = 32767

	fvor(:,:)  =   Pp(idx1,7,:,:);;;
    copy_VarMeta(Pp(0,0,:,:),fvor(:,:))
    fvor(:,londex:)=fvor@_FillValue;;;;;
 printVarSummary(fvor)
    ; do j=0,dims(2)-1
     clat = new(25,float)    ; Create arrays to hold great circle.
     clon = new(25,float)    ; As of NCL 6.5.0, can be float or double
  
    nggcog(tlat_p(o,yy(i)),tlon_p(o,yy(i)),5.0,clat,clon);;;
	min_lat = min(clat)
    min_lon = min(clon)
    max_lat = max(clat)	
    max_lon = max(clon)
	    
    newTRMMyearly := fvor({min_lat:max_lat},{min_lon:max_lon})
    delete(fvor)
 printVarSummary(newTRMMyearly)
    lat2d  = conform(newTRMMyearly,newTRMMyearly&latitude,0)
    lon2d  = conform(newTRMMyearly,newTRMMyearly&longitude,1)
	in_circle  = gc_inout(lat2d,lon2d,clat,clon)
    newTRMMyearly = where(in_circle,newTRMMyearly,newTRMMyearly@_FillValue)
    x1D = ndtooned(newTRMMyearly)     ; only do this once
	xMax = max(newTRMMyearly)
    xind = ind(x1D.eq.xMax)
    inds = ind_resolve(xind(0),dimsizes(newTRMMyearly)) 
    ilat = inds(0,0)        ; select the latitude index where the X array is at its' maximum  
    ilon = inds(0,1)        ; select the longitude index where the X array is at its' maximum
    lat_max = newTRMMyearly&latitude(ilat)   ; insert the latitude index into the lat coordinate variable
    lon_max = newTRMMyearly&longitude(ilon)   ; insert the longitude index into the lon coordinate variable

	path_lat(o,yy(i),3) = lat_max;;
	path_lon(o,yy(i),3) = lon_max;;
	
   delete(newTRMMyearly)
   delete(idx1)
  
  delete(lat2d)
   delete(lon2d)
  delete(clat)
   delete(clon)
   delete(x1D)
   delete(in_circle)
   delete(xind)
  end do	
     delete(y1)
   delete(y2)
   delete(yy)
 end do
 
 ;;;-------------700hPa--------------------------   

 do o=0 ,dim-1;;
   y1=ind((tlon_p(o,:).ge.75).and.(tlon_p(o,:).le.145))
  y2=ind((tlat_p(o,:).ge.25).and.(tlat_p(o,:).le.65))
  yy=venn2_intersection(y1,y2)
 
  do i=0,dimsizes(yy)-1;;;;;
    idx1	 =	ind(Tt(o,yy(i)).eq.Time2(:));;;;	
    londex   =  ind( toint(path_lon(o,yy(i),3) ).eq.lon1d(:)   )
   
    fvor = new((/dims(2),dims(3)/),float)
	fvor@_FillValue = 32767
	fvor(:,:)  =   Pp(idx1,6,:,:);;;
    copy_VarMeta(Pp(0,0,:,:),fvor(:,:))
    fvor(:,londex:)=fvor@_FillValue;;;;;
 printVarSummary(fvor)
    ; do j=0,dims(2)-1
     clat = new(25,float)    ; Create arrays to hold great circle.
     clon = new(25,float)    ; As of NCL 6.5.0, can be float or double
  
    nggcog(path_lat(o,yy(i),3),path_lon(o,yy(i),3),5.0,clat,clon);;;半径
	min_lat = min(clat)
    min_lon = min(clon)
    max_lat = max(clat)	
    max_lon = max(clon)
	    
    newTRMMyearly := fvor({min_lat:max_lat},{min_lon:max_lon})
    delete(fvor)
 printVarSummary(newTRMMyearly)
    lat2d  = conform(newTRMMyearly,newTRMMyearly&latitude,0)
    lon2d  = conform(newTRMMyearly,newTRMMyearly&longitude,1)
	in_circle  = gc_inout(lat2d,lon2d,clat,clon)
    newTRMMyearly = where(in_circle,newTRMMyearly,newTRMMyearly@_FillValue)
    x1D = ndtooned(newTRMMyearly)     ; only do this once
	xMax = max(newTRMMyearly)
   xind = ind(x1D.eq.xMax)
    inds = ind_resolve(xind(0),dimsizes(newTRMMyearly)) 
    ilat = inds(0,0)        ; select the latitude index where the X array is at its' maximum  
    ilon = inds(0,1)        ; select the longitude index where the X array is at its' maximum
    lat_max = newTRMMyearly&latitude(ilat)   ; insert the latitude index into the lat coordinate variable
    lon_max = newTRMMyearly&longitude(ilon)   ; insert the longitude index into the lon coordinate variable

	path_lat(o,yy(i),2) = lat_max;;
	path_lon(o,yy(i),2) = lon_max;;
	
	
   delete(newTRMMyearly)
   delete(idx1)
  
   delete(londex)
  delete(lat2d)
   delete(lon2d)
  delete(clat)
   delete(clon)
   delete(x1D)
   delete(in_circle)
   delete(xind)
  end do		
    delete(y1)
   delete(y2)
   delete(yy) 
 end do
  
 ;;;-------------600hPa--------------------------   

  do o=0 ,dim-1;;
    y1=ind((tlon_p(o,:).ge.75).and.(tlon_p(o,:).le.145))
  y2=ind((tlat_p(o,:).ge.25).and.(tlat_p(o,:).le.65))
  yy=venn2_intersection(y1,y2)
 
  do i=0,dimsizes(yy)-1;;;;;
    idx1	 =	ind(Tt(o,yy(i)).eq.Time2(:));;;;	
		
    londex   =  ind( toint(path_lon(o,yy(i),2) ).eq.lon1d(:)   )

    fvor = new((/dims(2),dims(3)/),float)
	fvor@_FillValue = 32767
	fvor(:,:)  =   Pp(idx1,5,:,:);;;
    copy_VarMeta(Pp(0,0,:,:),fvor(:,:))
	 print(londex)
    fvor(:,londex:)=fvor@_FillValue;;;;;
 printVarSummary(fvor)
    ; do j=0,dims(2)-1
     clat = new(25,float)    ; Create arrays to hold great circle.
     clon = new(25,float)    ; As of NCL 6.5.0, can be float or double
  
    nggcog(path_lat(o,yy(i),2),path_lon(o,yy(i),2),5.0,clat,clon);;;半径
	min_lat = min(clat)
    min_lon = min(clon)
    max_lat = max(clat)	
    max_lon = max(clon)
	    
    newTRMMyearly := fvor({min_lat:max_lat},{min_lon:max_lon})
    delete(fvor)
 printVarSummary(newTRMMyearly)
    lat2d  = conform(newTRMMyearly,newTRMMyearly&latitude,0)
    lon2d  = conform(newTRMMyearly,newTRMMyearly&longitude,1)
	in_circle  = gc_inout(lat2d,lon2d,clat,clon)
    newTRMMyearly = where(in_circle,newTRMMyearly,newTRMMyearly@_FillValue)
    x1D = ndtooned(newTRMMyearly)     ; only do this once
	xMax = max(newTRMMyearly)
    xind = ind(x1D.eq.xMax)
    inds = ind_resolve(xind(0),dimsizes(newTRMMyearly)) 
    ilat = inds(0,0)        ; select the latitude index where the X array is at its' maximum  
    ilon = inds(0,1)        ; select the longitude index where the X array is at its' maximum
    lat_max = newTRMMyearly&latitude(ilat)   ; insert the latitude index into the lat coordinate variable
    lon_max = newTRMMyearly&longitude(ilon)   ; insert the longitude index into the lon coordinate variable
    
	path_lat(o,yy(i),1) = lat_max;;
	path_lon(o,yy(i),1) = lon_max;;
	
   delete(newTRMMyearly)
   delete(idx1)
  
  delete(lat2d)
   delete(lon2d)
  delete(clat)
   delete(clon)
   delete(x1D)
   delete(in_circle)
   delete(xind)
  end do	
     delete(y1)
   delete(y2)
   delete(yy)
 end do
  
 ;;;-------------500hPa--------------------------   

 do o=0 ,dim-1;;
    y1=ind((tlon_p(o,:).ge.75).and.(tlon_p(o,:).le.145))
  y2=ind((tlat_p(o,:).ge.25).and.(tlat_p(o,:).le.65))
  yy=venn2_intersection(y1,y2)
 
  do i=0,dimsizes(yy)-1;;;;;
    idx1	 =	ind(Tt(o,yy(i)).eq.Time2(:));;;;	
    londex   =  ind( toint(path_lon(o,yy(i),1) ).eq.lon1d(:)   )

    fvor = new((/dims(2),dims(3)/),float)
	fvor@_FillValue = 32767
	fvor(:,:)  =   Pp(idx1,4,:,:);;;
    copy_VarMeta(Pp(0,0,:,:),fvor(:,:))
    fvor(:,londex:)=fvor@_FillValue;;;;;
 printVarSummary(fvor)
    ; do j=0,dims(2)-1
     clat = new(25,float)    ; Create arrays to hold great circle.
     clon = new(25,float)    ; As of NCL 6.5.0, can be float or double
  
    nggcog(path_lat(o,yy(i),1),path_lon(o,yy(i),1),5.0,clat,clon);;;
	min_lat = min(clat)
    min_lon = min(clon)
    max_lat = max(clat)	
    max_lon = max(clon)
	    
    newTRMMyearly := fvor({min_lat:max_lat},{min_lon:max_lon})
    delete(fvor)
 printVarSummary(newTRMMyearly)
    lat2d  = conform(newTRMMyearly,newTRMMyearly&latitude,0)
    lon2d  = conform(newTRMMyearly,newTRMMyearly&longitude,1)
	in_circle  = gc_inout(lat2d,lon2d,clat,clon)
    newTRMMyearly = where(in_circle,newTRMMyearly,newTRMMyearly@_FillValue)
    x1D = ndtooned(newTRMMyearly)     ; only do this once
	xMax = max(newTRMMyearly)
    xind = ind(x1D.eq.xMax)
    inds = ind_resolve(xind(0),dimsizes(newTRMMyearly)) 
    ilat = inds(0,0)        ; select the latitude index where the X array is at its' maximum  
    ilon = inds(0,1)        ; select the longitude index where the X array is at its' maximum
    lat_max = newTRMMyearly&latitude(ilat)   ; insert the latitude index into the lat coordinate variable
    lon_max = newTRMMyearly&longitude(ilon)   ; insert the longitude index into the lon coordinate variable
    
	path_lat(o,yy(i),0) = lat_max;;
	path_lon(o,yy(i),0) = lon_max;;
	
   delete(newTRMMyearly)
   delete(idx1)
  
   delete(londex)
  delete(lat2d)
   delete(lon2d)
  delete(clat)
   delete(clon)
   delete(x1D)
   delete(in_circle)
   delete(xind)
  end do	
     delete(y1)
   delete(y2)
   delete(yy)
 end do
  
 printVarSummary(path_lat)

	outfil1=addfile("/public/home/wangji19/program3/vor/mgc_pressure_vor_latind5_Strict_Westward_2011-2022.nc","c")  
    outfil1->path_lat=path_lat
    outfil2=addfile("/public/home/wangji19/program3/vor/mgc_pressure_vor_lonind5_Strict_Westward_2011-2022.nc","c")  
    outfil2->path_lon=path_lon

 end