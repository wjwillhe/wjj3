load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************

undef("t_test2")
function t_test2(x1[*][*][*]:float,x2[*][*][*]:float)
local ave_x1,ave_x2,n1,n2,s1_ts,s2_ts,s1s2_ts,t_ts,prob_phase,df,confi
begin
	ave_x1 = dim_avg_n_Wrap(x1(:,:,:),0);样本均值
	ave_x2 = dim_avg_n_Wrap(x2(:,:,:),0);样本均值
	    n1 = dimsizes(x1(:,0,0))
	    n2 = dimsizes(x2(:,0,0))
	df=n1+n2-2
	
	s1_ts       =dim_variance_n(x1(:,:,:),0) ;样本方差
    s2_ts       =dim_variance_n(x2(:,:,:),0)  ;样本方差
    s1s2_ts    =(n1-1)*s1_ts+(n2-1)*s2_ts   

	prob_phase = new((/dimsizes(x1&latitude),dimsizes(x1&longitude)/),float)
	prob_phase!0 = "latitude"
	prob_phase!1 = "longitude"
	prob_phase&latitude = x1&latitude
	prob_phase&longitude = x1&longitude
	
	do j = 0,dimsizes(x1&latitude)-1
		do jj = 0,dimsizes(x1&longitude)-1		     
				if(ismissing(s1s2_ts(j,jj)))then						
				        s1s2_ts(j,jj)=0
				end if
		end do
	end do
	
	do j = 0,dimsizes(x1&latitude)-1
		do jj = 0,dimsizes(x1&longitude)-1		     
				if(s1s2_ts(j,jj).ne.0)then
						prob_phase(j,jj) =  (ave_x1(j,jj)-ave_x2(j,jj))/(sqrt(s1s2_ts(j,jj)/df)*sqrt(1.0/n1+1.0/n2))
				
				end if
		end do
	end do
    confi=ave_x1
	confi=student_t(prob_phase,df)
	;print(ave_sample)
	; print(confi)
	;ave_sample(:,:) = mask(ave_sample(:,:),confi(:,:).gt.sig,False)
	; ave_sample(:,:) = where(confi(:,:).gt.sig,0,ave_sample(:,:))
    ;print(ave_sample)
	return(confi)
	
end


begin
; ;;;
fw=addfile("/public/home/wangji19/data/era5_daily/300w.daily.1950-2023.nc","r")
  w= dble2flt(fw->w(:,:,{0:180}))

ft=addfile("/public/home/wangji19/data/era5_daily/300t.daily.1950-2023.nc","r")
  air=ft->t(:,:,{0:180})

fwl=addfile("/public/home/wangji19/data/era5_daily/300wClmDay.daily.1950-2023.nc","r")
  wl=fwl->wClmDay(:,:,{0:180})

ftl=addfile("/public/home/wangji19/data/era5_daily/300tClmDay.daily.1950-2023.nc","r")
  airl=ftl->tClmDay(:,:,{0:180})

delete(fw)
delete(ft)
delete(fwl)
delete(ftl)
printVarSummary(w)
printVarSummary(air)

dims = dimsizes(w)
level = 850
latt = w&latitude 
lont = w&longitude
  time=w&time
Time3  = cd_calendar(time, -2)
Time3@_FillValue=-32767
delete(time)
 ;;;; daily annual cycle
  t_calend = runave_n_Wrap ( airl, 31,0, 0)
  w_calend = runave_n_Wrap ( wl, 31,0, 0)
	mmdd    = Time3%10000
     md     = airl&time 
	 d1      = dimsizes(md)
delete(airl)
delete(wl)
printVarSummary(t_calend) 

; print(u_calend(20:25,2:3,0))
;calculate uv anomaly
w_ano = new((/dims(0),dims(1),dims(2)/),float)
w_ano@_FillValue=-32767
t_ano = new((/dims(0),dims(1),dims(2)/),float)
t_ano@_FillValue=-32767

 do  i = 0,dims(0)-1
      ide1                 = ind(mmdd(i).eq.md(:))
      w_ano(i,:,:)           = w(i,:,:)-w_calend(ide1,:,:)
      t_ano(i,:,:)           = air(i,:,:)-t_calend(ide1,:,:)

	   delete(ide1)
  end do 

copy_VarMeta(w(:,:,:),w_ano(:,:,:))
copy_VarMeta(air(:,:,:),t_ano(:,:,:))
; print(u_ano(20:25,2:3,0))
delete(w)
delete(air)

;filter high-frequency 
nwt = 13
ihp = 2;band pass
fca = 1./8.
fcb = 1./2.5
nsigma = 1
wgt = filwgts_lanczos(nwt, ihp, fca, fcb, nsigma)
w_h_0 = wgt_runave_n(w_ano, wgt, 0, 0)
t_h_0 = wgt_runave_n(t_ano, wgt, 0, 0)
copy_VarMeta(w_ano(:,:,:),w_h_0(:,:,:))
copy_VarMeta(t_ano(:,:,:),t_h_0(:,:,:))
delete(w_ano)
delete(t_ano)
delete(w_calend)
printVarSummary(w_h_0) 

;calculate kp 2 ke
bc = w_h_0 * t_h_0
;summarize
para = -1.*((1000./level)^0.72)*287/9.8
p2k_cli = conform(bc,para,1) * bc

copy_VarMeta(w_h_0(:,:,:),p2k_cli(:,:,:))
printVarSummary(p2k_cli) 

 setfileoption("nc","Format","LargeFile") 
outfil=addfile("/public/home/wangji19/program3/pv/300p2k_cli.daily.1950-2023.nc","c")  
  outfil->p2k_cli=p2k_cli

  
dims=dimsizes(p2k_cli)
  time = p2k_cli&time
Time3  = cd_calendar(time, -2)
Time3@_FillValue=-32767

cy = addfile("/public/home/wangji19/data/track_data/TRACK_ff_trs_pos_ERA5_GRID1_1950101_20230720_NH_CYC.nc","r")
;;;;;;;;;;;气旋数据,各类路径气旋的位置
lon=cy->traj_lon
lat=cy->traj_lat
vor=cy->traj_rel_vor_pos

lon@_FillValue=-999.0
lat@_FillValue=-999.0 
 vor@_FillValue=-999.0
 
time1=cy->traj_time 
time1_FillValue=-2147483647
time1=time1*6
time1@units="hours since 1950-01-01 00:00:0.0"
time1@long_name="Time of the trajectory"
Tim1=cd_calendar(time1,-2)
delete(time1)


 ;;;;;;;;;;;;ALL
 data=asciiread("/public/home/wangji19/program3/vor/mg_id.txt",-1,"integer")
xx=data;
nnumba=dimsizes(xx)
patha=ispan(0,nnumba-1,1)
 
    auw=new((/nnumba,dims(1),dims(2)/),float)
  do i = 0,nnumba-1
   y1   =  ind((lon(xx(i),:).ge.75).and.(lon(xx(i),:).le.145))
   y2   =  ind((lat(xx(i),:).ge.25).and.(lat(xx(i),:).le.65))
   yy   =  venn2_intersection(y1,y2)
   mmax =  maxind(vor(xx(i),yy))
   mmax =  yy(mmax)
   
    nmind =  ind( Tim1(xx(i),mmax).eq.Time3(:))
  auw(i,:,:) = p2k_cli(nmind,:,:)
  
   delete(y1)
   delete(y2)
   delete(yy)
   delete(mmax)
   delete(nmind)
  end do	 
  
  auw!0="path"
  auw&path=patha
copy_VarMeta(p2k_cli(0,:,:),auw(0,:,:))

  u_compa=new((/dims(1),dims(2)/),float)
  u_compa(:,:) =dim_avg_n_Wrap(auw(:,:,:),0)
copy_VarMeta(p2k_cli(0,:,:),u_compa(:,:))
print(max(u_compa))
print(min(u_compa))
 
 ;;;;;;;;;;;;strong
 dat1=asciiread("/public/home/wangji19/program3/vor/mg_id_maxvor_33.txt",-1,"integer")
 x1=dat1;
nnumb1=dimsizes(x1)
path1=ispan(0,nnumb1-1,1)
 
    uw12=new((/nnumb1,dims(1),dims(2)/),float)
  do i = 0,nnumb1-1
   y1   =  ind((lon(x1(i),:).ge.75).and.(lon(x1(i),:).le.145))
   y2   =  ind((lat(x1(i),:).ge.25).and.(lat(x1(i),:).le.65))
   yy   =  venn2_intersection(y1,y2)
   mmax =  maxind(vor(x1(i),yy))
   mmax =  yy(mmax)
   
    nmind =  ind( Tim1(x1(i),mmax).eq.Time3(:))
  uw12(i,:,:) = p2k_cli(nmind,:,:)
  
   delete(y1)
   delete(y2)
   delete(yy)
   delete(mmax)
   delete(nmind)
  end do	 
  
  uw12!0="path"
  uw12&path=path1
copy_VarMeta(p2k_cli(0,:,:),uw12(0,:,:))

  u_comp1=new((/dims(1),dims(2)/),float)
  u_comp1(:,:) =dim_avg_n_Wrap(uw12(:,:,:),0)
copy_VarMeta(p2k_cli(0,:,:),u_comp1(:,:))
print(max(u_comp1))
print(min(u_comp1))
;;;;;;;;;;;;;;;;;;;;;

confi_com1 = new((/dims(1),dims(2)/),float)
confi_com1 = t_test2(uw12,auw)
copy_VarMeta(u_comp1(:,:),confi_com1(:,:))
print(max(confi_com1))
print(min(confi_com1))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;----------------------------

 ;;;;;;;;;;;;weak
dat2=asciiread("/public/home/wangji19/program3/vor/mg_id_maxvor_-33.txt",-1,"integer")
x2=dat2;
nnumb2=dimsizes(x2)
path2=ispan(0,nnumb2-1,1)
 
    uw22=new((/nnumb2,dims(1),dims(2)/),float)
  do i = 0,nnumb2-1
   y1   =  ind((lon(x2(i),:).ge.75).and.(lon(x2(i),:).le.145))
   y2   =  ind((lat(x2(i),:).ge.25).and.(lat(x2(i),:).le.65))
   yy   =  venn2_intersection(y1,y2)
   mmax =  maxind(vor(x2(i),yy))
   mmax =  yy(mmax)
  nmind =  ind( Tim1(x2(i),mmax).eq.Time3(:))
  uw22(i,:,:) = p2k_cli(nmind,:,:)
   delete(y1)
   delete(y2)
   delete(yy)
   delete(mmax)
   delete(nmind)
  end do	 
 
  uw22!0="path"
  uw22&path=path2
copy_VarMeta(p2k_cli(0,:,:),uw22(0,:,:))

  u_comp2=new((/dims(1),dims(2)/),float)
  u_comp2(:,:) =dim_avg_n_Wrap(uw22(:,:,:),0)
copy_VarMeta(p2k_cli(0,:,:),u_comp2(:,:))

print(max(u_comp2))
print(min(u_comp2))
; ;;;;;;;;;;;;;;;;;;------------------------------
;;;;;;;;;;;;Composite field significance test;;;;;;;;;;;

confi_com2 = new((/dims(1),dims(2)/),float)
confi_com2 = t_test2(uw22,auw)
copy_VarMeta(u_comp2(:,:),confi_com2(:,:))
print(max(confi_com2))
print(min(confi_com2))
;;;;;;;;;;;;Significance test of difference field;;;;;;;;;;;

dif_uwnd=new((/dims(1),dims(2)/),float)
dif_uwnd(:,:)=u_comp1(:,:)-u_comp2(:,:)
copy_VarMeta(u_comp2(:,:),dif_uwnd(:,:))
print(max(dif_uwnd))
print(min(dif_uwnd))

s1_ts       =dim_variance_n(uw12(:,:,:),0) ;
s2_ts       =dim_variance_n(uw22(:,:,:),0)  ;
iflag = False
tval_opt = False
confi_dif =new((/dims(1),dims(2)/),float)
confi_dif=ttest(u_comp1,s1_ts,nnumb1,u_comp2,s2_ts,nnumb2,iflag,False)
copy_VarMeta(u_comp1(:,:),confi_dif(:,:))
print(max(confi_dif))
print(min(confi_dif))

wks = gsn_open_wks("eps","/public/home/wangji19/program3/The strong mgc_maxt_300hPa_ep2ek")  

gsn_define_colormap(wks,"cmp_b2r" )

 res=True
res@gsnAddCyclic = False;;;True
;res@gsnMaximize  = True
res@gsnDraw=False
res@gsnFrame=False
;res@pmTickMarkDisplayMode="Always"
res@gsnLeftString=""
res@gsnRightString=""
res@tmXBLabelFontHeightF    = 0.02    ;set the size of x-Axis words
res@tmYLLabelFontHeightF    = 0.02   ;set the size of y-Axis words
; res@gsnPanelLabelBar=False
mpres=res
 mpres@mpOutlineOn=True
mpres@mpFillOn=True;;False
mpres@mpFillColor="gray";;False
mpres@mpPerimOn = 	True
mpres@mpPerimDrawOrder = "PostDraw"
mpres@mpPerimLineThicknessF = 0.1
mpres@mpGeophysicalLineThicknessF=0.1
mpres@mpLimitMode="LatLon"
mpres@mpCenterLonF=90
mpres@mpMinLatF = 0
mpres@mpMaxLatF = 90
mpres@mpMinLonF = 0
mpres@mpMaxLonF = 180
; mpres@pmTickMarkDisplayMode="Always"
;;;;;;;;;;;xy经纬度间隔
mpres@tmXTOn=False
mpres@tmYROn=False
mpres@tmXBMode="Explicit"
mpres@tmXBValues=(/0,60,120,180/)
mpres@tmXBLabels=(/"0~S~o~N~","60~S~o~N~E","120~S~o~N~E","180~S~o~N~"/);
mpres@tmYLMode="Explicit"
mpres@tmYLValues=(/90,60,30,0/);;;"40~S~o~N~N",,"10~S~o~N~N"
mpres@tmYLLabels=(/"90~S~o~N~N","60~S~o~N~N","30~S~o~N~N","0~S~o~N~"/);
mpres@gsnMinorLonSpacing=10
mpres@gsnMajorLonSpacing=30
mpres@gsnMinorLatSpacing=10
; mpres@gsnMajorLatSpacing=20
map=gsn_csm_map(wks,mpres)

resc=res
resc@gsnAddCyclic = False;;;
resc@cnFillOn=True
resc@cnLinesOn=False
; resc@cnLineDashPattern=0
; resc@cnLineThicknessF=1
resc@cnLineLabelsOn=False
;;;;;;;;;;;;;;;;;;;;;;色标设置
resc@lbLabelBarOn         =  True
resc@cnInfoLabelOn          = False;;
resc@lbLabelFontHeightF   =0.01
resc@pmLabelBarOrthogonalPosF  =0.1;latlon0.1;;nh;0.01
resc@pmLabelBarHeightF    =0.1
resc@pmLabelBarWidthF     =0.6
  ; resc@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
  ; resc@cnMinLevelValF       = -30;;;-40;;-3.0;;-3.0;;              ; set min contour level
  ; resc@cnMaxLevelValF       =  30;;;40  ;;3.0;;3.0;;            ; set max contour level
  ; resc@cnLevelSpacingF      = 6;;5  ;;0.5;;0.5;
   resc@cnLevelSelectionMode="ExplicitLevels"
   resc@cnLevels=(/-5,-1,-0.5,-0.2,0,0.2,0.5,1,5/);;-30,-20,-10,6,10,20,,30
   resc@cnFillColors=(/5,9,14,20,0,0,42,48,53,58/);3,9,58,,64
   
plot1=gsn_csm_contour(wks,u_comp1,resc)
; plot1=gsn_csm_contour(wks,u_comp2,resc)
; plot1=gsn_csm_contour(wks,dif_uwnd,resc)
overlay(map,plot1)

;;;;;;;;;;;;;;;显著区域打点 ------------------
res2 = res
res2@cnFillOn                = True
res2@cnLinesOn               = False
res2@cnLevelSelectionMode   = "ExplicitLevels"
res2@cnLevels               = (/0.01/) ; significance level
res2@cnFillPalette          ="GMT_gray"
res2@cnFillColors           =(/2,-1/)
res2@cnFillPatterns         = (/17,-1/);;;;填充形状
res2@cnFillScales           = (/0.5,1/)
res2@cnFillDotSizeF         = 0.001;0.003
; res2@cnFillScales           = 2
res2@cnMonoFillPattern      = False
res2@cnMonoFillScale        = False
res2@cnInfoLabelOn          = False
res2@cnLineLabelsOn         = False
res2@lbLabelBarOn           = False
plot_p=gsn_csm_contour(wks,confi_com1(:,:),res2)
; plot_p=gsn_csm_contour(wks,confi_com2(:,:),res2)
; plot_p=gsn_csm_contour(wks,confi_dif(:,:),res2)

overlay(map,plot_p)

cnres=res
cnres@gsnAddCyclic = False;;; True;; 
cnres@cnLinesOn=True
cnres@cnFillOn=False
; cnres@gsnLeftString=letter(n)+" "+"Day"+" "+point(n)
cnres@cnLevelSelectionMode="AutomaticLevels"
cnres@cnLevelSpacingF=4;;2;;10;;10;;
cnres@cnLineColor="black"
cnres@cnLineThicknessF=1 

cnres@cnLineLabelsOn=True;;False;
cnres@cnLineLabelPerimOn =True
cnres@cnLabelDrawOrder     = "PostDraw"
cnres@cnLineLabelFontHeightF=0.008
cnres@cnLineLabelBackgroundColor = "white"
cnres@cnLineLabelInterval=2
cnres@cnLineLabelDensityF=0.5
cnres@cnLabelMasking=True
cnres@cnLineLabelPerimColor = "white"
cnres@cnLineLabelPerimThicknessF = 0.5
cnres@cnLineLabelPlacementMode="Constant";“Computed”
cnres@cnInfoLabelOn=False

; cnres@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
; cnres@cnMinLevelValF       = -30               ; set min contour level
; cnres@cnMaxLevelValF       = 30             ; set max contour level
; cnres@cnLevelSpacingF      = 6 

   cnres@cnLevelSelectionMode="ExplicitLevels"
   cnres@cnLevels=(/-1,1/)
   ; cnres@cnFillColors=(/3,9,12,17,20,26,34,42,48,53,58,64/);
   

; ; plot=gsn_csm_contour(wks,u_comp1,cnres)
 ; plot=gsn_csm_contour(wks,u_comp2,cnres)
; ; plot=gsn_csm_contour(wks,dif_uwnd,cnres)
; ; plot=ColorNegDashZeroPosContour(plot,"blue","white","red")
; overlay(map,plot)


draw(map)
frame(wks)
end
