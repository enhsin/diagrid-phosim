;;
;; @package phosim
;; @file validation_1H.pro
;; @brief validation task 1H
;;
;; @brief Created by:
;; @author En-Hsin Peng (Purdue)
;;
;; @brief Modified by:
;;
;; @warning This code is not fully validated
;; and not ready for full release.  Please
;; treat results with caution.
;;

;FUNCTION calculate_structure_function, img, pixelsize
;  dim=size(img)
;  N=dim(1)
;  Dstr=DBLARR(N)
;  npx=DBLARR(N)
;  maxidx=FLOOR(ALOG(N)/ALOG(1.42)-1)
;  NLoops=100000
;  FOR xidx=0L, maxidx DO BEGIN
;    dx=ROUND(1.42^xidx)
;    FOR yidx=0L, maxidx DO BEGIN
;      dy=ROUND(1.42^yidx)
;      d=ROUND(SQRT(1.42^(2*xidx)+1.42^(2*yidx)))
;      FOR l=0L, Nloops-1 DO BEGIN
;        ix=FLOOR(RANDOMU(S)*(N-dx))
;        iy=FLOOR(RANDOMU(S)*(N-dy))
;        Dstr[d]=Dstr[d]+(img[ix,iy]-img[ix+dx,iy+dy])^2
;        npx[d]=npx[d]+1.0
;      ENDFOR
;    ENDFOR
;  ENDFOR
;  posid=WHERE(npx GT 0)
;  r=FINDGEN(N)*pixelsize
;  RETURN, [[r[posid]], [Dstr[posid]/npx[posid]]]
;END

FUNCTION calculate_structure_function, img, pixelsize
  seed=100L
  dim=size(img)
  N=dim(1)
  Dstr=DBLARR(N)
  npx=DBLARR(N)
  maxr=FLOOR(ALOG(N)/ALOG(sqrt(2.0)))-1
  NLoops=1000
  FOR rd=0L,maxr DO BEGIN
     d=floor(sqrt(2)^(rd))
     FOR l=0L, Nloops-1 DO BEGIN
        ix=FLOOR(RANDOMU(seed)*N)
        iy=FLOOR(RANDOMU(seed)*N)
        phi=randomu(seed)*2*!Pi
        ixn=floor(ix+d*cos(phi))
        iyn=floor(iy+d*sin(phi))
        if (ixn ge 0 and ixn lt N and iyn ge 0 and iyn lt N) then begin
           Dstr[d]=Dstr[d]+(img[ix,iy]-img[ixn,iyn])^2
           npx[d]=npx[d]+1.0
        endif
     endfor
  endfor
  posid=WHERE(npx GT 0)
  r=FINDGEN(N)*pixelsize
  RETURN, [[r[posid]], [Dstr[posid]/npx[posid]]]
end



PRO validation_1H,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison
  print,'Task 1H'
  
  obsids=['1991','1992','1993','1994','1995','1996','1997','1998','1999']
  mm=5
  mm_all=N_ELEMENTS(obsids)
  mm=mm_all
  savec=fltarr(mm)
  ;Adams & Skrutskie 2MASS
  adams_fig4=[[0.09, 0.10, 0.10, 0.11, 0.12, 0.13, 0.13, 0.14, 0.15, 0.17, 0.18, 0.19, 0.21, 0.23, 0.25, 0.25, 0.27, 0.29, 0.29, 0.32, 0.32, 0.35, 0.38, 0.40, 0.42, 0.45, 0.47, 0.52, 0.55, 0.56, 0.63, 0.72, 0.82, 0.92, 0.97, 1.22, 1.33, 1.80, 2.08, 2.70, 3.70, 5.39], [1.64, 1.77, 1.70, 1.75, 1.83, 1.90, 1.79, 1.97, 2.15, 2.15, 2.25, 2.39, 2.48, 2.51, 2.39, 2.73, 2.45, 2.91, 3.16, 2.42, 2.67, 2.57, 2.70, 2.98, 2.80, 3.70, 2.87, 3.05, 2.84, 3.89, 2.84, 3.20, 2.87, 3.09, 3.09, 3.89, 3.66, 5.46, 4.89, 7.94, 11.71, 44.40]] 

  ;Ivezic et al. 2007, AJ 134, 3, 973.
  ivezic_fig22_cir=[[0.42,0.84,1.26,1.67,2.09],[0.019,0.025,0.030,0.035,0.043]]
  ivezic_fig22_tri=[[0.33,0.66,1.00,1.34,1.66,1.99,2.33],[0.007,0.014,0.022,0.030,0.039,0.047,0.053]]
  ivezic_fig22_sqr=[[0.41,0.84,1.26,1.68,2.09],[0.012,0.020,0.026,0.032,0.039]]

  !p.multi=[0,2,2]
  ;airglow screen
  print, ' airglow'
  strfun=DBLARR(mm,2)
  FOR m=0, mm-1 DO BEGIN
    obsid=obsids[m]
    image=mrdfits('airglowscreen_'+obsid+'.fits',0,/silent)
    pixelsize=15.0/3600  ; in degrees
    D_airglow=calculate_structure_function(image,pixelsize)
    D_airglow[*,1]=SQRT(D_airglow[*,1])*100
    IF m EQ 0 THEN begin
       PLOT, adams_fig4[*,0], adams_fig4[*,1], LINESTYLE=0, title='Airglow', ytitle='sqrt[D(r)] (%)', xtitle='r (deg)', xr=[0.01,10], yr=[0.2,100], /xstyle,/ystyle, /XLOG, /YLOG
       xyouts,1.0,70.0,'Airglow Screens',/data
    endif

    OPLOT, D_airglow[*,0], D_airglow[*,1], LINESTYLE=1,color=250
    strfun[m,0]=INTERPOL(D_airglow[*,1],D_airglow[*,0],0.1)
    strfun[m,1]=INTERPOL(D_airglow[*,1],D_airglow[*,0],1.0)
  ENDFOR
  legend,linestyle=[1,0],['phosim','Adams & Skrutskie 1996'],box=0,charsize=0.8,/top,/left,color=[250,0]

  k=0
  value(nnn,k)=MEAN(strfun[*,0])
  model=INTERPOL(adams_fig4[*,1],adams_fig4[*,0],0.1)
  tolerance_high(nnn,k)=model+0.1
   tolerance_low(nnn,k)=model-0.1
  unit(nnn,k)=' %'
  name(nnn,k)='Airglow Fluctuation at 0.1 deg'
  comparison(nnn,k)='Adams & Skrutskie'

  k=1
  value(nnn,k)=MEAN(strfun[*,1])
  model=INTERPOL(adams_fig4[*,1],adams_fig4[*,0],1.0)
  tolerance_high(nnn,k)=model+1.0
  tolerance_low(nnn,k)=model-1.0
  unit(nnn,k)=' %'
  name(nnn,k)='Airglow Fluctuation at 1.0 deg'
  comparison(nnn,k)='Adams & Skrutskie'

;  print,max(image),min(image)
  ximage=findgen(1024)*15./3600.
  yimage=findgen(1024)*15./3600.
  loadct,3
  contour,image,ximage,yimage,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (degrees)',ytitle='Position (degrees)'
  loadct,39

  ;cloud screen
;  print, ' cloud'
  N=1024
  v=20.0  ; wind velocity (m/s) ~ 3-15 deg/min
  ap=2.5  ; SDSS aperture size (m)   
  expt=55 ; exposure time (s)
  pixelscale=256.0 ; (cm)
  wx=ROUND(SQRT(ap*(ap+v*expt))*100/pixelscale)    
  wy=wx
  totalcloudmean=0.0
  FOR m=0, mm_all-1 DO BEGIN
    obsid=obsids[m]
    SPAWN, "grep cloudmean atmosphere_"+obsid+".pars| awk '{print $3}'", cloudmean
    totalcloudmean=totalcloudmean+cloudmean[0]+cloudmean[1]
  ENDFOR
  totalcloudmean=totalcloudmean/mm_all
;  print, '  average cloudmean', totalcloudmean
 
  strfun=DBLARR(mm_all,2)
  FOR m=0, mm_all-1 DO BEGIN
    obsid=obsids[m]
    SPAWN, "grep height atmosphere_"+obsid+".pars| awk '{print $3}'", height
    SPAWN, "grep cloudvary atmosphere_"+obsid+".pars| awk '{print $3}'", cloudvary
    
    pixelsize=pixelscale/((height[1]+height[2])/2*1e5)/!PI*180  ; in degrees

    image=mrdfits('cloudscreen_'+obsid+'_1.fits',0,/silent)
    image1=[[image,image,image],[image,image,image],[image,image,image]]
    image=mrdfits('cloudscreen_'+obsid+'_2.fits',0,/silent)
    image2=[[image,image,image],[image,image,image],[image,image,image]]
    simage=SMOOTH(cloudvary[0]*image1+cloudvary[1]*image2,[wx,wy])
    D_cloud=calculate_structure_function(simage[N:2*N-1,N:2*N-1],pixelsize)
    D_cloud[*,1]=SQRT(D_cloud[*,1])
    
    IF m EQ 0 THEN begin
       PLOT, ivezic_fig22_cir[*,0], totalcloudmean/1.3*ivezic_fig22_cir[*,1], PSYM=1, ytitle='sqrt[D(r)] (mag)', xtitle='r (deg)', xr=[0.01,10], yr=[1e-3,0.2], /xstyle,/ystyle, /XLOG, /YLOG
       xyouts,1.0,0.15,'Cloud Screens',/data
    endif

    OPLOT, D_cloud[*,0], D_cloud[*,1], LINESTYLE=1,color=250

    strfun[m,0]=INTERPOL(D_cloud[*,1],D_cloud[*,0],ivezic_fig22_cir[0,0])
    strfun[m,1]=INTERPOL(D_cloud[*,1],D_cloud[*,0],ivezic_fig22_cir[4,0])
  ENDFOR
  legend,linestyle=[1],['phosim'],box=0,charsize=0.8,color=[250],position=[0.012,0.11]
  legend,psym=[1],['Ivezic et al. 2007'],box=0,charsize=0.8,position=[0.012,0.08]

  k=2
  value(nnn,k)=MEAN(strfun[*,1])
  model=ivezic_fig22_cir[4,1]
  tolerance_high(nnn,k)=model*2
  tolerance_low(nnn,k)=model*0.1
  unit(nnn,k)=' mag'
  name(nnn,k)='Cloud SF at 2 deg'
  comparison(nnn,k)='Ivezic et al. 2007'
;  print, tolerance_high(nnn,k), value(nnn,k), model, tolerance_low(nnn,k)

  cloudimage=10.^(-0.4*(cloudvary[0]*image1+cloudvary[1]*image2+totalcloudmean))
  ximage=findgen(1024)*2.56
  yimage=findgen(1024)*2.56
  loadct,0
  contour,cloudimage(0:1023,0:1023),ximage,yimage,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (m)',ytitle='Position (m)'
  loadct,39

  ;turbulence screen
  print, ' turbulence'
  ARCSEC=0.000004841369
  natmospherefile=7
  totalseeing=0.67
  wavelength=0.62
  alt=89.0
  zenith=(90.0-alt)*!pi/180
  zenithfactor=(1/cos(zenith))^0.6
  wavelengthFactor=(wavelength/0.5)^(-0.2)
  seeing=(totalseeing+1e-6)*zenithfactor*wavelengthFactor
  strfun=DBLARR(mm,2)
  FOR m=0, mm-1 DO BEGIN
    obsid=obsids[m]
    SPAWN, "grep seeing atmosphere_"+obsid+".pars | grep -v total | awk '{print $3}'", seefactor
    SPAWN, "grep outer atmosphere_"+obsid+".pars | awk '{print $3}'", outerscales
    seefactor=seefactor*zenithfactor
    density_norm=DBLARR(natmospherefile)
    see_norm=DBLARR(natmospherefile)
    F=5 ; 5x larger than the fine screen
    N=1024*F
    image=DBLARR(N,N)
;    mx=64 & my=64 & nm=64
;    cx=256 & cy=256 & nc=4

    mx=0 & my=0 & nm=128
    cx=0 & cy=0 & nc=16
    lx=0 & ly=0 & nl=2

    outerScaleCorrection=0.0
    FOR i=0L, natmospherefile-1 DO BEGIN
      imagep=mrdfits(strcompress('atmospherescreen_'+obsid+'_'+string(i)+'_largep.fits',/remove),0,/silent) ; 512cm
      imagec=mrdfits(strcompress('atmospherescreen_'+obsid+'_'+string(i)+'_coarsep.fits',/remove),0,headerc,/silent) ; 64cm
      imagem=mrdfits(strcompress('atmospherescreen_'+obsid+'_'+string(i)+'_mediump.fits',/remove),0,/silent) ; 8cm
      imagef=mrdfits(strcompress('atmospherescreen_'+obsid+'_'+string(i)+'_finep.fits',/remove),0,/silent)   ; 1cm
      imagetmp=mrdfits(strcompress('atmospherescreen_'+obsid+'_'+string(i)+'_coarsex.fits',/remove),0,headerx,/silent)
      density_norm[i]=FXPAR(headerc,'NORM')
      see_norm[i]=FXPAR(headerx,'NORM')
;      print,'a',see_norm[i],density_norm[i]
      outerScaleCorrection=outerScaleCorrection+(seefactor[i]/(zenithfactor*(totalseeing/2.35+1e-6)))^(5./3.)*density_norm[i]/1.0e-6
      image1=REBIN(imagem[mx:mx+nm*F-1,my:my+nm*F-1],N,N)+REBIN(imagec[cx:cx+nc*F-1,cy:cy+nc*F-1],N,N)+REBIN(imagep[lx:lx+nl*F-1,ly:ly+nl*F-1],N,N)
      FOR j=0L, F-1 DO BEGIN
        FOR q=0L, F-1 DO BEGIN
          image1[j*1024:(j+1)*1024-1,q*1024:(q+1)*1024-1]=image1[j*1024:(j+1)*1024-1,q*1024:(q+1)*1024-1]+imagef
        ENDFOR
      ENDFOR
      image=image+image1*density_norm[i]*seefactor[i]/(totalseeing*zenithfactor/2.35)
   ENDFOR
;    outerScaleCorrection=3e-3
;    outerscaleCorrection=4.0
;    print,'b',outerScaleCorrection
;    outerScaleCorrection=1.0/(outerScaleCorrection);
    image=image*sqrt(0.0229*(0.98*(1e-4*wavelength)/(seeing*ARCSEC))^(-5./3.))*4.0e8
    pixelsize=1/100.0  ; in meters
    D_turb=calculate_structure_function(image,pixelsize)
    if m eq 0 then turbimage1=image
    if m eq 1 then turbimage2=image
    if m eq 2 then turbimage3=image


    IF m EQ 0 THEN BEGIN
      rr=0.01*EXP(ALOG(30.0/0.01)/20.0)^(FINDGEN(20)+1)
      r0=202140*wavelength/1e6/totalseeing
;      print,'r0 ',r0
      y1=6.88*(rr/r0)^(5./3.)
      plot,rr,y1, LINESTYLE=0, ytitle='D(r)', xtitle='r (m)', xr=[0.02,40], yr=[1e-2,1e5], /xstyle,/ystyle, /XLOG, /YLOG
      xyouts,2.0,3e4,'Turbulence Screens',/data

   endif
    outerscale=total(seefactor^(5./3.)*outerscales)/total(seefactor^(5./3.))
    outerscale=10.0^(1.0+float(m)/(float(mm-1))*1.0)
;    print,outerscale
    k0=(2*!Pi)/outerscale
    u=findgen(100000)/100000.*1e3
    y2=fltarr(N_elements(rr))
    for i=0L,N_elements(rr)-1 do begin
       integrand=u*(1-beselj(u,0))/((u*u+k0*k0*rr(i)*rr(i))^(11./6.))
       f=int_tabulated(u,integrand)
       y2(i)=0.173*(outerscale/r0)^(5./3.)*5./3.*(k0*rr(i))^(5./3.)*f
    endfor
    ;Lucke & Young 2007, AO 256.
    oplot,rr,y2,linestyle=2
    OPLOT, D_turb[*,0], D_turb[*,1],color=250,linestyle=1
;    print,'c',D_turb[0,1]
    savec(m)=D_turb[0,1]
;    print,'x',D_turb[0,0]
    strfun[m,0]=INTERPOL(D_turb[*,1],D_turb[*,0],0.02)
    strfun[m,1]=INTERPOL(D_turb[*,1],D_turb[*,0],0.1)
  ENDFOR
  legend,linestyle=[1,2,0],['phosim','von Karman (Lucke & Young 2007)','Kolmogorov (Fried 1965)'],box=0,charsize=0.8,/top,/left,color=[250,0,0]
  k=3
  value(nnn,k)=MEAN(strfun[*,1])
  model=6.88*(0.1/r0)^(5.0/3)
  tolerance_high(nnn,k)=model*1.2
  tolerance_low(nnn,k)=model*0.8
  unit(nnn,k)=' rad!U2!N'
  name(nnn,k)='Phase Structure Function at 0.1 m'
  comparison(nnn,k)='Kolmogorov Spectrum'
;  print, tolerance_high(nnn,k), value(nnn,k), tolerance_low(nnn,k)
  k=4
  value(nnn,k)=ALOG(MEAN(strfun[*,1])/MEAN(strfun[*,0]))/ALOG(0.1/0.02)
  tolerance_high(nnn,k)=5.0/3*1.1
  tolerance_low(nnn,k)=5.0/3*0.9
  unit(nnn,k)=' '
  name(nnn,k)='Phase Structure Function Slope Index'
  comparison(nnn,k)='Kolmogorov Spectrum'
;  print, tolerance_high(nnn,k), value(nnn,k), tolerance_low(nnn,k)

  sss=size(turbimage1)
  xxx=findgen(sss(1))*pixelsize
  yyy=findgen(sss(2))*pixelsize
  loadct,5
  contour,turbimage1,xxx,yyy,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (m)',ytitle='Position (m)'
  contour,turbimage2,xxx,yyy,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (m)',ytitle='Position (m)'
  contour,turbimage3,xxx,yyy,nlevels=11,/fill,/xstyle,/ystyle,xtitle='Position (m)',ytitle='Position (m)'
  loadct,39

  ss='Structure Functions'
  XYOUTS,0.1,0.98,ss,/normal
  ss='Validation Task 1H; '+vers
  XYOUTS,0.7,0.98,ss,/normal

  task(nnn,0)='1H Structure Functions'

;  print,median(savec)
END
