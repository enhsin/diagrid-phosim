;;
;; @package phosim
;; @file validation_1G.pro
;; @brief validation task 1G
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

FUNCTION star_centroid, image
  dim=SIZE(image)
  maximg=MAX(image,I) 
  IX=I MOD dim(1)
  IY=I/dim(1)
  xlo=IX-50
  xhi=IX+50
  IF xlo LT 0 THEN xlo=0
  IF xhi GE dim(1) THEN xhi=dim(1)-1 
  ylo=IY-50
  yhi=IY+50
  IF ylo LT 0 THEN ylo=0
  IF yhi GE dim(2) THEN yhi=dim(2)-1 
  subdata=image[xlo:xhi,ylo:yhi]
  measurepsf,subdata,rms,e1,e2,medx,medy,flux  
  medx=medx+(xlo+xhi)/2.0
  medy=medy+(ylo+yhi)/2.0
  RETURN, [medx,medy]
END

FUNCTION airRefraction, wavelength, temperature, pressure, water_pressure
  n=64.328+29498.1/(146-1/wavelength/wavelength)+255.4/(41-1/wavelength/wavelength)
  n=n*pressure*(1+(1.049-0.0157*temperature)*1e-6*pressure)/720.883/(1+0.003661*temperature)
  n=n-((0.0624-0.000680/wavelength/wavelength)/(1+0.003661*temperature)*water_pressure)
  n=1e-6*n+1
  RETURN, n
END

PRO validation_1G,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison
  print,'Task 1G'
  !p.multi=0
  chipids=['R02_S11','R20_S11','R24_S11','R42_S11']
  alt=[45,60,75,89]
  zen=90-alt
  dx=DBLARR(4,4)
  dy=DBLARR(4,4)
  FOR j=0, 3 DO BEGIN
    FOR i=0, 3 DO BEGIN
      image=mrdfits(strcompress('lsst_e_18'+string(alt[j])+'_f1_'+chipids[i]+'_E000.fits.gz',/remove),0,/silent)
      xy_red=star_centroid(image)
      image=mrdfits(strcompress('lsst_e_19'+string(alt[j])+'_f1_'+chipids[i]+'_E000.fits.gz',/remove),0,/silent)
      xy_blue=star_centroid(image)
      dx[i,j]=xy_red[0]-xy_blue[0]
      dy[i,j]=xy_red[1]-xy_blue[1]
    ENDFOR
  ENDFOR   
  dr=SQRT(dx^2+dy^2)*0.2  ; in arcsec
  PLOT, zen, dr[0,*], PSYM=1, ytitle='separation (")', xtitle='zenith angle (deg)', yr=[0,0.63], xr=[-0.8,53],/xstyle,/ystyle
  OPLOT, zen+1.5, dr[1,*], PSYM=2
  OPLOT, ABS(zen-1.5), dr[2,*], PSYM=4
  OPLOT, zen, dr[3,*], PSYM=5

  water_pressure=8
  pressure=520
  temperature=20
  z=FINDGEN(54)
  model=TAN(z*!pi/180)*206265*(airRefraction(0.42,temperature,pressure,water_pressure)-airRefraction(0.52,temperature,pressure,water_pressure))
  OPLOT, z, model, LINESTYLE=0
  XYOUTS, 2, 0.6, 'wavelengths: 420nm, 520nm'
  XYOUTS, 2, 0.57, 'line: Filippenko 1982'

  ss='Differential Chromatic Refraction'
  XYOUTS,0.1,0.98,ss,/normal
  ss='Validation Task 1G; '+vers
  XYOUTS,0.7,0.98,ss,/normal

  k=0
  FOR j=0, 2, 2 DO BEGIN
    FOR i=0, 2 DO BEGIN
      value(nnn,k)=dr[i,j]
      IF i EQ 0 THEN z=zen[j]
      IF i EQ 1 THEN z=zen[j]+1.5
      IF i EQ 2 THEN z=ABS(zen[j]-1.5)
      model=TAN(z*!pi/180)*206265*(airRefraction(0.42,temperature,pressure,water_pressure)-airRefraction(0.52,temperature,pressure,water_pressure)) 
      tolerance_high(nnn,k)=model*1.1
      tolerance_low(nnn,k)=model*0.9
      unit(nnn,k)=' arcsec'
      name(nnn,k)='Differential Refraction'
      comparison(nnn,k)='Filippenko 1982'
      k=k+1
    ENDFOR
  ENDFOR

  task(nnn,0)='1G Differential Chromatic Refraction'
END






