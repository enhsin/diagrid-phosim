;;
;; @package phosim
;; @file validation_4A.pro
;; @brief validation task 4A
;;
;; @brief Created by:
;; @author John R. Peterson (Purdue)
;;
;; @brief Modified by:
;;
;; @warning This code is not fully validated
;; and not ready for full release.  Please
;; treat results with caution.
;;

pro validation_4A,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 4A'
  !p.multi=0

  readcol,'mag_dist',mm,point,diffuse,/silent
  plot,[-1,40],[1.0,1e11],/nodata,/ylog,/xstyle,/ystyle,xtitle='Magnitude',ytitle='Number'
;number of source per sq. degree
  oplot,mm,point/(!Pi*2.1*2.1),color=250
  oplot,mm,diffuse/(!pi*2.1*2.1),color=50
  back=fltarr(N_elements(mm))+1e-12
  g1=where(mm eq 22)
  back(g1)=3600.*3600.
  oplot,mm,back,color=200

  calcphotons=(point+diffuse)/(!pi*2.1*2.1)*(13.3/60.)^2*10.^(0.4*(34.25-mm))

  legend,linestyle=[0,0,0],['Stars mag!U-1!N deg!U-2!N','Galaxies mag!U-1!N deg!U-2!N','Equiv. Arcsec Background'],color=[250,50,200],/top,/right

  legend,psym=[4,5],['Photons mag!U-1!N chip!U-1!N','Time mag!U-1!N chip!U-1!N (seconds)'],color=[20,0],/top,/left




  readcol,'speed_out_0',cc,bb,dd,aa,format=('A,A,A,A'),/silent
  time=float(dd)*3600.0*24.0+float(strmid(aa,0,2))*3600.0+float(strmid(aa,3,2))*60.0+float(strmid(aa,6,2))
  ttime=time

  spawn,'/bin/ls -lrtT lsst_e_4001*.fits.gz | awk ''{print $7,$8}'' > temp.txt'
  readcol,'temp.txt',dd,aa,format=('A,A'),/silent
  time=float(dd)*3600.0*24.0+float(strmid(aa,0,2))*3600.0+float(strmid(aa,3,2))*60.0+float(strmid(aa,6,2))
  ttime=[ttime,time]

  readcol,'speed_out_1',cc,bb,dd,aa,format=('A,A,A,A'),/silent
  time=float(dd)*3600.0*24.0+float(strmid(aa,0,2))*3600.0+float(strmid(aa,3,2))*60.0+float(strmid(aa,6,2))
  ttime=[ttime,time]

  spawn,'/bin/ls -lrtT *_4000*.fits.gz | awk ''{print $7,$8}'' > temp.txt'
  readcol,'temp.txt',dd,aa,format=('A,A'),/silent
  time=float(dd)*3600.0*24.0+float(strmid(aa,0,2))*3600.0+float(strmid(aa,3,2))*60.0+float(strmid(aa,6,2))
  ttime=[ttime,time]

  readcol,'speed_out_2',cc,bb,dd,aa,format=('A,A,A,A'),/silent
  time=float(dd)*3600.0*24.0+float(strmid(aa,0,2))*3600.0+float(strmid(aa,3,2))*60.0+float(strmid(aa,6,2))
  ttime=[ttime,time]


  parsrays=(ttime(1)-ttime(0))

  spawn,'more speed_out | grep Dark | awk ''{print $4,$10}'' | sed -e "s/,//" | sed -e "s/,//" | sed -e "s/,//" > speed_background'
  readcol,'speed_background',bb,aa,format='(F,F)',/silent

  background=total(aa)/2.0
  backgroundrate=total(bb)/total(aa)
  backphotons=fltarr(N_elements(mm))+1e-12
  g1=where(mm eq 22)
  backphotons(g1)=total(bb)/2.0

  parstrim=(ttime(3)-ttime(2))-(ttime(20)-ttime(19))
  e2adc=0.5*((ttime(19)-ttime(3))+(ttime(36)-ttime(20)))
  raysray=(ttime(20)-ttime(19))

  spawn,'more speed_out | grep Astro | grep -v m=3 | awk ''{print $4,$10}'' | sed -e "s/,//" | sed -e "s/,//" | sed -e "s/,//" > speed_background'

  readcol,'speed_background',bb,aa,format='(F,F)',/silent
  ray=total(aa)/2.0+background

  spawn,'more speed_out | grep Astro | grep -v m=3 | awk ''{print $2,$4,$10}'' | sed -e "s/,//" | sed -e "s/,//" | sed -e "s/,//" | sed -e "s/m=//" > speed_background'

  readcol,'speed_background',cc,bb,aa,format='(F,F)',/silent
  astrophotons=fltarr(N_elements(mm))
  astrotime=fltarr(N_elements(mm))
  astrophotons(cc)=bb
  astrotime(cc)=aa
  g1=where(astrophotons eq 0)
  astrophotons(g1)=calcphotons(g1)
  oplot,mm,astrophotons+backphotons,color=20,psym=4


  spawn,'more speed_out | grep Astro | grep -E ''m=29|m=28|m=27|m=26|m=25|m=24|m=23|m=22|m=21|m=20|m=19|m=18|m=17|m=16''| awk ''{print $4,$10}'' | sed -e "s/,//" | sed -e "s/,//" | sed -e "s/,//" > speed_background'

  readcol,'speed_background',bb,aa,format='(F,F)',/silent
  rayunsat=total(aa)/2.0
  unsatrate=total(bb)/total(aa)

  spawn,'more speed_out | grep Astro | grep m=15 | awk ''{print $4,$10}'' | sed -e "s/,//" | sed -e "s/,//" | sed -e "s/,//" > speed_background'
  readcol,'speed_background',bb,aa,format='(F,F)',/silent
  ray15=total(aa)/2.0
  ray15rate=total(bb)/total(aa)

  spawn,'more speed_out | grep Astro | grep m=14 | awk ''{print $4,$10}'' | sed -e "s/,//" | sed -e "s/,//" | sed -e "s/,//" > speed_background'
  readcol,'speed_background',bb,aa,format='(F,F)',/silent
  ray14=total(aa)/2.0
  ray14rate=total(bb)/total(aa)

  spawn,'more speed_out | grep Astro | grep m=13 | awk ''{print $4,$10}'' | sed -e "s/,//" | sed -e "s/,//" | sed -e "s/,//" > speed_background'
  readcol,'speed_background',bb,aa,format='(F,F)',/silent
  ray13=total(aa)/2.0
  ray13rate=total(bb)/total(aa)

  spawn,'more speed_out | grep Astro | grep m=12 | awk ''{print $4,$10}'' | sed -e "s/,//" | sed -e "s/,//" | sed -e "s/,//" > speed_background'
  readcol,'speed_background',bb,aa,format='(F,F)',/silent
  ray12=total(aa)/2.0
  ray12rate=total(bb)/total(aa)

  xyouts,0.1,0.70,'Atm+Opt Setup '+string(round((parsrays-(raysray-ray))/6.0)/10.0,format='(F6.1)'),/norm
  xyouts,0.1,0.67,'Trim            '+string(round((parstrim-(parsrays-(raysray-ray)))/6.0)/10.0,format='(F6.1)'),/norm
  xyouts,0.1,0.64,'Raytrace Setup '+string(round((raysray-ray)/6.0)/10.0,format='(F6.1)'),/norm
  xyouts,0.1,0.61,'Raytrace       '+string(round((ray/6.0))/10.0,format='(F6.1)'),/norm
  xyouts,0.1,0.58,'   Unsaturated     '+string(round((rayunsat)/6.0)/10.0,format='(F6.1)'),/norm
  xyouts,0.1,0.55,'   Background     '+string(round((background)/6.0)/10.0,format='(F6.1)'),/norm
  xyouts,0.1,0.52,'   Saturated       '+string(round((ray-rayunsat-background)/6.0)/10.0,format='(F6.1)'),/norm
  xyouts,0.1,0.49,'e2adc          '+string(round((e2adc/6.0))/10.0,format='(F6.1)'),/norm
  xyouts,0.1,0.45,'Total           '+string(round(((parsrays-(raysray-ray))+((parstrim-(parsrays-(raysray-ray))))+((raysray-ray))+ray+e2adc)/6.0)/10.0,format='(F6.1)'),/norm
  

  xyouts,0.1,0.40,' Unsaturated rate   '+string(round(unsatrate/1000.0),format='(I5)')+' kphoton/s',/norm
  xyouts,0.1,0.37,' Background speedup '+string(round(backgroundrate/unsatrate),format='(I5)'),/norm
  xyouts,0.1,0.34,' m=15 speedup       '+string(round(ray15rate/unsatrate),format='(I5)'),/norm
  xyouts,0.1,0.31,' m=14 speedup       '+string(round(ray14rate/unsatrate),format='(I5)'),/norm
  xyouts,0.1,0.28,' m=13 speedup       '+string(round(ray13rate/unsatrate),format='(I5)'),/norm
  xyouts,0.1,0.25,' m=12 speedup       '+string(round(ray12rate/unsatrate),format='(I5)'),/norm


  value(nnn,0)=round(((parsrays-(raysray-ray))+((parstrim-(parsrays-(raysray-ray))))+((raysray-ray))+ray+e2adc)/6.0)/10.0
  value(nnn,1)=round(unsatrate/1000.0)
  value(nnn,2)=round(backgroundrate/unsatrate)
  value(nnn,3)=round(ray12rate/unsatrate)


  rate=fltarr(N_elements(mm))+unsatrate
  g1=where(mm le 12)
  rate(g1)=(ray12rate/unsatrate+ray12rate/ray13rate*(12-mm(g1)))*unsatrate
  g1=where(mm eq 13)
  rate(g1)=ray13rate
  g1=where(mm eq 14)
  rate(g1)=ray14rate
  g1=where(mm eq 15)
  rate(g1)=ray15rate


  oplot,mm,astrotime+backphotons/backgroundrate,psym=5


  ss='Speed Test'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 4A; '+vers
  xyouts,0.7,0.98,ss,/normal

  tolerance_low(nnn,0)=0.0
  tolerance_high(nnn,0)=0.0
  tolerance_low(nnn,0)=0.0
  tolerance_high(nnn,0)=300.0

  tolerance_low(nnn,1)=100
  tolerance_high(nnn,1)=1e30

  tolerance_low(nnn,2)=1.0
  tolerance_high(nnn,2)=1e30

  tolerance_low(nnn,3)=1.0
  tolerance_high(nnn,3)=1e30
  name(nnn,0)='Total simulation time'
  name(nnn,1)='Unsaturation simulation rate'
  name(nnn,2)='Dark Sky Speed Factor'
  name(nnn,3)='m=12 Speed Factor'
  unit(nnn,0)=' minutes'
  comparison(nnn,0)='N/A'

  unit(nnn,1)=' kphotons/s'
  comparison(nnn,1)='N/A'

  unit(nnn,2)=' '
  comparison(nnn,2)='N/A'

  unit(nnn,3)=' '
  comparison(nnn,3)='N/A'
  task(nnn,0)='4A Chip w/ Stars and Galaxies'

END
