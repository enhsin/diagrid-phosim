;;
;; @package phosim
;; @file validation_3C.pro
;; @brief validation task 3C
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

pro validation_3C,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 3C'
  !p.multi=[0,1,2]


  tolerance_low(nnn,0)=2.9*0.0
  tolerance_high(nnn,0)=2.9*2.0

  buffer=200

  N=14L
  flux=fltarr(N)
  variance=fltarr(N)
  autovert=fltarr(N)
  autohori=fltarr(N)

  for i=0L,N-1 do begin

     if i lt 10 then file='320'+string(i,format='(I1)')+'_f2_R22_S11'
     if i ge 10 then file='32'+string(i,format='(I2)')+'_f2_R22_S11'
     j=i+20
     if j lt 10 then file2='320'+string(j,format='(I1)')+'_f2_R22_S11'
     if j ge 10 then file2='32'+string(j,format='(I2)')+'_f2_R22_S11'

     filename='lsst_e_'+file+'_E000.fits.gz'
     data=mrdfits(filename,0,header,/silent)
     filename='lsst_e_'+file2+'_E000.fits.gz'
     data2=mrdfits(filename,0,header,/silent)

     print,total(data),total(data2)
     print,file
     print,file2


     dataa=data+data2
     datab=data-data2

     resulta=moment(dataa((1820-buffer):(1820+buffer),(2036-buffer):(2036+buffer)))
     resultb=moment(datab((1820-buffer):(1820+buffer),(2036-buffer):(2036+buffer)))

     flux(i)=resulta(0)/2.0
     variance(i)=resultb(1)/2.0

     shortbuffer=buffer-1
     aa=0.0 & bb=0.0 & cc=0.0 
     for k=1820-shortbuffer,1820+shortbuffer do begin
        for l=2036-shortbuffer,2036+shortbuffer do begin
           aa=aa+(datab(k,l)*datab(k,l+1))
           bb=bb+(datab(k,l)*datab(k+1,l))
           cc=cc+1.0
        endfor
     endfor
     autovert(i)=aa/cc/variance(i)
     autohori(i)=bb/cc/variance(i)


;     value(nnn,i)=sqrt(clippedresult(1)-noise^2)*100.0
;     name(nnn,i*2)='PRNU at '+string(wave,format='(I3)')+' nm'
;     unit(nnn,i*2)='%'
;     comparison(nnn,i*2)='Prototype Devices at BNL'



   endfor

  varianceerr=flux/sqrt((2.0*buffer)^2)*sqrt(2.0)

  plot,flux,variance,psym=4,xr=[1e1,1e6],yr=[1e1,1e6],xtitle='Signal',ytitle='Variance',/xstyle,/ystyle,/xlog,/ylog

  errplot,flux,variance-varianceerr,variance+varianceerr,width=0.0

  xxx=10.^(1.0+5.0*findgen(1000)/1000.)
  oplot,xxx,xxx

  bnlflux=[20000.0,40000.0,60000.0,80000.,100000.0,120000.,140000.,160000.]/140000.*100000.
  bnlvariance=[20000.0,40000.0,60000.0*32.0/33.0,80000.0*40./45.,100000.*50./55.,100000.*40./55.,100000.*5./55.,100000.*5./55.]/140000.*100000.


  oplot,bnlflux,bnlvariance,psym=5,color=250



  nu=variance/flux
  pcterr=varianceerr/flux

  plot,flux,nu,psym=4,/xlog,xr=[1e1,1e6],xtitle='Signal',ytitle='Variance Non-Linearity',yr=[0,1.2],/xstyle

  errplot,flux,nu-pcterr,nu+pcterr,width=0.0


   bnlnu=bnlvariance/bnlflux

  oplot,bnlflux,bnlnu,psym=5,color=250


  oplot,xxx,1.0+fltarr(N_elements(xxx)),linestyle=2

;  oplot,xxx,1.0-0.04*(xxx/1e5),linestyle=2

  ss='Variance vs. Signal on Difference Images'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 3C; '+vers
  xyouts,0.7,0.98,ss,/normal

  plot,flux,autovert,xtitle='Signal',ytitle='Column Autocorrelation',xr=[1e1,1e6],/xlog,psym=4,yr=[0,0.1]
  plot,flux,autohori,xtitle='Signal',ytitle='Row Autocorrelation',xr=[1e1,1e6],/xlog,psym=4,yr=[0,0.1]

  ss='Variance vs. Signal on Difference Images'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 3C; '+vers
  xyouts,0.7,0.98,ss,/normal

  task(nnn,0)='3C Variance vs. Signal'
  print,flux
  print,variance


END
