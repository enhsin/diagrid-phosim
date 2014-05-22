;;
;; @package phosim
;; @file validation_2B.pro
;; @brief validation task 2B
;;
;; @brief Created by:
;; @author Nathan Todd (Purdue)
;;
;; @brief Modified by:
;; @author John R. Peterson (Purdue)
;;
;; @warning This code is not fully validated
;; and not ready for full release.  Please
;; treat results with caution.
;;

pro validation_2B,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison


  print,'Task 2B'
  !p.multi=[0,5,3]


  !EXCEPT=0
  sigma = fltarr(2,10)
  ellip1 = fltarr(2,10)
  ellip2 = fltarr(2,10)
  centx=fltarr(2,10)
  centy=fltarr(2,10)
  for n=0,1 do begin
     if n eq 0 then high=9
     if n eq 1 then high=4

     for m=0,high do begin

        temp = strcompress('lsst_e_21'+string(n)+string(m)+'_f2_R32_S21_E000.fits.gz',/remove)

        FITS_READ, temp, data, header

        measurepsf,data,rms,e1,e2,medx,medy,flux

        sigma[n,m]=rms
        ellip1(n,m)=e1
        ellip2(n,m)=e2
        centx(n,m)=medx
        centy(n,m)=medy

     endfor
  endfor


;make contour plots

  temp = strarr(15)
  ii=0
  for n=0,1 do begin
     if n eq 0 then high=9 
     if n eq 1 then high=4

     for m=0,high do begin

        temp(ii) = strcompress('lsst_e_21'+string(n)+string(m)+'_f2_R32_S21_E000.fits.gz',/remove)
        ii++
     endfor
  endfor
  data0=mrdfits(temp(0),0,/silent)
  data1=mrdfits(temp(1),0,/silent)
  data2=mrdfits(temp(2),0,/silent)
  data3=mrdfits(temp(3),0,/silent)
  data4=mrdfits(temp(4),0,/silent)
  data5=mrdfits(temp(5),0,/silent)
  data6=mrdfits(temp(6),0,/silent)
  data7=mrdfits(temp(7),0,/silent)
  data8=mrdfits(temp(8),0,/silent)
  data9=mrdfits(temp(9),0,/silent)
  data10=mrdfits(temp(10),0,/silent)
  data11=mrdfits(temp(11),0,/silent)
  data12=mrdfits(temp(12),0,/silent)
  data13=mrdfits(temp(13),0,/silent)
  data14=mrdfits(temp(14),0,/silent)

  image0=data0[1800:2200,1836:2236]
  image1=data1[1800:2200,1836:2236]
  image2=data2[1800:2200,1836:2236]
  image3=data3[1800:2200,1836:2236]
  image4=data4[1800:2200,1836:2236]
  image5=data5[1800:2200,1836:2236]
  image6=data6[1800:2200,1836:2236]
  image7=data7[1800:2200,1836:2236]
  image8=data8[1800:2200,1836:2236]
  image9=data9[1800:2200,1836:2236]
  image10=data10[1800:2200,1836:2236]
  image11=data11[1800:2200,1836:2236]
  image12=data12[1800:2200,1836:2236]
  image13=data13[1800:2200,1836:2236]
  image14=data14[1800:2200,1836:2236]


  sss=size(image0)
  xx=(findgen(sss[1])-sss[1]/2)*0.5
  yy=(findgen(sss[2])-sss[2]/2)*0.5


  contour,alog(image0>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11
  xyouts,-80.0,80.0,'aber+pert '+string(sigma[0,0]*2*0.01*sqrt(2*alog(2)),format='(F4.2)')+'"',/data
  contour,alog(image1>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11
  xyouts,-80.0,80.0,'aber+pert '+string(sigma[0,1]*2*0.01*sqrt(2*alog(2)),format='(F4.2)')+'"',/data
  contour,alog(image2>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11
  xyouts,-80.0,80.0,'aber+pert '+string(sigma[0,2]*2*0.01*sqrt(2*alog(2)),format='(F4.2)')+'"',/data
  contour,alog(image3>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11
  xyouts,-80.0,80.0,'aber+pert '+string(sigma[0,3]*2*0.01*sqrt(2*alog(2)),format='(F4.2)')+'"',/data
  contour,alog(image4>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11
  xyouts,-80.0,80.0,'aber+pert '+string(sigma[0,4]*2*0.01*sqrt(2*alog(2)),format='(F4.2)')+'"',/data
  contour,alog(image5>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11
  xyouts,-80.0,80.0,'aber+pert '+string(sigma[0,5]*2*0.01*sqrt(2*alog(2)),format='(F4.2)')+'"',/data
  contour,alog(image6>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11
  xyouts,-80.0,80.0,'aber+pert '+string(sigma[0,6]*2*0.01*sqrt(2*alog(2)),format='(F4.2)')+'"',/data
  contour,alog(image7>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11
  xyouts,-80.0,80.0,'aber+pert '+string(sigma[0,7]*2*0.01*sqrt(2*alog(2)),format='(F4.2)')+'"',/data
  contour,alog(image8>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11
  xyouts,-80.0,80.0,'aber+pert '+string(sigma[0,8]*2*0.01*sqrt(2*alog(2)),format='(F4.2)')+'"',/data
  contour,alog(image9>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11
  xyouts,-80.0,80.0,'aber+pert '+string(sigma[0,9]*2*0.01*sqrt(2*alog(2)),format='(F4.2)')+'"',/data

  contour,alog(image10>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11
  xyouts,-80.0,80.0,'aber+track '+string(sigma[1,0]*2*0.01*sqrt(2*alog(2)),format='(F4.2)')+'"',/data
  contour,alog(image11>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11
  xyouts,-80.0,80.0,'aber+track '+string(sigma[1,1]*2*0.01*sqrt(2*alog(2)),format='(F4.2)')+'"',/data
  contour,alog(image12>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11
  xyouts,-80.0,80.0,'aber       '+string(sigma[1,2]*2*0.01*sqrt(2*alog(2)),format='(F4.2)')+'"',/data
  contour,alog(image13>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11
  xyouts,-80.0,80.0,'ab+dm see '+string(sigma[1,3]*2*0.01*sqrt(2*alog(2)),format='(F4.2)')+'"',/data
  contour,alog(image14>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11
  xyouts,-80.0,80.0,'ab+det    '+string(sigma[1,4]*2*0.01*sqrt(2*alog(2)),format='(F4.2)')+'"',/data


  sigmaopt=sigma*0.5
  centxopt=centx*0.5
  centyopt=centy*0.5
  ellip1opt=ellip1
  ellip2opt=ellip2
  ss='Optics Perturbations & Tracking'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 2B; '+vers
  xyouts,0.7,0.98,ss,/normal

  tolerance_low(nnn,0)=0.0
  tolerance_high(nnn,0)=0.381

  tolerance_low(nnn,1)=0.0
  tolerance_low(nnn,2)=0.0
  tolerance_low(nnn,3)=0.0
  tolerance_low(nnn,4)=0.0
  tolerance_high(nnn,1)=0.097
  tolerance_high(nnn,2)=0.090
  tolerance_high(nnn,3)=0.259
  tolerance_high(nnn,4)=0.246

  tolerance_high(nnn,5)=0.21
  tolerance_high(nnn,6)=0.21
  tolerance_high(nnn,7)=0.21
  tolerance_high(nnn,8)=0.21
  tolerance_high(nnn,9)=0.21

  tolerance_high(nnn,10)=1.0
  tolerance_high(nnn,11)=1.0
  tolerance_high(nnn,12)=1.0
  tolerance_high(nnn,13)=1.0
  tolerance_high(nnn,14)=1.0

  name(nnn,0)='Total non-atm PSF size'
  name(nnn,1)='  Optics Design'
  name(nnn,2)='  Internal Seeing'
  name(nnn,3)='  Perturbations/Misalignments'
  name(nnn,4)='  Charge Diffusion'

  name(nnn,5)='PSF Sqrt(ellip)*FWHM'
  name(nnn,6)='  Optics Design'
  name(nnn,7)='  Internal Seeing'
  name(nnn,8)='  Perturbations/Misalignments'
  name(nnn,9)='  Charge Diffusion'

  name(nnn,10)='PSF Centroid'
  name(nnn,11)='  Optics Design'
  name(nnn,12)='  Internal Seeing'
  name(nnn,13)='  Perturbations/Misalignments'
  name(nnn,14)='  Charge Diffusion'

  task(nnn,0)='2B Instrumental PSF Properties'

  unit(nnn,0)=' arcsec'
  comparison(nnn,0)='LSST Design'

  unit(nnn,1)=' arcsec'
  comparison(nnn,1)='LSST Design'
  unit(nnn,2)=' arcsec'
  comparison(nnn,2)='LSST Design'
  unit(nnn,3)=' arcsec'
  comparison(nnn,3)='LSST Design'
  unit(nnn,4)=' arcsec'
  comparison(nnn,4)='LSST Design'

  unit(nnn,5)=' arcsec'
  comparison(nnn,5)='LSST Design'
  unit(nnn,6)=' arcsec'
  comparison(nnn,6)='LSST Design'
  unit(nnn,7)=' arcsec'
  comparison(nnn,7)='LSST Design'
  unit(nnn,8)=' arcsec'
  comparison(nnn,8)='LSST Design'
  unit(nnn,9)=' arcsec'
  comparison(nnn,9)='LSST Design'
  unit(nnn,10)=' arcsec'
  comparison(nnn,10)='LSST Design'
  unit(nnn,11)=' arcsec'
  comparison(nnn,11)='LSST Design'
  unit(nnn,12)=' arcsec'
  comparison(nnn,12)='LSST Design'
  unit(nnn,13)=' arcsec'
  comparison(nnn,13)='LSST Design'
  unit(nnn,14)=' arcsec'
  comparison(nnn,14)='LSST Design'

   value(nnn,0)=sqrt((mean(sigmaopt(0,*))*2.35)^2+(mean(sigmaopt(1,0:1))*2.35)^2+(sigmaopt(1,4)*2.35)^2-1.0*(sigmaopt(1,2)*2.35)^2)/50.0
   value(nnn,1)=sigmaopt(1,2)*2.35/50.0
   value(nnn,2)=sqrt((mean(sigmaopt(1,3))*2.35)^2-(sigmaopt(1,2)*2.35)^2)/50.0
   value(nnn,3)=sqrt((mean(sigmaopt(0,*)*2.35))^2+(mean(sigmaopt(1,0:1))*2.35)^2-2.0*(sigmaopt(1,2)*2.35)^2)/50.0
   value(nnn,4)=sigmaopt(1,4)*2.35/50.0

   value(nnn,6)=(mean(ellip1opt(1,2)^4+ellip2opt(1,2)^4))^(0.25)*value(nnn,1)
   value(nnn,7)=(mean(ellip1opt(1,3)^4+ellip2opt(1,3)^4))^(0.25)*value(nnn,2)
   value(nnn,8)=(mean(ellip1opt(0,*)^4+ellip2opt(0,*)^4))^(0.25)*value(nnn,3)
   value(nnn,9)=(mean(ellip1opt(1,4)^4+ellip2opt(1,4)^4))^(0.25)*value(nnn,4)
   value(nnn,5)=sqrt(total((value(nnn,6:9))^2))

   medx=centxopt(1,2) & medy=centyopt(1,2)
   value(nnn,11)=sqrt((centxopt(1,2)-medx)^2+(centyopt(1,2)-medy)^2)/50.0
   value(nnn,12)=sqrt(mean((centxopt(1,3)-medx)^2+(centyopt(1,3)-medy)^2))/50.0
   value(nnn,13)=sqrt(mean((centxopt(0,*)-medx)^2+(centyopt(0,*)-medy)^2))/50.0
   value(nnn,14)=sqrt(mean(centxopt(1,4)^2+centyopt(1,4)^2))/50.0
   value(nnn,10)=sqrt(total((value(nnn,11:14))^2))


END
