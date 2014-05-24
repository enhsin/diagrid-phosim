;;
;; @package phosim
;; @file validation_4D.pro
;; @brief validation task 4D
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

pro validation_4D,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 4D'
  !p.multi=[0,2,1]

  ll=0.01+findgen(11)*8.0

  data=mrdfits('lsst_e_4300_f2_R10_S00_E000.fits.gz',0,/silent)
  image=rebin(data(0:1999,0:3999),200,400)
  sss=size(image)
  xx=(findgen(sss(1))-(sss(1)-2)/2.0)*100.0
  yy=(findgen(sss(2))-(sss(2)-2)/2.0)*100.0
  contour,image,xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11,levels=ll
  data1=rebin(data,500,509)*64.0

  data=mrdfits('lsst_e_4301_f2_R10_S00_E000.fits.gz',0,/silent)
  image=rebin(data(0:1999,0:3999),200,400)
  sss=size(image)
  xx=(findgen(sss(1))-(sss(1)-2)/2.0)*100.0
  yy=(findgen(sss(2))-(sss(2)-2)/2.0)*100.0
  contour,image,xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11,levels=ll
  data2=rebin(data,500,509)*64.0

  ss='Background Optimization Accuracy'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 4D; '+vers
  xyouts,0.7,0.98,ss,/normal

  good=where(data1 ge 1 or data2 ge 1)
  chi2=total((data1(good)-data2(good))^2/((data1(good)+data2(good))))/float(N_elements(good))
  value(nnn,0)=chi2
  tolerance_low(nnn,0)=0.0
  tolerance_high(nnn,0)=2.0

  name(nnn,0)='Pixel Comparison'
  task(nnn,0)='4D Corner chip Back Opt on/off'
  unit(nnn,0)=' (!4V!3!U2!N/dof)'
  comparison(nnn,0)='Exact Calculation'

END
