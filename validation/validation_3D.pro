;;
;; @package phosim
;; @file validation_3D.pro
;; @brief validation task 3D
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

pro validation_3D,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 3D'
  !p.multi=[0,1,2]


  tolerance_low(nnn,0)=2.9*0.0
  tolerance_high(nnn,0)=2.9*2.0


  ss='Variance vs. Signal on Difference Images'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 3D; '+vers
  xyouts,0.7,0.98,ss,/normal


  nn=20
  w=fltarr(20)
  x=fltarr(20)
  y=fltarr(20)
  z=fltarr(20)
  v=fltarr(20,nn)
  er=fltarr(20)

  data=mrdfits('lsst_e_3300_f2_R22_S11_E000.fits.gz')

  mag=14.5
  for ii=0,19 do begin
     mag=mag+0.5
     for jj=0,nn-1 do begin

        x0=296+180*jj
        y0=3742-180*ii

        datas=data((x0-50):(x0+50),(y0-50):(y0+50))

        measurepsf,datas,a,b,c,d,e,f

        x(ii)=mag
        y(ii)=y(ii)+a*2.35*0.2
        print,a
        z(ii)=z(ii)+max(datas)
        w(ii)=w(ii)+f
        v(ii,jj)=a*2.35*0.2

        for kk=0,19 do begin
           if nn gt 1 then result=moment(v(kk,*))
           if nn gt 1 then er(kk)=sqrt(result(1))
        endfor


        g=where(z/float(nn) lt 5000.0)
        y1=median(y(g))/float(nn)
        xx=findgen(1000)/1000.*1e5
        y1=0.724
        in=0.04/5.0

        plot,x,y/float(nn),psym=-4,/xstyle,/ystyle,yr=[0.72,0.73]
        if nn gt 1 then errplot,x,y/float(nn)-er/float(nn),y/float(nn)+er/float(nn),width=0.0

        mm=findgen(1000)/1000.*10.0+15.
        pp=2143.*10.^(0.4*(20.0-mm))
        oplot,mm,y1*(1.0+in*(pp/1e5)),linestyle=2
        oplot,mm,y1*(1.0+0.0*(pp/1e5)),linestyle=2


        plot,z/float(nn),y/float(nn),psym=-4,/ystyle,yr=[0.72,0.73]
        if nn gt 1 then errplot,z/float(nn),y/float(nn)-er/float(nn),y/float(nn)+er/float(nn),width=0.0

        oplot,xx,y1*(1.0+in*(xx/1e5)),linestyle=2
        oplot,xx,y1*(1.0+0.0*(xx/1e5)),linestyle=2

     endfor
  endfor


  task(nnn,0)='3D PSF size vs. intensity'
  print,flux
  print,variance


END
