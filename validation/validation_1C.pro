;;
;; @package phosim
;; @file validation_1C.pro
;; @brief validation task 1C
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

pro validation_1C,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  for qqq=0,2 do begin

     if qqq eq 0 then print,'Task 1C'
     if qqq eq 1 then print,'Task 1E'
     if qqq eq 2 then print,'Task 1F'

     lsigmaatm=fltarr(5,3,5,5,2)
     le1atm=fltarr(5,3,5,5,2)
     le2atm=fltarr(5,3,5,5,2)
     lposx=fltarr(5,3,5,5,2)
     lposy=fltarr(5,3,5,5,2)
     lmedx=fltarr(5,3,5,5,2)
     lmedy=fltarr(5,3,5,5,2)
     maxnn=4

     for oo=0,1 do begin
        for nn=0,maxnn do begin
           print,' Observation '+string(nn,format='(I1)')+' '+string(oo,format='(I1)')

           for ii=0,49 do begin

              if ii lt 25 then begin
                 xx=floor(ii/5)
                 yy=ii mod 5
                 lposx(nn,0,xx,yy,oo)=-1.41+0.705*xx
                 lposy(nn,0,xx,yy,oo)=-1.41+0.705*yy

                 if (xx eq 0 and yy eq 0) or $
                    (xx eq 4 and yy eq 0) or $
                    (xx eq 0 and yy eq 4) or $
                    (xx eq 4 and yy eq 4) then goto,skipraft

                 if qqq eq 0 then filename='lsst_e_120'+string(nn,format='(I1)')+'_f2_R'+string(xx,format='(I1)')+string(yy,format='(I1)')+'_S11_E00'+string(oo,format='(I1)')+'.fits.gz'
                 if qqq eq 1 then filename='lsst_e_140'+string(nn,format='(I1)')+'_f2_R'+string(xx,format='(I1)')+string(yy,format='(I1)')+'_S11_E00'+string(oo,format='(I1)')+'.fits.gz'
                 if qqq eq 2 then filename='lsst_e_150'+string(nn,format='(I1)')+'_f2_R'+string(xx,format='(I1)')+string(yy,format='(I1)')+'_S11_E00'+string(oo,format='(I1)')+'.fits.gz'

                 data=mrdfits(filename,0,/silent)
                 if xx eq 2 and yy eq 2 then begin
                    for kk=0,2 do begin
                       for ll=0,2 do begin
                          xl=2000.+(kk-1)*1200.-600.
                          xh=2000.+(kk)*1200.-600.
                          yl=2036.+(ll-1)*1200.-600.
                          yh=2036.+(ll)*1200.-600.
                          lposx(nn,2,kk,ll,oo)=-0.033+0.033*kk
                          lposy(nn,2,kk,ll,oo)=-0.033+0.033*ll
                          subdata=data(xl:xh,yl:yh)
                          measurepsf,subdata,rms,e1,e2,medx,medy,flux
                          lsigmaatm(nn,2,kk,ll,oo)=rms
                          le1atm(nn,2,kk,ll,oo)=e1
                          le2atm(nn,2,kk,ll,oo)=e2
                          lmedx(nn,2,kk,ll,oo)=medx
                          lmedy(nn,2,kk,ll,oo)=medy
                       endfor
                    endfor
                 endif else begin
                    measurepsf,data,rms,e1,e2,medx,medy,flux
                    lsigmaatm(nn,0,xx,yy,oo)=rms
                    le1atm(nn,0,xx,yy,oo)=e1
                    le2atm(nn,0,xx,yy,oo)=e2
                    lmedx(nn,0,xx,yy,oo)=medx
                    lmedy(nn,0,xx,yy,oo)=medy
                 endelse


              endif else begin
                 xx=floor((ii-25)/5)
                 yy=(ii-25) mod 5
                 lposx(nn,1,xx,yy,oo)=-0.47+0.235*xx
                 lposy(nn,1,xx,yy,oo)=-0.47+0.235*yy

                 if xx eq 0 or xx eq 4 or yy eq 0 or yy eq 4 then goto,skipraft
                 if xx eq 2 and yy eq 2 then goto,skipraft

                 if qqq eq 0 then filename='lsst_e_120'+string(nn,format='(I1)')+'_f2_R22_S'+string(xx-1,format='(I1)')+string(yy-1,format='(I1)')+'_E00'+string(oo,format='(I1)')+'.fits.gz'
                 if qqq eq 1 then filename='lsst_e_140'+string(nn,format='(I1)')+'_f2_R22_S'+string(xx-1,format='(I1)')+string(yy-1,format='(I1)')+'_E00'+string(oo,format='(I1)')+'.fits.gz'
                 if qqq eq 2 then filename='lsst_e_150'+string(nn,format='(I1)')+'_f2_R22_S'+string(xx-1,format='(I1)')+string(yy-1,format='(I1)')+'_E00'+string(oo,format='(I1)')+'.fits.gz'

                 data=mrdfits(filename,0,/silent)
                 measurepsf,data,rms,e1,e2,medx,medy,flux
                 lsigmaatm(nn,1,xx,yy,oo)=rms
                 le1atm(nn,1,xx,yy,oo)=e1
                 le2atm(nn,1,xx,yy,oo)=e2
                 lmedx(nn,1,xx,yy,oo)=medx
                 lmedy(nn,1,xx,yy,oo)=medy
              endelse



skipraft:

           endfor
        endfor
     endfor


;copy arrays to not get confused about other exposure
     sigmaatm=fltarr(5,3,5,5)
     e1atm=fltarr(5,3,5,5)
     e2atm=fltarr(5,3,5,5)
     posx=fltarr(5,3,5,5)
     posy=fltarr(5,3,5,5)
     for nn=0,maxnn do begin
        for ii=0,4 do begin
           for jj=0,4 do begin
              for mm=0,2 do begin
                 sigmaatm(nn,mm,ii,jj)=lsigmaatm(nn,mm,ii,jj,0)
                 e1atm(nn,mm,ii,jj)=le1atm(nn,mm,ii,jj,0)
                 e2atm(nn,mm,ii,jj)=le2atm(nn,mm,ii,jj,0)
                 posx(nn,mm,ii,jj)=lposx(nn,mm,ii,jj,0)
                 posy(nn,mm,ii,jj)=lposy(nn,mm,ii,jj,0)
              endfor
           endfor
        endfor
     endfor

     for rrr=0,1 do begin
        !p.multi=[0,2,2]

;subtract off common mode ellipticity fit
;NEED TO DO EACH ONE INDIVIDUALLY
        if rrr eq 1 then begin
           medx=median(posx)
           medy=median(posy)
           x=posx-medx
           y=posy-medy
           xs=(max(x)-min(x))/2.0
           ys=(max(y)-min(y))/2.0
           x=x/xs
           y=y/ys
           for nn=0,maxnn do begin
              minchi=1e30
              av=fltarr(12)
              bav=av
              for iter=0,10 do begin
                 for k=0,11 do begin
                    av=bav
                    for aa=-0.1,0.1,0.001 do begin
                       av(k)=aa
                       e1t=av(0)+av(2)*x+av(3)*y+av(6)*x*y+av(7)*x^2+av(8)*y^2
                       e2t=av(1)+av(4)*x+av(5)*y+av(9)*x*y+av(10)*x^2+av(11)*y^2
                       chi=0.0
                       for ii=0,4 do begin
                          for jj=0,4 do begin
                             for mm=0,2 do begin
                                if sigmaatm(nn,mm,ii,jj) ne 0 then begin
                                   chi=chi+(e1t(nn,mm,ii,jj)-e1atm(nn,mm,ii,jj))^2+(e2t(nn,mm,ii,jj)-e2atm(nn,mm,ii,jj))^2
                                endif
                             endfor
                          endfor
                       endfor
                       if chi lt minchi then begin
                          bav(k)=av(k)
                          minchi=chi
                          e1c=e1atm-e1t
                          e2c=e2atm-e2t
                       endif
                    endfor
                 endfor
              endfor
              for ii=0,4 do begin
                 for jj=0,4 do begin
                    for mm=0,2 do begin
                       e1atm(nn,mm,ii,jj)=e1c(nn,mm,ii,jj)
                       e2atm(nn,mm,ii,jj)=e2c(nn,mm,ii,jj)
                    endfor
                 endfor
              endfor

           endfor
        endif


;used to be 0.25
        q1=histogram(sqrt((sigmaatm(where(sigmaatm ne 0))*2.35/10.0)^2-0.00^2),min=0.,bin=0.02,max=1.3)
        xx=findgen(N_elements(q1))*0.02
        plot,xx,q1/total(q1),xtitle='PSF FWHM (arcseconds)',psym=10,yr=[0,0.3],/ystyle,ytitle='Fraction/bin'
        result=moment(sqrt((sigmaatm(where(sigmaatm ne 0))*2.35/10.0)^2-0.00^2))
        ss=string(result(0),format='(F4.2)')+' +/- '+string(sqrt(result(1)),format='(F4.2)')
        xyouts,0.1,0.1,ss,/data
        pass1=abs(result(0)-0.65)/sqrt(result(1))
        value(nnn,6*qqq)=result(0)

        for nn=0,maxnn do begin
           kk=0L & sigmaarr=fltarr(50)
           for ii=0,4 do begin
              for jj=0,4 do begin
                 for mm=0,2 do begin
                    if (sigmaatm(nn,mm,ii,jj) gt 0) then begin
                       sigmaarr(kk)=sigmaatm(nn,mm,ii,jj)
                       kk=kk+1
                    endif
                 endfor
              endfor
           endfor
           sigmaarr=sigmaarr(0:(kk-1))


           q1=histogram(sqrt((sigmaarr*2.35/10.0)^2-0.00^2),min=0.,bin=0.02,max=1.3)
           xx=findgen(N_elements(q1))*0.02
           oplot,xx,q1/total(q1)/5.0,color=(nn+1)*50,psym=10,linestyle=0
        endfor


        xyouts,0.12,0.9,'Input seeing was 0.65" (atmospheric)',/normal
        xyouts,0.12,0.87,'Stars in grid pattern',/normal
        xyouts,0.12,0.84,'for five different atmospheres',/normal


        q2=histogram(sqrt(e1atm(where(sigmaatm ne 0))^2+e2atm(where(sigmaatm ne 0))^2),min=0.,bin=0.01,max=0.15)
        xx=findgen(N_elements(q2))*0.01
        plot,xx,q2/total(q2),xtitle='Ellipticity',psym=10,yr=[0,0.4],/ystyle,ytitle='Fraction/bin'
        result=moment(sqrt(e1atm(where(sigmaatm ne 0))^2+e2atm(where(sigmaatm ne 0))^2))
        ss='phoSim:  '+string(result(0),format='(F5.3)')+' +/- '+string(sqrt(result(1)),format='(F5.3)')
        xyouts,0.07,0.2,ss,/data
        resulti=result

        if rrr eq 1 then filename='subaru_1.txt' else filename='subaru_0.txt'
        readcol,filename,aa,ab,ac,ad,ae,/silent
        ell_su=sqrt(ac^2+ad^2)
        qq=histogram(ell_su,min=0.,bin=0.01,max=0.15)
        xx=findgen(N_elements(qq))*0.01
        oplot,xx,qq/total(qq),psym=10,linestyle=2
        if rrr eq 1 then filename='cfht_1.txt' else filename='cfht_0.txt'
        readcol,filename,aa,ab,ac,ad,ae,/silent
        ell_cf=sqrt(ac^2+ad^2)
        qq=histogram(ell_cf,min=0.,bin=0.01,max=0.15)
        xx=findgen(N_elements(qq))*0.01
        oplot,xx,qq/total(qq),psym=10,linestyle=1

        legend,linestyle=[0,2,1],['phoSim','Subaru','CFHT'],/top,/right
        result=moment(ell_cf)
        ss='CFHT:   '+string(result(0),format='(F5.3)')+' +/- '+string(sqrt(result(1)),format='(F5.3)')
        xyouts,0.07,0.10,ss,/data
        result=moment(ell_su)
        ss='Subaru: '+string(result(0),format='(F5.3)')+' +/- '+string(sqrt(result(1)),format='(F5.3)')
        xyouts,0.07,0.15,ss,/data


        pass2=abs(result(0)-resulti(0))/sqrt(result(1)+resulti(1))
        value(nnn,1+6*qqq)=resulti(0)
        tolerance_low(nnn,1+6*qqq)=0.0
        tolerance_high(nnn,1+6*qqq)=result(0)

        value(nnn,2+6*qqq)=sqrt(resulti(1))
        tolerance_low(nnn,2+6*qqq)=0.0
        tolerance_high(nnn,2+6*qqq)=sqrt(result(1))


;kstwo,sqrt(e1atm(where(sigmaatm ne 0))^2+e2atm(where(sigmaatm ne 0))^2),ell_su,D,prob
;ss='KS p-value: '+string(prob,format='(F5.3)')
;xyouts,0.07,19.0,ss,/data


        for nn=0,maxnn do begin
           kk=0L & earr=fltarr(50)
           for ii=0,4 do begin
              for jj=0,4 do begin
                 for mm=0,2 do begin
                    if (sigmaatm(nn,mm,ii,jj) gt 0) then begin
                       earr(kk)=e1atm(nn,mm,ii,jj)*e1atm(nn,mm,ii,jj)+e2atm(nn,mm,ii,jj)*e2atm(nn,mm,ii,jj)
                       kk=kk+1
                    endif
                 endfor
              endfor
           endfor
           earr=earr(0:(kk-1))
           q2=histogram(sqrt(earr),min=0.,bin=0.01,max=0.15)
           xx=findgen(N_elements(q2))*0.01
           oplot,xx,q2/total(q2)/5.0,color=(nn+1)*50,psym=10,linestyle=0
        endfor



        corre=fltarr(5,20)
        corrs=fltarr(5,20)
        corree=fltarr(5,20)
        corrse=fltarr(5,20)
        correa=fltarr(5,20)
        corrsa=fltarr(5,20)
        corrn=fltarr(5,20)

        for nn=0,maxnn do begin

           for ii=0,4 do begin
              for jj=0,4 do begin
                 for kk=0,4 do begin
                    for ll=0,4 do begin
                       for mm=0,2 do begin

                          if (sigmaatm(nn,mm,ii,jj) gt 0.0 and $
                              sigmaatm(nn,mm,kk,ll) gt 0.0) then begin

                             dy=posy(nn,mm,ii,jj)-posy(nn,mm,kk,ll)
                             dx=posx(nn,mm,ii,jj)-posx(nn,mm,kk,ll)
                             r2=dx*dx+dy*dy
;    cos2phi=(dx*dx-dy*dy)/r2
;    sin2phi=2.0*dx*dy/r2
;    if (r2 eq 0.0) then begin
;        cos2phi=1.0
;        sin2phi=0.0
;    endif
;    epar1=e1atm(ii,jj)*cos2phi+e2atm(ii,jj)*sin2phi
;    epar2=e1atm(kk,ll)*cos2phi+e2atm(kk,ll)*sin2phi
;    eper1=e2atm(ii,jj)*cos2phi-e1atm(ii,jj)*sin2phi
;    eper2=e2atm(kk,ll)*cos2phi-e1atm(kk,ll)*sin2phi
;    corr=epar1*epar2+eper1*eper2
                             corr=e1atm(nn,mm,ii,jj)*e1atm(nn,mm,kk,ll)+e2atm(nn,mm,ii,jj)*e2atm(nn,mm,kk,ll)

                             xi=round(alog10((sqrt(r2))>0.01)*4.0+8.0)
                             corre(nn,xi)=corre(nn,xi)+corr
                             corree(nn,xi)=corree(nn,xi)+corr*corr

                             corr=sigmaatm(nn,mm,ii,jj)*sigmaatm(nn,mm,kk,ll)
                             corrs(nn,xi)=corrs(nn,xi)+corr
                             corrse(nn,xi)=corrse(nn,xi)+corr*corr
                             corrn(nn,xi)=corrn(nn,xi)+1.0
                             correa(nn,xi)=correa(nn,xi)+0.5*(e1atm(nn,mm,ii,jj)*e1atm(nn,mm,ii,jj)+e2atm(nn,mm,ii,jj)*e2atm(nn,mm,ii,jj)+e1atm(nn,mm,kk,ll)*e1atm(nn,mm,kk,ll)+e2atm(nn,mm,kk,ll)*e2atm(nn,mm,kk,ll))
                             corrsa(nn,xi)=corrsa(nn,xi)+0.5*(sigmaatm(nn,mm,ii,jj)+sigmaatm(nn,mm,kk,ll))
                          endif

                       endfor
                    endfor
                 endfor
              endfor
           endfor
        endfor


        for nn=0,maxnn do begin
           meansigma=corrsa(nn,*)/(corrn(nn,*)>1)
           corrse(nn,*)=sqrt(corrse(nn,*)/(corrn(nn,*)>1)-(corrs(nn,*))^2/(corrn(nn,*)>1)^2)
           corrse(nn,*)=corrse(nn,*)/meansigma/meansigma
           corrs(nn,*)=corrs(nn,*)/meansigma/meansigma
           corrs(nn,*)=corrs(nn,*)/(corrn(nn,*)>1)
        endfor

        xx=10.^((-8.0+findgen(20))/4.0)

        corrst=corrs(0,*)+corrs(1,*)+corrs(2,*)+corrs(3,*)+corrs(4,*)
        corrsnt=corrn(0,*)+corrn(1,*)+corrn(2,*)+corrn(3,*)+corrn(4,*)
        corrsst=sqrt((corrs(0,*)-corrst/5.0)^2+(corrs(1,*)-corrst/5.0)^2+(corrs(2,*)-corrst/5.0)^2+(corrs(3,*)-corrst/5.0)^2+(corrs(4,*)-corrst/5.0)^2)/sqrt(5.0)

        gg=where(corrsnt gt 0)
        plot,xx(gg),corrst(gg)/5.0,psym=4,xtitle='Angle (degrees)',ytitle='PSF Size Correlation',xr=[0.01,20.0],/xstyle,yr=[0.95,1.05],/xlog
        errplot,xx(gg),corrst(gg)/5.0-corrsst(gg),corrst(gg)/5.0+corrsst(gg),width=0.0

        if rrr eq 1 then filename='subaru_corr_1.txt' else filename='subaru_corr_0.txt'
        readcol,filename,sa,sb,sc,sd,se,/silent
        oplot,sa,sb,psym=6
        if rrr eq 1 then filename='cfht_corr_1.txt' else filename='cfht_corr_0.txt'
        readcol,filename,sa,sb,sc,sd,se,/silent
        oplot,sa,sb,psym=5

        for ttt=3,35 do begin
           if rrr eq 1 and ttt lt 10 then filename='cfht_corr_10'+string(ttt,format='(I1)')+'_1'+'.txt'
           if rrr eq 1 and ttt ge 10 then filename='cfht_corr_1'+string(ttt,format='(I2)')+'_1'+'.txt'
           if rrr eq 0 and ttt lt 10 then filename='cfht_corr_10'+string(ttt,format='(I1)')+'_0'+'.txt'
           if rrr eq 0 and ttt ge 10 then filename='cfht_corr_1'+string(ttt,format='(I2)')+'_0'+'.txt'
           readcol,filename,sa,sb,sc,sd,se,/silent
           oplot,sa,sb,psym=-5,symsize=0.5,linestyle=1
        endfor



        legend,psym=[4,6,5],['phoSim','Subaru','CFHT'],/top,/right


        minchi=1e30
        for corrl=0.0,10000.0,10.0 do begin
           y=1.0/((xx(gg)-0.01)/corrl+1.0)
           chi=total((y-corrst(gg)/5.0)^2)
           if chi lt minchi then begin
              minchi=chi
              bestcorrl=corrl
           endif
        endfor
        xx=findgen(20000)/1000.0
        oplot,xx,1.0/(xx/bestcorrl+1.0)


;xx=findgen(8)/2.0
        xx=10.^((-8.0+findgen(20))/4.0)

        for nn=0,maxnn do begin
           gg=where(corrn(nn,*) gt 0)
           oplot,xx(gg),corrs(nn,gg),psym=-4,color=50*(nn+1),symsize=0.5
        endfor

        for nn=0,maxnn do begin
           meane=sqrt(correa(nn,*)/(corrn(nn,*)>1))
           corree(nn,*)=sqrt(corree(nn,*)/(corrn(nn,*)>1)-(corre(nn,*))^2/(corrn(nn,*)>1)^2)
           corree(nn,*)=corree(nn,*)/meane/meane
           corre(nn,*)=corre(nn,*)/meane/meane
           corre(nn,*)=corre(nn,*)/(corrn(nn,*)>1)
        endfor


;xx=findgen(8)/2.0
        xx=10.^((-8.0+findgen(20))/4.0)
        corret=corre(0,*)+corre(1,*)+corre(2,*)+corre(3,*)+corre(4,*)
        corrent=corrn(0,*)+corrn(1,*)+corrn(2,*)+corrn(3,*)+corrn(4,*)
        correst=sqrt((corre(0,*)-corret/5.0)^2+(corre(1,*)-corret/5.0)^2+(corre(2,*)-corret/5.0)^2+(corre(3,*)-corret/5.0)^2+(corre(4,*)-corret/5.0)^2)/sqrt(5.0)

        gg=where(corrent gt 0)
        plot,xx(gg),corret(gg)/5.0,psym=4,xtitle='Angle (degrees)',ytitle='Ellipticity Vector Correlation',xr=[0.01,20.0],/xstyle,yr=[-0.5,1.5],/xlog
        errplot,xx(gg),corret(gg)/5.0-correst(gg),corret(gg)/5.0+correst(gg),width=0.0

        if rrr eq 1 then filename='subaru_corr_1.txt' else filename='subaru_corr_0.txt'
        readcol,filename,sa,sb,sc,sd,se,/silent
        oplot,sa,sc,psym=6
        if rrr eq 1 then filename='cfht_corr_1.txt' else filename='cfht_corr_0.txt'
        readcol,filename,sa,sb,sc,sd,se,/silent
        oplot,sa,sc,psym=5

        for ttt=3,35 do begin
           if rrr eq 1 and ttt lt 10 then filename='cfht_corr_10'+string(ttt,format='(I1)')+'_1'+'.txt'
           if rrr eq 1 and ttt ge 10 then filename='cfht_corr_1'+string(ttt,format='(I2)')+'_1'+'.txt'
           if rrr eq 0 and ttt lt 10 then filename='cfht_corr_10'+string(ttt,format='(I1)')+'_0'+'.txt'
           if rrr eq 0 and ttt ge 10 then filename='cfht_corr_1'+string(ttt,format='(I2)')+'_0'+'.txt'
           readcol,filename,sa,sb,sc,sd,se,/silent
           oplot,sa,sc,psym=-5,symsize=0.5,linestyle=1
        endfor

        legend,psym=[4,6,5],['phoSim','Subaru','CFHT'],/top,/right

        minchi=1e30
        for corrl=0.0,100.0,0.1 do begin
           y=1.0/((xx(gg)-0.01)/corrl+1.0)
           chi=total((y-corret(gg)/5.0)^2)
           if chi lt minchi then begin
              minchi=chi
              bestcorrl=corrl
           endif
        endfor
        xx=findgen(20000)/1000.0
        oplot,xx,1.0/(xx/bestcorrl+1.0)

        xx=10.^((-8.0+findgen(20))/4.0)
        for nn=0,maxnn do begin
           gg=where(corrn(nn,*) gt 0)
           oplot,xx(gg),corre(nn,gg),psym=-4,color=50*(nn+1),symsize=0.5
        endfor




        plots,findgen(100.)/100.*(3.0-0.03)+0.03,fltarr(100)+0.5,linestyle=1
        plots,fltarr(100)+0.03,findgen(100)/100.*0.2+0.4,linestyle=0
        plots,fltarr(100)+3.0,findgen(100)/100.*0.2+0.4,linestyle=0

        if qqq eq 0 then ss='Atmospheric PSF'
        if qqq eq 1 then ss='High Airmass Atmospheric PSF'
        if qqq eq 2 then ss='Optics PSF'
        if rrr eq 1 then ss=ss+' with Common Mode Subtracted'
        xyouts,0.1,0.98,ss,/normal
        if qqq eq 0 then ss='Validation Task 1C; '+vers
        if qqq eq 1 then ss='Validation Task 1E; '+vers
        if qqq eq 2 then ss='Validation Task 1F; '+vers
        xyouts,0.7,0.98,ss,/normal

        gg=where(sigmaatm gt 0.0)



        pass3=abs(alog10(bestcorrl)-alog10(0.3))/1.0

        value(nnn,3+6*qqq)=bestcorrl
        if qqq eq 0 then begin
           tolerance_low(nnn,3+6*qqq)=0.03
           tolerance_high(nnn,3+6*qqq)=3.03
        endif
        if qqq ge 0 then begin
           tolerance_low(nnn,3+6*qqq)=0.03
           tolerance_high(nnn,3+6*qqq)=1e30
        endif


;value(2,0)=(pass1+pass2+pass3)/3.0


;-------------------------------
;BONUS TASK:  DIFFERENTIAL ATMOSPHERE ASTROMETRY
;-------------------------------
        !p.multi=[0,1,2]

        daposx=fltarr(5,3,5,5)
        daposy=fltarr(5,3,5,5)

        for nn=0,maxnn do begin
           for ii=0,4 do begin
              for jj=0,4 do begin
                 for mm=0,2 do begin
                    daposx(nn,mm,ii,jj)=lmedx(nn,mm,ii,jj,1)-lmedx(nn,mm,ii,jj,0)
                    daposy(nn,mm,ii,jj)=lmedy(nn,mm,ii,jj,1)-lmedy(nn,mm,ii,jj,0)
                 endfor
              endfor
           endfor
        endfor


        corrda=fltarr(5,20)
        corrn=fltarr(5,20)

        for nn=0,maxnn do begin

           for ii=0,4 do begin
              for jj=0,4 do begin
                 for kk=0,4 do begin
                    for ll=0,4 do begin
                       for mm=0,2 do begin

                          if (sigmaatm(nn,mm,ii,jj) gt 0.0 and $
                              sigmaatm(nn,mm,kk,ll) gt 0.0) then begin

                             dy=posy(nn,mm,ii,jj)-posy(nn,mm,kk,ll)
                             dx=posx(nn,mm,ii,jj)-posx(nn,mm,kk,ll)
                             r2=dx*dx+dy*dy

                             xi=round(alog10((sqrt(r2))>0.01)*4.0+8.0)
                             
                             corr=sqrt((daposx(nn,mm,ii,jj)-daposx(nn,mm,kk,ll))^2+$
                                       (daposy(nn,mm,ii,jj)-daposy(nn,mm,kk,ll))^2)
                             
                             corrda(nn,xi)=corrda(nn,xi)+corr*0.1*1000.0
                             corrn(nn,xi)=corrn(nn,xi)+1.0

                          endif

                       endfor
                    endfor
                 endfor
              endfor
           endfor
        endfor


        xx=10.^((-8.0+findgen(20))/4.0)

        corrdat=corrda(0,*)+corrda(1,*)+corrda(2,*)+corrda(3,*)+corrda(4,*)
        corrnt=corrn(0,*)+corrn(1,*)+corrn(2,*)+corrn(3,*)+corrn(4,*)

        corrdate=sqrt((corrda(0,*)/(corrn(0,*)>1)-corrdat/(corrnt>1))^2+$
                      (corrda(1,*)/(corrn(1,*)>1)-corrdat/(corrnt>1))^2+$
                      (corrda(2,*)/(corrn(2,*)>1)-corrdat/(corrnt>1))^2+$
                      (corrda(3,*)/(corrn(3,*)>1)-corrdat/(corrnt>1))^2+$
                      (corrda(4,*)/(corrn(4,*)>1)-corrdat/(corrnt>1))^2)/sqrt(5.0)

        gg=where(corrdat gt 0)
        plot,xx(gg),corrdat(gg)/(corrnt(gg)>1),psym=4,xtitle='Angle (degrees)',ytitle='d(Exp Astr Shift)  (milliarcsecs)',xr=[0.01,4.0],/xstyle,yr=[0.0,50.0],/xlog
        errplot,xx(gg),corrdat(gg)/(corrnt(gg)>1)-corrdate(gg),corrdat(gg)/(corrnt(gg)>1)+corrdate(gg),width=0.0

        for nn=0,maxnn do begin
           gg=where(corrn(nn,*) gt 0)
           oplot,xx(gg),corrda(nn,gg)/(corrn(nn,gg)>1),psym=4,color=50*(nn+1)
        endfor

        for xicounter=1,7,2 do begin
           counter=0L
           for nn=0,maxnn do begin
              for ii=0,4 do begin
                 for jj=0,4 do begin
                    for kk=0,4 do begin
                       for ll=0,4 do begin
                          for mm=0,2 do begin

                             if (sigmaatm(nn,mm,ii,jj) gt 0.0 and $
                                 sigmaatm(nn,mm,kk,ll) gt 0.0) then begin

                                dy=posy(nn,mm,ii,jj)-posy(nn,mm,kk,ll)
                                dx=posx(nn,mm,ii,jj)-posx(nn,mm,kk,ll)
                                r2=dx*dx+dy*dy
                                xi=round(alog10((sqrt(r2))>0.01)*4.0+8.0)
                                if xi eq xicounter or xi eq xicounter+1 then begin

                                   corr=sqrt((daposx(nn,mm,ii,jj)-daposx(nn,mm,kk,ll))^2+$
                                             (daposy(nn,mm,ii,jj)-daposy(nn,mm,kk,ll))^2)

                                   if counter eq 0 then shar=[corr]*0.1*1000.0
                                   if counter gt 0 then shar=[shar,corr*0.1*1000.0]
                                   counter=counter+1L
                                endif

                             endif

                          endfor
                       endfor
                    endfor
                 endfor
              endfor
           endfor

           qq=histogram(shar,min=0,bin=5,max=100)
           xxx=findgen(N_elements(qq))*5.0
           if xicounter eq 1 then plot,xxx,qq/total(qq),psym=10,yr=[0,0.5],xtitle='d(Exposure Astrometric shift) (milliarcsecs)',ytitle='Fraction/bin'
           if xicounter ne 1 then oplot,xxx,qq/total(qq),color=(xicounter-1)*40,psym=10
           result=moment(shar)
           ss=string(result(0),format='(F5.1)')+' +/- '+string(sqrt(result(1)),format='(F5.1)')
           xyouts,0.7,0.3-xicounter/2.0*0.03,ss,color=(xicounter-1)*40.0,/normal
           if xicounter eq 1 then value(nnn,4+6*qqq)=result(0)
           if xicounter eq 1 then value(nnn,5+6*qqq)=sqrt(result(1))

        endfor

        legend,linestyle=[0,0,0,0],['<3''','3'' to 6''','6'' to 20''','>20'''],color=[0,80,160,240],/right

        ss='Differential Atmospheric Astrometry'
        xyouts,0.1,0.98,ss,/normal
        if qqq eq 0 then ss='Validation Task 1C; '+vers
        if qqq eq 1 then ss='Validation Task 1E; '+vers
        if qqq eq 2 then ss='Validation Task 1F; '+vers
        xyouts,0.7,0.98,ss,/normal

     endfor

  endfor


  name(nnn,0)='Atm X=1 PSF size'
  name(nnn,1)='  Average Ellipticity'
  name(nnn,2)='  Stdev Ellipticity'
  name(nnn,3)='  Ellipticity Decorrelation Length'
  name(nnn,4)='  Avg Diff Astrometry, 3 arcmin'
  name(nnn,5)='  Stdev Diff Astrometry, 3 arcmin'
  name(nnn,6)='Atm X=1.4 PSF size'
  name(nnn,7)='  Average Ellipticity'
  name(nnn,8)='  Stdev Ellipticity'
  name(nnn,9)='  Ellipticity Decorrelation Length'
  name(nnn,10)='  Avg Diff Astrometry, 3 arcmin'
  name(nnn,11)='  Stdev Diff Astrometry, 3 arcmin'

  name(nnn,12)='Atm+Opt PSF size'
  name(nnn,13)='  Average Ellipticity'
  name(nnn,14)='  Stdev Ellipticity'
  name(nnn,15)='  Ellipticity Decorrelation Length'
  name(nnn,16)='  Avg Diff Astrometry, 3 arcmin'
  name(nnn,17)='  Stdev Diff Astrometry, 3 arcmin'

  tolerance_low(nnn,0)=0.0
  tolerance_high(nnn,0)=2.0
  tolerance_high(nnn,4)=10.0
  tolerance_high(nnn,5)=6.0
  tolerance_high(nnn,10)=10.0
  tolerance_high(nnn,11)=6.0
  tolerance_high(nnn,16)=10.0
  tolerance_high(nnn,17)=6.0
  tolerance_low(nnn,0)=0.65-0.05
  tolerance_high(nnn,0)=0.65+0.05

  tolerance_low(nnn,6)=0.65-0.05
  tolerance_high(nnn,6)=0.65+0.05

  tolerance_low(nnn,12)=sqrt(0.65^2+0.39^2)-0.05
  tolerance_high(nnn,12)=sqrt(0.65^2+0.39^2)+0.05

  unit(nnn,0)=' arcsec'
  unit(nnn,1)=' '
  unit(nnn,2)=' '
  unit(nnn,3)=' degrees'
  unit(nnn,4)=' milliarcsec'
  unit(nnn,5)=' milliarcsec'

  unit(nnn,6)=' arcsec'
  unit(nnn,7)=' '
  unit(nnn,8)=' '
  unit(nnn,9)=' degrees'
  unit(nnn,10)=' milliarcsec'
  unit(nnn,11)=' milliarcsec'

  unit(nnn,12)=' arcsec'
  unit(nnn,13)=' '
  unit(nnn,14)=' '
  unit(nnn,15)=' degrees'
  unit(nnn,16)=' milliarcsec'
  unit(nnn,17)=' milliarcsec'
  comparison(nnn,0)='Input value'
  comparison(nnn,1)='CFHT (some opt included)'
  comparison(nnn,2)='CFHT (some opt included)'
  comparison(nnn,3)='Estimate (of order degree)'
  comparison(nnn,4)='Monet worst case estimate for 15s'
  comparison(nnn,5)='Monet !4r!3 of !4r!3 Estimate'

  comparison(nnn,6)='Input value'
  comparison(nnn,7)='CFHT (some opt included)'
  comparison(nnn,8)='CFHT (some opt included)'
  comparison(nnn,9)='Estimate (of order degree)'
  comparison(nnn,10)='Monet worst case estimate for 15s'
  comparison(nnn,11)='Monet !4r!3 of !4r!3 Estimate'

  comparison(nnn,12)='Input value+LSST Design'
  tolerance_high(nnn,13)=0.1
  tolerance_high(nnn,14)=0.1
  comparison(nnn,13)='LSST Design'
  comparison(nnn,14)='LSST Design'
  comparison(nnn,15)='Estimate (should be ~degree)'
  comparison(nnn,16)='Monet worst case estimate for 15s'
  comparison(nnn,17)='Monet !4r!3 of !4r!3 Estimate'
  task(nnn,0)='1C Grids of Stars'

END
