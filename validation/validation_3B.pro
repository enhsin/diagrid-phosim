;;
;; @package phosim
;; @file validation_3B.pro
;; @brief validation task 3B
;;
;; @brief Created by:
;; @author Mark Hannel (Purdue)
;;
;; @brief Modified by:
;; @author John R. Peterson (Purdue)
;;
;; @warning This code is not fully validated
;; and not ready for full release.  Please
;; treat results with caution.
;;

PRO BNLmeasurement
  ix0=260
  iy0=999
  w=20
  bw=10
  idxmap=intarr(542,2022)
  idxmap[ix0-w-bw:ix0+w+bw,iy0+w:iy0+w+bw]=1
  idxmap[ix0-w-bw:ix0+w+bw,iy0-w-bw:iy0-w]=1
  idxmap[ix0-w-bw:ix0-w,iy0-w-bw:iy0+w+bw]=1
  idxmap[ix0+w:ix0+w+bw,iy0-w-bw:iy0+w+bw]=1
  backidx=where(idxmap EQ 1)
  
  fld=['0V/','5V/','10V/','30V/','50V/','70V/']
  voltage=[0,5,10,30,50,70]
  nv=n_elements(voltage)
  sigma=dblarr(nv)
  sigma_err=dblarr(nv)
  
  for ii=0, nv-1 do begin
    spawn,'ls '+fld[ii]+'*.fits.gz > temp.txt'
    readcol,'temp.txt',files,format='(A)'
    spawn, 'rm temp.txt'
  
    sigma1=dblarr(n_elements(files))
    for jj=0, n_elements(files)-1 do begin
      data0=mrdfits(files[jj],6,header,/silent)
      bzero=FXPAR(header,'BZERO')
      data0=data0+bzero
      data1=data0[ix0-w:ix0+w,iy0-w:iy0+w]
      background=mean(data0[backidx])
      measurepsf,data1-background,rms,e1,e2,medx,medy,flux
      yfit = GAUSS2DFIT(data1, A)
      grms=sqrt(A[2]^2+A[3]^2)/sqrt(2)
      if abs(grms/rms-1) gt 0.05 or FINITE(rms,/NAN) then begin
        print, ii, jj, rms, flux, medx, medy, e1, e2, grms/rms
      endif
      sigma1[jj]=grms*2.0*sqrt(2.0*alog(2.0))*10  ; pixsize=10 microns
    endfor
    sigma[ii]=mean(sigma1)
    sigma_err[ii]=stddev(sigma1)
  endfor
  
  OPENW, 10, 'charge_diffusion_BNL.txt'
  for ii=0, nv-1 do begin
    print, voltage[ii], sigma[ii], sigma_err[ii]
    PRINTF, 10, FORMAT='(F,F,F)', voltage[ii], sigma[ii], sigma_err[ii]
  endfor
  CLOSE, 10

END


pro validation_3B,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 3B'
  !p.multi=0

  !EXCEPT=0

  readcol,'charge_diffusion_BNL.txt',dvoltage,dsigma,dsigma_err,format='(F,F,F)'
  dsigma_err=sqrt(dsigma_err^2+5.0^2)  ; 2 microns sys
  dsigma_sub=sqrt(dsigma^2-14.6^2)


  ix0=2000
  iy0=2036
  w=100
  cc=[0,250,40,229,100]
  ll=[0,1,2,3,4]

  for jj=1, 5 do begin 
    spawn,'/bin/ls lsst_e_31'+strcompress(string(jj),/remove)+'* > temp.txt'
    readcol,'temp.txt',files,format='(A)'
    spawn, 'rm temp.txt'
    voltage=float(strmid(files,strpos(files(0),'3'),5))-(31000+jj*100)
    nv=n_elements(voltage)
    sigma=dblarr(nv)
    
    for ii=0, nv-1 do begin
      data=mrdfits(files[ii],0,header,/silent)
      measurepsf,data[ix0-w:ix0+w,iy0-w:iy0+w],rms,e1,e2,medx,medy,flux
      ;yfit = GAUSS2DFIT(data[ix0-w:ix0+w,iy0-w:iy0+w], A)
      ;rms=sqrt(A[2]^2+A[3]^2)/sqrt(2)
      sigma[ii]=rms*2.0*sqrt(2.0*alog(2.0))     ;pixsize=1 micron
    endfor
    voltage=voltage[0:nv-2]
    sigma=sqrt(sigma[0:nv-2]^2-sigma[nv-1]^2)
  
    if jj eq 1 then begin 
      PLOT, voltage, sigma, LINESTYLE=ll[jj-1], color=cc[jj-1], xtitle='Over Depletion Bias (V)', ytitle='Charge Diffusion (microns FWHM)', xr=[-1,75], yr=[0,31], /xstyle,/ystyle
      value(nnn,0)=sigma(0)
    endif else begin
      OPLOT,  voltage, sigma, LINESTYLE=ll[jj-1], color=cc[jj-1] 
    endelse

  endfor
  OPLOTERR, dvoltage, dsigma, dsigma_err
  OPLOT, dvoltage, dsigma_sub, PSYM=1
  legend,linestyle=ll,['PhoSim r band', '       u band', '       y band','       163.15K', '       183.15K'],box=0,charsize=1.0,/top,/right,color=cc
  legend,psym=[7,1],["BNL lab measurement","BNL sub. 14.6 microns"],box=0,charsize=1.0,/top,/left,color=[0,0]

  ss='Electron Diffusion'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 3B; '+vers
  xyouts,0.7,0.98,ss,/normal

  tolerance_low(nnn,0)=dsigma(0)-dsigma_err(0)
  tolerance_high(nnn,0)=dsigma(0)+dsigma_err(0)

  name(nnn,0)='Diffusion size at 0V'
  unit(nnn,0)=' microns'
  comparison(nnn,0)='BNL Lab Measurements'
  task(nnn,0)='3B Charge Diffusion'


END
