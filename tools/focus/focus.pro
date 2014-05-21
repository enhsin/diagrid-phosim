focus=fltarr(6,40)
qe=fltarr(6,40)
xfocus=fltarr(6,40)

for filter=0,5 do begin
for ii=0L,39 do begin
defocus=-0.020+0.001*ii
openw,1,'object_test1'
if filter eq 0 then printf,1,'filter 0'
if filter eq 1 then printf,1,'filter 1'
if filter eq 2 then printf,1,'filter 2'
if filter eq 3 then printf,1,'filter 3'
if filter eq 4 then printf,1,'filter 4'
if filter eq 5 then printf,1,'filter 5'
printf,1,'chipid R22_S11'

sstring='body 11 5 '+string(defocus)
printf,1,sstring
;printf,1,'source ../ancillary/standardatm'
printf,1,'natmospherefile 0'
printf,1,'object 0.0 0.0 18 ../ancillary/sky/sed_flat.txt 0.0 star'
printf,1,'outputfilename test'
printf,1,'lsst'
close,1

spawn,'lsst < object_test1'

data=mrdfits('test.fits',0)
ssize=100.0
meanx=2048.0
meany=2048.0
tot=total(data)
minchi=1e30
for sigma=0.0,2.0,0.01 do begin
chi=0.0
for i=2000L,2100 do begin
for j=2000L,2100 do begin
        gg=tot*exp(-((float(i)-meanx)^2+(float(j)-meany)^2)/2.0/sigma/sigma)/(2*!Pi*sigma*sigma)
    chi=chi+(gg-data(i,j))^2/(data(i,j)>1) 
endfor
endfor

if chi lt minchi then begin
    minchi=chi
    ssize=sigma
endif
;print,sigma,chi,ssize
endfor


xfocus(filter,ii)=defocus
focus(filter,ii)=ssize
qe(filter,ii)=tot

set_plot,'PS'
device,filename='focus.ps'
device,/color
loadct,39
plot,xfocus(0,*),focus(0,*),color=0,yr=[0,1],xtitle='Defocus Position (mm)',ytitle='1 Sigma PSF Size (pixels)'
oplot,xfocus(1,*),focus(1,*),color=60
oplot,xfocus(2,*),focus(2,*),color=100
oplot,xfocus(3,*),focus(3,*),color=140
oplot,xfocus(4,*),focus(4,*),color=180
oplot,xfocus(5,*),focus(5,*),color=220
device,/close

endfor
endfor


END
