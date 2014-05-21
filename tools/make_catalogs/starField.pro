pro starField, spacing
;spacing giving in arcseconds

;parameters
Unrefracted_RA = 0 ;pointingra
Unrefracted_Dec = 0 ;pointingdec
Unrefracted_Azimuth = 0 ;azimuth * 3.14159 / 180
Unrefracted_Altitude = 1.57080 ;(90 - zenith) * 3.14159 / 180
Slalib_date = '1994/7/19/0.298822999997' ;month
SIM_SEED = 17219418 ;seed
Opsim_rotskypos = 0 ;rotationalangle
Opsim_obshistid = 42 ;number ; only needed for header
Opsim_rottelpos = 0 ;spiderangle
Opsim_moondec = -0.348541 ;moon_dec
Opsim_rawseeing = 0.958908 ;constrainedseeing * 2.35482
Opsim_moonra = 4.29155 ;moon_ra
Opsim_expmjd = 49552.298823 ;number ; only needed for header
Opsim_moonalt = 0.199561 ;moon_alt
Opsim_sunalt = -0.993843 ;(3.14159 / 2.0) - sunzen
Opsim_filter = 2 ;filter
Opsim_dist2moon = 2.00046 ;moon_dist

;star spacing calculations
dis = spacing / 206264.
maxdis = FIX(0.035 / dis)
print, dis
print, maxdis
mag = 22 ;magnitude

;printing catalog
openw,1,'star_grid'
printf,1,'Unrefracted_RA '+string(Unrefracted_RA)
printf,1,'Unrefracted_Dec '+string(Unrefracted_Dec)
printf,1,'Unrefracted_Azimuth '+string(Unrefracted_Azimuth)
printf,1,'Unrefracted_Altitude '+string(Unrefracted_Altitude)
printf,1,'Slalib_date '+string(Slalib_date)
printf,1,'SIM_SEED '+string(SIM_SEED)
printf,1,'Opsim_rotskypos '+string(Opsim_rotskypos)
printf,1,'Opsim_obshistid '+string(Opsim_obshistid)
printf,1,'Opsim_rottelpos '+string(Opsim_rottelpos)
printf,1,'Opsim_moondec '+string(Opsim_moondec)
printf,1,'Opsim_rawseeing '+string(Opsim_rawseeing)
printf,1,'Opsim_moonra '+string(Opsim_moonra)
printf,1,'Opsim_expmjd '+string(Opsim_expmjd)
printf,1,'Opsim_moonalt '+string(Opsim_moonalt)
printf,1,'Opsim_sunalt '+string(Opsim_sunalt)
printf,1,'Opsim_filter '+string(Opsim_filter)
printf,1,'Opsim_dist2moon '+string(Opsim_dist2moon)

printf,1,'object 0 -0.0 0.0 '+string(mag)+' ../ancillary/sky/sed_flat.txt 0 0 0 0 0 0 star none none'

for i = 1, maxdis do begin
   printf,1,'object 0 '+string(-0.0 + i * dis)+' 0.0 '+string(mag)+' ../ancillary/sky/sed_flat.txt 0 0 0 0 0 0 star none none'
   printf,1,'object 0 '+string(-0.0 - i * dis)+' 0.0 '+string(mag)+' ../ancillary/sky/sed_flat.txt 0 0 0 0 0 0 star none none'
endfor

for j = 1, maxdis do begin
   printf,1,'object 0 -0.0 '+string(0.0 + j * dis)+' '+string(mag)+' ../ancillary/sky/sed_flat.txt 0 0 0 0 0 0 star none none'
   printf,1,'object 0 -0.0 '+string(0.0 - j * dis)+' '+string(mag)+' ../ancillary/sky/sed_flat.txt 0 0 0 0 0 0 star none none'
endfor

for i = 1, maxdis do begin
   for j = 1, maxdis do begin
      printf,1,'object 0 '+string(-0.0 + i * dis)+' '+string(0.0 + j * dis)+' '+string(mag)+' ../ancillary/sky/sed_flat.txt 0 0 0 0 0 0 star none none'
      printf,1,'object 0 '+string(-0.0 - i * dis)+' '+string(0.0 - j * dis)+' '+string(mag)+' ../ancillary/sky/sed_flat.txt 0 0 0 0 0 0 star none none'
      printf,1,'object 0 '+string(-0.0 + i * dis)+' '+string(0.0 - j * dis)+' '+string(mag)+' ../ancillary/sky/sed_flat.txt 0 0 0 0 0 0 star none none'
      printf,1,'object 0 '+string(-0.0 - i * dis)+' '+string(0.0 + j * dis)+' '+string(mag)+' ../ancillary/sky/sed_flat.txt 0 0 0 0 0 0 star none none'
   endfor
endfor
close,1


end
