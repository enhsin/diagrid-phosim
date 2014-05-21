PRO JOIN_FITS, home_dir

;-----------------------------------------------
;
; NAME: JOIN_FITS
;
; PURPOSE: To combine 25 fits files to make one
;          larger fits file of the central focalplane of lsst
;
; INPUTS: The 25 .fits files to be joined
;
; OUTPUTS: .fits file of the combined files
;
; FUNCTIONS USED: NONE
;
; WRITTEN FOR: John Peterson, Purdue University
;
; BY: Alan Meert
;
; DATE: 11 JUNE 2008
;
;-----------------------------------------------

COMPILE_OPT IDL2

data = lonarr(20500,20500)
FOR x = -2, 2 DO BEGIN
    FOR y = -2, 2 DO BEGIN
        
        x_num = num2str(x)
        y_num = num2str(y)
        file_name = home_dir + '/' + 'focalplane_' + x_num + '_' + y_num + '.fits'

        input = mrdfits(file_name)

        x_offset = (x+2)*4101
        y_offset = (y+2)*4101

        FOR countx = 0, 4095 DO BEGIN
            FOR county = 0, 4095 DO BEGIN
                data[countx + x_offset, county + y_offset] = input[countx, county]
            ENDFOR
        ENDFOR
    ENDFOR
ENDFOR

WRITEFITS, 'total.fits', data

END
