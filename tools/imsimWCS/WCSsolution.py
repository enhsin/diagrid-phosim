import os
import pyfits
import numpy
from numpy import *
from numpy.linalg import *
import pylab as N
from math import pi
from math import *
import sys

# Jim Pizagno 7 May 2009
## parameters that can be adjusted:
#defsexlines = 14
#Dmag_limit = 0.5  # DC3 run: stars mag changes by 0.2-0.4 mag from airmass=1 -> arimass=2. for conservative=0.75
#radius_limit=47.0  #originally 50.0, 75.0 may be better for DC3_final/Wide*fits.
#radius_min=0.0
#sigmaclip=0.2
#dc3info=open('/astro/net/lsst1/jpizagno/DC3_final/dc3.txt','r')
#file = '/astro/net/lsst1/jpizagno/DC3/Catalogs/StarCatalogFormatTrim.dat'
#file = '/astro/net/lsst1/jpizagno/DC3_final/Catalogs/StarCatalogFormatTrim_4x.dat'
## the names of any files being read in or written out.
##
## run this program by typing:   "%python WCSsolution.py image.fits"

def SigmaClip(Dx,Dy,sigmaclip,Xmeasured,Ymeasured,rainput,decinput,Dmag):
	ibad=[]
	for i in range(0,len(Dx)):
		distfrommeanSigX = fabs(Dx[i] - mean(Dx))/std(Dx)
		distfrommeanSigY = fabs(Dy[i] - mean(Dy))/std(Dy)
		if distfrommeanSigX > sigmaclip or distfrommeanSigY > sigmaclip:
			ibad.append(i)
	for m in range(0,len(ibad)):
		# print 'removing i=',ibad[m]-m,' Dx=',Dx[ibad[m]-m],' Dy=',Dy[ibad[m]-m],' \n'
		Dx.remove(Dx[ibad[m]-m])
		Dy.remove(Dy[ibad[m]-m])
		Xmeasured.remove(Xmeasured[ibad[m]-m])			
		Ymeasured.remove(Ymeasured[ibad[m]-m])
		rainput.remove(rainput[ibad[m]-m])
		decinput.remove(decinput[ibad[m]-m])	
		Dmag.remove(Dmag[ibad[m]-m])			
	return


# Description
# Given a catalog with x y ra dec and a fits image fit a wcs and update the fits header information.
# We assume that the input catalogs have the form
#    x and y with lower left as the origin
#    ra and dec are in degrees.
# Originally fitWCS.v1.py
def fitWCS(Xin,Yin,Xin2,Yin2,XYin,x2y2,RAin,Decin,fitsFile,writeoutfile,UpDateHeader):
        RA = RAin        
        Dec = Decin      
        #Create Z Matrices (1, x, y, x^2, y^2, xy) with x,y the CCD coordinates
        Z = N.ones((len(Xin),6))              #   original:  N.ones((Xin.shape[0],3))
        Z[:,1] = Xin                          #   original:  data[:,0]
        Z[:,2] = Yin                          #   original:  data[:,1]
	Z[:,3] = Xin2
	Z[:,4] = Yin2
	Z[:,5] = XYin
	converge_test=0
	number=0
	Beta1 = N.linalg.lstsq(Z,RA)   # FITTING MACHINERY DONE HERE, lstsq() does the (re)fitting
        Beta2 = N.linalg.lstsq(Z,Dec)
	Beta =  N.vstack((Beta1[0],Beta2[0]))	
# apply matrix to get positions  in x and y and RA Dec
       	SkyPostions =  N.dot(Beta,Z.transpose())
	
	if writeoutfile==1:
		outfile = fitsFile+'RADEC.HigherOrder.Residuals.dat'
		try:
			f=open(outfile,'w')
		except:
			print 'ERROR in WCS.  couldnt open file=',file
		f.write("#RA_input(deg)  Dec_input(deg)  RA_Solution-RA_input(arcsec) Dec_Solution-Dec_input(arcsec)  X(pix)  Y(pix) \n") 
		for i in range(0,len(RA)-1):
			f.write(str(RA[i])+" "+str(Dec[i])+" "+str((SkyPostions[0][i]-RA[i])*3600.0)+" "+str((SkyPostions[1][i]-Dec[i])*3600.0)+" "+str(Xin[i])+" "+str(Yin[i])+"  \n")
		f.close()
        RADEC = N.vstack((RA,Dec))
        CCDPositions = N.dot(N.linalg.pinv(Beta), RADEC.transpose()[0])
	# create WCS (origin and rotation matrix)
        matrix = N.vstack((Beta[:,1],Beta[:,2]))
        pos = RADEC.transpose() - Beta[:,0]
	# Open fits Read in fits file
	if UpDateHeader == 1:
		try:
           		hdulist = pyfits.open(fitsFile,mode="update")
        	except:
           		print "Cant read fits file",fitsFile
		#update fits header
        	header = hdulist[0].header
		header.update('EQUINOX',2000.0,'J2000 Equinox')  #change to float
        	header.update('RADESYS','FK5','')
        	header.update('CTYPE1','RA---TAN-SIP','')
		header.update('CTYPE2','DEC--TAN-SIP','')
		header.update('CRVAL1',Beta[0,0],'')
		header.update('CRVAL2',Beta[1,0],'')
		header.update('CD1_1',Beta[0,1],'')
        	header.update('CD1_2',Beta[0,2],'')
        	header.update('CD2_1',Beta[1,1],'')
        	header.update('CD2_2',Beta[1,2],'')		
        	header.update('CRPIX1',0.0,'')
		header.update('CRPIX2',0.0,'')
		cd11=Beta[0,1]
		cd12=Beta[0,2]
		cd21=Beta[1,1]
		cd22=Beta[1,2]
		A20=(Beta[0,3] - (cd12/cd22)*Beta[1,3])/(cd11 - (cd12/cd22)*cd21)
		B20=(Beta[1,3]-cd21*A20)/cd22
		A02=(Beta[0,4] - (cd12/cd22)*Beta[1,4])/(cd11 - (cd12/cd22)*cd21)
		B02=(Beta[1,4]-cd21*A02)/cd22
		A11=(Beta[0,5] - (cd12/cd22)*Beta[1,5])/(cd11 - (cd12/cd22)*cd21)
		B11=(Beta[1,5]-cd21*A11)/cd22
		header.update('A_ORDER',2,'polynomial order, axis 1')
		header.update('A_0_2',A02,'distortion coefficient')
		header.update('A_1_1',A11,'distortion coefficient')
		header.update('A_2_0',A20,'distortion coefficient')
		header.update('B_ORDER',2,'polynomial order, axis 2')
		header.update('B_0_2',B02,'distortion coefficient')
		header.update('B_1_1',B11,'distortion coefficient')
		header.update('B_2_0',B20,'distortion coefficient')
		header.update('CDELT1',0.2/3600.0,'place scale in degrees per pixel')  #LSST=0.2arcsec/pix*/3600(arcsec/deg)
		header.update('CDELT2',0.2/3600.0,'plate scale in degrees per pixel')
        	# Calculate Inverse Transform.
		ZP = N.ones((len(Xin),6))
		factor=(cd11*cd22-cd12*cd21)
		cdP11=cd22/factor
		cdP12=-cd12/factor
		cdP21=-cd21/factor
		cdP22=cd11/factor
		for i in range(0,len(RA)):
			U = cdP11*(RA[i]-Beta[0,0]) + cdP12*(Dec[i]-Beta[1,0])
			V = cdP21*(RA[i]-Beta[0,0]) + cdP22*(Dec[i]-Beta[1,0])
			ZP[i,1] = U                       #   original:  data[:,0]
        		ZP[i,2] = V                         #   original:  data[:,1]
			ZP[i,3] = U*U
			ZP[i,4] = V*V
			ZP[i,5] = U*V
			
        	# set up linear equations and solve for rotation matrix
        	Beta1P = N.linalg.lstsq(ZP,Xin)   
        	Beta2P = N.linalg.lstsq(ZP,Yin)
        	BetaP = N.vstack((Beta1P[0],Beta2P[0]))
		AP20=BetaP[0,3]
		BP20=BetaP[1,3]
		AP02=BetaP[0,4]
		BP02=BetaP[1,4]
		AP11=BetaP[0,5]
		BP11=BetaP[1,5]
		header.update('AP_ORDER',2,'polynomial order, axis 1')
		header.update('AP_0_2',AP02,'distortion coefficient')
		header.update('AP_1_1',AP11,'distortion coefficient')
		header.update('AP_2_0',AP20,'distortion coefficient')
		header.update('BP_ORDER',2,'polynomial order, axis 2')
		header.update('BP_0_2',BP02,'distortion coefficient')
		header.update('BP_1_1',BP11,'distortion coefficient')
		header.update('BP_2_0',BP20,'distortion coefficient')
		hdulist.close()
        return

# get FITS image header information for each chip.
# get bounds.
# add 1arcmin buffer to bounds when searching catalog.
def SimpleFit(image,ra_star,dec_star,mag_star,sed_star,x_pred,y_pred):
    	try:
        	hdulist = pyfits.open(image)
   	except:
        	print "Cant read fits file"
	prihdr=hdulist[0].header
	naxis1 = hdulist[0].header['NAXIS1']
	naxis2 = hdulist[0].header['NAXIS2']
	crpix1 = float(hdulist[0].header['CRPIX1'])
	crval1 = float(hdulist[0].header['CRVAL1'])  # ra
	crpix2 = float(hdulist[0].header['CRPIX2'])
	crval2 = float(hdulist[0].header['CRVAL2'])  # dec
	cd1_1 = float(hdulist[0].header['CD1_1'])
	cd1_2 = float(hdulist[0].header['CD1_2'])
	cd2_1 = float(hdulist[0].header['CD2_1'])
	cd2_2 = float(hdulist[0].header['CD2_2'])
	hdulist.close()
	print 'header info = ',cd1_1,cd1_2,cd2_1,cd2_2
	tmpra=zeros(4)
	tmpdec=zeros(4)
	tmpra[0] = cd1_1*(1.0-crpix1) + cd1_2*(1.0-crpix2) + crval1  # ra
	tmpdec[0] = cd2_1*(1.0-crpix1) + cd2_2*(1.0-crpix2) + crval2  # dec
	tmpra[1] = cd1_1*(1.0-crpix1) + cd1_2*(naxis2-crpix2) + crval1
	tmpdec[1] = cd2_1*(1.0-crpix1) + cd2_2*(naxis2-crpix2) + crval2
	tmpra[2] = cd1_1*(naxis1-crpix1) + cd1_2*(naxis2-crpix2) + crval1
	tmpdec[2] = cd2_1*(naxis1-crpix1) + cd2_2*(naxis2-crpix2) + crval2
	tmpra[3] = cd1_1*(naxis1-crpix1) + cd1_2*(1.0-crpix2) + crval1
	tmpdec[3] = cd2_1*(naxis1-crpix1) + cd2_2*(1.0-crpix2) + crval2
	ra_low=min(tmpra) - 1.0/60.0
	ra_high=max(tmpra) + 1.0/60.0
	dec_low=min(tmpdec) - 1.0/60.0
	dec_high=max(tmpdec) + 1.0/60.0  #add an arcmin to all
# get sources of stars within bounds. 
# calculate PREDICTED x/y from input ra/dec and WCS solution.
# output or save Mag, ID, ra, dec, X_predicted, Y_predicted.
# output only bright objects??
	mag_input_limit=25.0  # only stars brighter than this mag
	GAL_mag_input_limit=100.0   # get ALL Galaxies, output CCDChip Catalog, output regions, Do NOT use in fit.
	print 'Calculating predicted XY'
#	file=image+'XY.pred.dat'
#	f=open(file,'w')
	file2=image+'RADEC.InputCatalog.Stars.reg'
	f2=open(file2,'w')
	f2.write('# Regions predicted by INPUT STAR Catalog:  /astro/net/lsst1/jpizagno/DC3/Catalogs/StarCatalogFormatTrim.dat \n')
	f2.write('global color=green font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source \n')
	f2.write('fk5 \n')	
	for i in range(0,len(ra_star)):
		if ra_star[i] < ra_high and ra_star[i] > ra_low:
			if dec_star[i] < dec_high and dec_star[i] > dec_low:
	                        if mag_star[i] < mag_input_limit:  #only get bright stars, maybe need more
	                                x_pred[i] = (cd2_2*(ra_star[i]-crval1)-cd1_2*(dec_star[i]-crval2))/(cd2_2*cd1_1-cd1_2*cd2_1) + crpix1
	                                y_pred[i] = ((dec_star[i]-crval2) - cd2_1*(x_pred[i]-crpix1))/(cd2_2)+crpix2
					f2.write("point("+str(ra_star[i])+","+str(dec_star[i])+") # point=circle \n")
	f2.close()
# run SEXtractor with os module.
# set parameter file to only look at objects with S/N > 20, or bright.
	numsexlines=0
	print 'Running SExtracton on ',image
	os.system("rm test.cat")
	os.system("sex "+image)
	os.system("wc -l test.cat > delme")
	f=open("delme","r")
	line=f.readline()
	a=line.split()
	f.close()
	os.system("rm delme")
	defsexlines=14
	numsexlines=int(a[0])-defsexlines  #=number of lines to Skip.   may have to change 14.
	print 'read ',numsexlines,' lines in test.cat for SExtractor result run on ',image
# from SEX run get mag, x_IMAGE, y_image, ra_image,dec_image.
	print 'Getting SExtractor results image = ',image
	file='test.cat'
	f=open(file,'r')
	mag_sex=zeros(numsexlines)
	ra_sex=zeros(numsexlines)
	dec_sex=zeros(numsexlines)
	x_sex=zeros(numsexlines)
	y_sex=zeros(numsexlines)
	xwin_sex=zeros(numsexlines)
	ywin_sex=zeros(numsexlines)
	i=-1
	for line in f.readlines():
		a=line.split()
		if a[0] != '#':
			if a[15] =='0':  #only use objects without flag.
				i=i+1
				mag_sex[i] = float(a[2])
				x_sex[i] = float(a[10])
				y_sex[i] = float(a[11])
				ra_sex[i] = float(a[12])
				dec_sex[i] = float(a[13])
				xwin_sex[i] = float(a[8])
				ywin_sex[i] = float(a[9])
	f.close()

# Calculate Simple solution.
# search for all objects within 300 pixels  within 0.3 magnitudes., haveing SNR>20 (specified default.sex)
# calculate simple linear solution.  Dx Dy.
	Dx=[]
	Dy=[]
	Xmeasured = []
	print 'size of Xmeasured before loop = ',len(Xmeasured),' \n'
	Ymeasured = []
	rainput = []
	decinput = []
	DatOut=[]
	Dmag=[]
	Dmag_limit = 0.5  # DC3 run: stars mag changes by 0.2-0.4 mag from airmass=1 -> arimass=2. for conservative=0.75
	radius_limit=47.0  #originally 50.0, 75.0 may be better for DC3_final/Wide*fits.
	radius_min=0.0
	print 'Making simple comparison between positions and mags for image = ',image,' \n'
	for i in range(0,len(ra_star)):
		catches=float(0.0)
		if x_pred[i] > 0:  #only search stars within this chip
			for j in range(0,numsexlines):
        	                Dmagtmp=fabs(mag_sex[j] - mag_star[i])  #for testing only
				if Dmagtmp < Dmag_limit:
					r=sqrt(pow(x_pred[i]-xwin_sex[j],2)+pow(y_pred[i]-ywin_sex[j],2))
					if r < radius_limit and r > radius_min:
						catches=catches+1.0
						if catches < 2.0: #print out stars with ONE hit.
							# avoids repeats. 
							test=0
							for m in range(0,len(Xmeasured)):
								if xwin_sex[j] == Xmeasured[m]:
									test=1  #found a repeat.
							if test==0:
	                                                	Dx.append(x_pred[i]-xwin_sex[j])
        	                                       		Dy.append(y_pred[i]-ywin_sex[j])
                	                                	Xmeasured.append(xwin_sex[j])
                        	                        	Ymeasured.append(ywin_sex[j])
                                	                	rainput.append(ra_star[i])
                                        	        	decinput.append(dec_star[i])
								Dmag.append(float(mag_star[i] - mag_sex[j]))  # M_AB-measured
	print 'size of Xmeasured AFTER loop = ',len(Xmeasured),' \n'
	print 'mean(Dx) = ',mean(Dx),' std(Dx) = ',std(Dx),' \n'
	print 'mean(Dmag) = ',mean(Dmag),' std(Dmag) = ',std(Dmag),' \n'
	
# prune outliers.  if len(Dx)>20.  
	sigmaclip=0.2  # clip everything outside this sigma limit.  
	if len(Dx) > 20:
		print mean(Dx),std(Dx),median(Dx),mean(Dy),std(Dy),median(Dy)	
		print 'calling SigmaClip()'
		SigmaClip(Dx,Dy,sigmaclip,Xmeasured,Ymeasured,rainput,decinput,Dmag)
	print 'size of Xmeasured after PRUNING = ',len(Xmeasured),' \n'
	print 'mean(Dx) = ',mean(Dx),' std(Dx) = ',std(Dx),' \n'
	print 'mean(Dmag) = ',mean(Dmag),' std(Dmag) = ',std(Dmag),' \n'
	file2=image+'SEx.StarsUsedinFit.reg'
	f2=open(file2,'w')
	f2.write('# Regions predicted by INPUT \n')
	f2.write('global color=blue font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source \n')
	f2.write("fk5 \n")
	for i in range(0,len(rainput)):	
		f2.write("point("+str(rainput[i])+","+str(decinput[i])+") # point=circle \n")
	f2.close()
		

#  Use ra/dec from input catalog, and x/y from SExtractor.
#  set up x^2 y^2 x*y arrays.
	x2=[] #x * x
	y2=[] #y * y
	xy=[] #x * y
	x2y2=[] #x^2 + y^2
	for i in range(0,len(Xmeasured)):
		x2.append(Xmeasured[i]*Xmeasured[i])
		y2.append(Ymeasured[i]*Ymeasured[i])
		xy.append(Xmeasured[i]*Ymeasured[i])
		x2y2.append(Xmeasured[i]*Xmeasured[i] + Ymeasured[i]*Ymeasured[i])
	if len(Xmeasured) < 20:
		print '***************************************************** WARNING     WARNING*********'
		print 'len(Xmeasured)<20 so NOT calling fitWCS'
		exit()
	else:
		print 'calling fitWCS on image ',image,' \n'
		fitWCS(Xmeasured,Ymeasured,x2,y2,xy,x2y2,rainput,decinput,image,1,1)
	return

def main():
    print 'This software requires Source Extractor'
    try:
        image = sys.argv[1]
    except:
        print "USAGE: fitWCS image"
        exit()

    #get sources of stars within bounds. get input ra/dec.
    print 'getting all stars at beginning. \n'
    # Figure out which file to get:
    dc3info=open('/astro/net/lsst1/jpizagno/DC3_final/dc3.txt','r')
    test=0
    for i in range(0,1117):
	    line=dc3info.readline()
	    a=line.split()
 	    tmp=a[0]+".fits"
	    if tmp==image:  #found a hit
		    filetester=a[1]
   		    if a[1] == '0':  # then get regular file
			    file = '/astro/net/lsst1/jpizagno/DC3/Catalogs/StarCatalogFormatTrim.dat'
			    test=1
		    if a[1] =='1':   # get file with 4x density of stars.
			    file = '/astro/net/lsst1/jpizagno/DC3_final/Catalogs/StarCatalogFormatTrim_4x.dat'
			    test=1
    dc3info.close()
    if test==0:
	    print '*******couldnt find matching image name in dc3.txt why? **********WARNING !!!!! \n'
    try:
    	f=open(file,'r')
    except:
	print 'ERROR.  couldnt open file=',file
    ra_star=[] 
    dec_star=[] 
    mag_star=[] 
    sed_star=[] 
    x_pred=[] 
    y_pred=[] 
    for line in f.readlines():
            a=line.split()
            ra_star.append(float(a[1]))
            dec_star.append(float(a[2]))
            mag_star.append(float(a[3]))
	    sed_star.append(a[4])
            x_pred.append(float(-1.E10))
            y_pred.append(float(-1.E10))
    f.close()


# Call routine for matching.
    SimpleFit(image,ra_star,dec_star,mag_star,sed_star,x_pred,y_pred)

if __name__ == "__main__":
    main()

