from __future__ import with_statement
import os, sys
import datetime
import subprocess


myFile = sys.argv[1]


logdir = open(myFile).readlines()

for log in logdir:
    log = log.strip()
    logs, filter = log.split('-')
    os.chdir('%s/logs' %(log))
    print 'Working on:', log 
    cmd = 'grep \'This[[:space:]]job[[:space:]]is[[:space:]]starting\' *.out | awk \'{print $11,$7, $8, $9}\' > start.txt'
    subprocess.check_call(cmd, shell=True)

    cmd = 'grep \'PBS[[:space:]]job[[:space:]]finished[[:space:]]at\' *.out | awk \'{print $10, $6, $7, $8}\' > end.txt'
    subprocess.check_call(cmd, shell=True)

    monthMap = {"Jan":"1", "Feb":"2", "Mar":"3", "Apr":"4", "May":"5", "Jun":"6", "Jul":"7", "Aug":"8", "Sep":"9", "Oct":"10", "Nov":"11", "Dec":"12"}

    startTime = open('start.txt').readlines()
    endTime = open('end.txt').readlines()

    sum = 0
    for time in range(0, len(startTime)-1):
        starttime = startTime[time].strip()
        syr, smo, sday, stime = starttime.split(' ')
        sm = monthMap[smo]
        shh, smm, sss = stime.split(':')
    
        endtime = endTime[time].strip()
        eyr, emo, eday, etime = endtime.split(' ')
        em = monthMap[emo]
        ehh, emm, ess = etime.split(':')

        start = datetime.datetime(int(syr), int(sm), int(sday), int(shh), int(smm), int(sss))
        #print 'Start:', start
        end = datetime.datetime(int(eyr), int(em), int(eday), int(ehh), int(emm), int(ess))
        #print 'End  :', end
        delta = end-start
        #print 'Delta:', str(delta)

        #print delta

        hh, mm, ss = str(delta).split(':')
        s1 = int(hh)*3600
        s2 = int(mm)*60
        tot = s1 + s2 + int(ss)
        sum += tot

    average = sum/(len(startTime))
    hours = average/3600.
    #print 'Average Job Time (sec):', average
    print 'Average Job Time (hrs) for %s: %.2f' %(log, hours)
    outfile = 'out.txt'
    os.chdir('../..')
    with file(outfile, 'a') as parFile:
        parFile.write('%s %.2f \n' %(filter, hours))

fuList = []
fgList = []
frList = []
fiList = []
fzList = []
fyList = []

outf = open(outfile).readlines()
for line in outf:
    line = line.strip()
    filt, hrs = line.split(' ')
    
    if filt == 'fu':
        fuList.append(hrs)
    elif filt == 'fg':
        fgList.append(hrs)
    elif filt == 'fr':
        frList.append(hrs)
    elif filt == 'fi':
        fiList.append(hrs)
    elif filt == 'fz':
        fzList.append(hrs)
    else:
        fyList.append(hrs)

#print fzList
#print fiList
filtList = [fuList, fgList, frList, fiList, fzList, fyList]
fList = ['fu', 'fg', 'fr', 'fi', 'fz', 'fy']
x = 0
sum = 0.0
ave = 0.0


## usum = 0
## gsum = 0
## rsum = 0
## isum = 0
## zsum = 0
## ysum = 0

## for fu in fuList:
##     usum += fuList[fu]
## for fg in fgList:
##     gsum += fgList[fg]
## for fr in frList:
##     rsum += frList[fr]
## for fi in fiList:
##     isum += fiList[fi]
## for fz in fzList:
##     zsum += fzList[fz]
## for fy in fyList:
##     ysum += fyList[fy]

## aveu = 0
## aveg = 0
## aver = 0
## avei = 0
## avez = 0
## avey = 0

## aveu = usum/len(fuList)
## print aveu
## aveg = gsum/len(fgList)
## print aveg
## aver = rsum/len(frList)
## print aver
## avei = isum/len(fiList)
## print avei
## avez = zsum/len(fzList)
## print avez
## avey = ysum/len(fyList)
## print avey

## for l in range(0, len(filtList)-1):
##     flist = filtList[l]
##     if len(flist) != 0:
##         for f in range(0, len(flist)-1):
##             fl = flist[f].strip()
##             intf = int(fl)
##             print sum
##             sum += flist[f]
##         ave = sum/len(list)
##         print 'Average %s job time: %.3f' %(fList[x], ave) 
##         x += 1
            
