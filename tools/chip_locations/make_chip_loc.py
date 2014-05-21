import sys
import math
def readOffset(id):
    fh = open("../../data/focal_plane/sta_misalignments/offsets/pars_%s"%(id))
    xoff = None
    yoff = None
    rot = None
    for l in fh:
        if l.startswith("body 11 3"):
            xoff = float(l.rstrip().split()[3])
        if l.startswith("body 11 4"):
            yoff = float(l.rstrip().split()[3])
        if l.startswith("chipangle"):
            rot = float(l.rstrip().split()[1])
    return xoff, yoff, rot
def readCorners(filename):
    fh = open(filename)
    corners = {}
    id = None
    vals = []
    for l in fh:
        if l.startswith("R"):
            if len(vals) == 8:
                positions = [(vals[0], vals[1]), (vals[2], vals[3]), (vals[4],\
                    vals[5]), (vals[6], vals[7])]
                corners[id] = positions
            id = l.rstrip()
            vals = []
        else:
            vals.append(float(l.strip()))
    positions = [(vals[0], vals[1]), (vals[2], vals[3]), (vals[4],\
         vals[5]), (vals[6], vals[7])]
    corners[id] = positions
    return corners



fh_chip = open(sys.argv[1])
pixscale = 0.01
corns = readCorners(sys.argv[2])
for l in fh_chip:
    flds = l.rstrip().split()
    xcent = float(flds[1])/1000.
    ycent = float(flds[2])/1000.
    xpix = int(flds[3])
    ypix = int(flds[4])
    dx = pixscale*xpix/2.
    dy = pixscale*ypix/2.
    xoff,yoff,rot = readOffset(flds[0])
    print flds[0].strip()
    id = flds[0].strip()
    xp = []
    yp = []
    xo = []
    yo = []
    for i in ((-1,-1),(-1,1),(1,1),(1,-1)):
        xp.append(i[0]*dx*math.cos(rot*math.pi/180.) -
                i[1]*dy*math.sin(rot*math.pi/180.)+xcent+xoff)
        yp.append(i[0]*dx*math.sin(rot*math.pi/180.) +
                i[1]*dy*math.cos(rot*math.pi/180.)+ycent+yoff)
        xo.append(i[0]*dx + xcent)
        yo.append(i[1]*dy + ycent)
        
    if not xp[2] == xp[3]:
        print "Rotation: ",(math.atan((yp[2] - yp[3])/(xp[2] - xp[3]))*180./math.pi)%90., rot%90., \
             (math.atan((corns[id][2][1] - corns[id][3][1])/(corns[id][2][0] - corns[id][3][0]))*180./math.pi)%90.

    else:
        print "Rotation: ",0., rot
    print "sx sy nx ny jx jy"
    for pos in zip(xp,yp,xo,yo,corns[id]):
        print " ".join([str(el) for el in pos])
