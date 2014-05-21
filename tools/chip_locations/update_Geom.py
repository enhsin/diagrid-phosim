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



fh_geom = open(sys.argv[1])
name = ''
for l in fh_geom:
    flds = l.rstrip().split()
    if l.strip().startswith("name:"):
        namestr = l.strip().split(":")
        if len(namestr) == 4:
            name = "R%s%s_S%s%s"%(namestr[2][0], namestr[2][2], namestr[3][0],
                    namestr[3][2])
        print l.rstrip()

    elif len(name) > 0 and l.strip().startswith("offset:"):
        xoff, yoff, rot = readOffset(name)
        flds = l.strip().split()
        print "        offset:",float(flds[1])+xoff,float(flds[2])+yoff

    elif len(name) > 0 and l.strip().startswith("orientation:"):
        xoff, yoff, rot = readOffset(name)
        flds = l.strip().split()
        print "        orientation: 0.000000 0.000000 %f"%(rot)
    else:
        print l.rstrip()
