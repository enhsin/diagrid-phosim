showq | grep nms | grep -v Running | awk '{print "canceljob "$1}' | sh
showq | grep nms | grep -v Idle | awk '{print "canceljob "$1}' | sh
showq | grep nms | grep -v Hold | awk '{print "canceljob "$1}' | sh
rocks run host command="rm -rf /state/partition1/nms/8*"
