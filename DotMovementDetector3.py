from ij import IJ
from ij.gui import Roi
from pta import PTA
from pta.gui import ShowPdata
from pta.track import DetectParticle
from pta.data import PtaParam
from java.util import ArrayList
from ij.plugin.frame import ThresholdAdjuster
from ij.process import ImageProcessor
import math
from ij.plugin.frame import RoiManager
from pta.data import FPoint
from ij.gui import ShapeRoi
from java.awt import Color


class FPext(object):
	hasMoved = False
	isMovingNext = False
	checked = False
	hasNext = False
	hasPrevious = False
	def __init__(self, fp):
		self.fp = fp

	def getFP(self):
		return self.fp

# set data path
#datafilepath = "/Users/miura/Dropbox/20130805_Osaka/data/transferrin-movement/tconc1_3.tif"
#datafilepath = "/Users/miura/Dropbox/20130805_Osaka/data/transferrin-movement/tconc10_16.tif"
#datafilepath = '/Users/miura/Dropbox/meetings/20130805_Osaka/data/transferrin-movement/tconc10_16.tif'
datafilepath = '/Users/miura/Dropbox/people/Tina/migrating_f11_13.tif'
#datafilepath = '/Users/miura/Dropbox/people/Tina/migration_f8_f15.tif'
datafilepath = '/Users/miura/Dropbox/people/Tina/migration_small_f1_f25.tif'
datafilepath = '/Users/miura/Dropbox/people/Tina/cell7_from5withROCK_1ulvirus_ON_cropvirus_part.tif'

##### detection parameters ######
scanROIsize = 10
dotMinIntensity = 80.0
dotMinSize = 3
###
criticalDistance = 0.5
####################

# load data as an ImagePlus object
imp = IJ.openImage(datafilepath)

# default ROI = full frame
scanAreaRoi = Roi(0,0,imp.getWidth(),imp.getHeight())

# set intensity threshold to the iamge
ip = imp.getProcessor()
ht = ip.getAutoThreshold()
#ip.setThreshold(0, ht, ImageProcessor.RED_LUT)
ip.setThreshold(ht, 255, ImageProcessor.RED_LUT)

# **** Set PTA parameters ****

PTA.setDetectionState(True)
#PTA.setDebugMode()
# set no GUI mode
PTA.setNoGUI(True)

# set Detection Parameters.
# PtaParam.Builder(int roiSizex,int roiSizey,boolean do2dGaussfit)
ptap = PtaParam.Builder(scanROIsize,scanROIsize, False).build()

### more ways to customize paramters
#ptap.setDo2dGaussfit(False)
#ptap.setDo2dGaussfit(True) 
ptap.setMinIntensity(dotMinIntensity) #default 100
#ptap.setSearchPointIncrement(1)
ptap.setMinSize( dotMinSize ) 

# **** Tracking ****

# instantiate a DetectParticle object
# DetectParticle(PtaParam ptap, ij.ImagePlus imp, ij.gui.Roi scanRoi, boolean nogui)
dp = DetectParticle(ptap,imp,scanAreaRoi, True)

# set range of frames to analyze
dp.setStackRange(1, imp.getStackSize())

# start threads. 
dp.start()
# wait till the processing finishes
dp.join()
# stop the thread. 
dp.interrupt()

PTA.setDetectionState(False)

# **** output results ****

# show detection results
imp.show()

# Particle lists
ll = dp.getalldplist()

### get ROI_Manager 
rm = RoiManager.getInstance()
if rm is None:
    rm = RoiManager()
else:
    rm.runCommand('reset')
    		
	
# track lists
tracks = dp.getLinkedPointList()
#print "=== tracks ==="
#for t in tracks:
#   if len(t) > 2: 
#      print t


plist = []
for t in tracks:
    if t[0].getFrame() == 1:
        for p in t:
            plist.append(p)
            print p.getFrame(), p.getSx(), p.getSy()

for i in range(imp.getStackSize()):
    rr = None
    for p in plist:
        if p.getFrame() == i + 1:
            if rr is None:
                rr = ShapeRoi(p.retRoi())
            else:
                rr = rr.or(ShapeRoi(p.retRoi()))
    if rr is not None:
		rr.setPosition( i + 1 )
		rr.setStrokeColor(Color.GREEN)
		rm.addRoi(rr)	


PTA.setNoGUI(False)
