# @File(label="Select a image stack") dotimage

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
from ij.plugin import Duplicator


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
#datafilepath = '/Users/miura/Dropbox/people/Tina/migrating_f11_13.tif'
#datafilepath = '/Users/miura/Dropbox/people/Tina/migration_f8_f15.tif'
#datafilepath = '/Users/miura/Dropbox/people/Tina/migration_small_f1_f25.tif'
#datafilepath = '/Users/miura/Dropbox/people/Tina/cell7_from5withROCK_1ulvirus_ON_cropvirus_part.tif'
datafilepath = dotimage.getAbsolutePath()

##### detection parameters ######
scanROIsize = 12
dotMinIntensity = 2.0
dotMinSize = 2
linkageFrame = 2
###
criticalDistance = 0.5
####################

# load data as an ImagePlus object
imp = IJ.openImage(datafilepath)
dupimp = Duplicator().run(imp)

# default ROI = full frame
scanAreaRoi = Roi(0,0,imp.getWidth(),imp.getHeight())

# set intensity threshold to the image
#ip = imp.getStack().getProcessor(imp.getStackSize())
ip = imp.getProcessor()
ht = ip.getAutoThreshold()
print "Threshold: ", ht
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
ptap.setLinkageFrame(linkageFrame) 

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
#imp.show()

# Particle lists
ll = dp.getalldplist()

    		
# track lists
tracks = dp.getLinkedPointList()
#print "=== tracks ==="
#for t in tracks:
#   if len(t) > 2: 
#      print t

IJ.run(dupimp, "RGB Color", "")
tlist = []
for t in tracks:
    if t[0].getFrame() == 1:
		tlist.append(t)
		for p in t:
			ip = dupimp.getStack().getProcessor(p.getFrame())
			if t[len(t)-1].getFrame() == dupimp.getStackSize():
				ip.setColor(Color.GREEN)
			else:
				ip.setColor(Color.RED)			
			ip.draw(p.retRoi())
        #print p.getFrame(), p.getSx(), p.getSy()
dupimp.show()

PTA.setNoGUI(False)
