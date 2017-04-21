'''
# @File(label="Select CSV file listing data folders and CellIDs") datalistfile
'''

# Dot detection based on MOSAIC feature point detector
# tracks detected dots over time with nearest neighbor method
# exports sum of "RippedOff" dots during 1-4th frame and 5-10th frame (drug addition at 4th frame)
# Kota Miura

from java.awt import Color
from ij import ImagePlus, IJ, ImageStack
from ij.io import RoiDecoder
from ij.process import StackStatistics
from ij.gui import Roi
from ij.plugin.frame import RoiManager
from mosaic.core.utils import MosaicUtils
from mosaic.core.detection import FeaturePointDetector
import os, csv, sys
import math


### key final variables ###

## particle detection parameters

# cutoff value for non-particle discrimination (float)
mosaic_cutoff = 0.0
# percentile for intensity threshold (float)
mosaic_percentile = 0.035
# estimated particle radius (int)
mosaic_radius = 2

## particle linking
# maximally allowed distance for particle being considered to be "stationary"
# unit = pixels. int value. 
criticalDist = 3
## frame numbers for sampling range. 
# As of 20170419 3 frames each from pre and post drug treatment 
# (drug is added just after 4th frame)
preDrug_Starts = 2
preDrug_Ends = 4
postDrug_Starts = 5
postDrug_Ends = 7
postDrug2_Starts = 8
postDrug2_Ends = 10


## particle class
class FP(object):
	
	def __init__(self, x, y, frame):
		self.x = x
		self.y = y
		self.frame = frame + 1

	def getX(self):
		return self.x
		
	def getY(self):
		return self.y
		
	def getFrame(self):
		return self.frame

	def retRoi(self):
		ww = 12
		tlx = self.x - math.floor(ww/2.0)
		tly = self.y - math.floor(ww/2.0)
		return Roi(tlx, tly, ww, ww)
		
def fetch_datalist(path):
	adict = {}
	gppath = ''
	with open(path, 'rb') as datalistcsv:
		datalist = csv.reader( datalistcsv )
		for i, row in enumerate(datalist):
			if i == 0:
				pass
			elif i ==1:
				gppath = row[0]
			else:
				templist =  filter(lambda x: x != '', row[2:])
				templist = map(int, templist)
				templist = map(lambda x: x -1, templist)
				adict[row[1]] = templist
	return gppath, adict

# p1 and p2 are particles
def dist(p1, p2):
	return math.sqrt(math.pow(p1.getX() - p2.getX(), 2) + math.pow(p1.getY() - p2.getY(), 2))


# searches a particle at similar position recursively
# fs frames, p a particle, track: a list
def searchNextParticle(fs, p, track):
	pn = None
	track.append(p)
	if  p.getFrame() < len(fs):
		ptcles = fs[p.getFrame()]  #particles in next frame
		dmin = 10
		for p2 in ptcles:
			dd = dist(p, p2)
			if dd < dmin:
				dmin = dd
				pn = p2
		if pn != None and dmin < criticalDist:
			searchNextParticle(fs, pn, track)
	else:
		print "search finished at frame", track[-1].getFrame() 



def core(cellroi, imagepath):
	#aroi = RoiDecoder.open(pp)
	aroi = cellroi
	print aroi
	maskhist = aroi.getMask().getHistogram( 256 )
	cellarea = maskhist[255]
	print "CellArea [pixels]:", cellarea
	
	imp = ImagePlus(imagepath)
	reallength = imp.getCalibration().pixelWidth
	totalframes = imp.getStackSize()
	particleArray = []
	frames = []
	for f in range(totalframes):
		particles = []
		currentframe = f + 1
		stack = MosaicUtils.GetSubStackInFloat(imp.getStack(), currentframe, currentframe) 
		detector = FeaturePointDetector(stack.getProcessor(1).getMax(), stack.getProcessor(1).getMin())
		detector.setDetectionParameters(mosaic_cutoff, mosaic_percentile, mosaic_radius, 0.5, False)
		detectedParticles = detector.featurePointDetection(stack)
		print 'particles:', len(detectedParticles)
		for p in detectedParticles:
			if aroi.contains(int(p.getX()), int(p.getY())):
				fp = FP(p.getX(), p.getY(), f)
				particles.append(fp)
		frames.append(particles)
		
	tracks = []
	trackingStartFrame = 4
	#for p in frames[ trackingStartFrame - 1 ]:
	for p in frames[ 0 ]:	
		track = []
		searchNextParticle(frames, p, track)
		tracks.append(track)
		
	for t in tracks:
		print len(t)
	
	# sum up number of particles ripped off during pre and post drug (from 20170419, 3 frames each)
	postcounts = 0
	precounts = 0
	postcounts2 = 0	
	for t in tracks:
		lastframe = t[-1].getFrame()
		if lastframe >= postDrug_Starts and lastframe <= postDrug_Ends:
			postcounts += 1
		elif lastframe >= preDrug_Starts and lastframe <= preDrug_Ends:
			precounts += 1
		elif lastframe >= postDrug2_Starts and lastframe <= postDrug2_Ends:
			postcounts2 += 1
	
	IJ.run(imp, "RGB Color", "")
	for t in tracks:
		#if len(t) == imp.getStackSize() - trackingStartFrame + 1:
		if len(t) == imp.getStackSize():
			cc = Color.GREEN
		else:
			cc = Color.RED
		for p in t:
			ip = imp.getStack().getProcessor(p.getFrame())
			ip.setColor(cc)
			ip.draw(p.retRoi())

	# exporting time-course of rip-off counts
	lastframe = totalframes
	countTrackLastPoints = [0] * ( lastframe )
	for t in tracks:
		if len(t) != lastframe :
			thislastFrame = t[-1].getFrame()
			countTrackLastPoints[thislastFrame -1] += 1
	# counting total counts in a frame
	countTrackTotalPoints = [0] * ( lastframe )
	for t in tracks:
		for ap in t:
			countTrackTotalPoints[ap.getFrame() - 1] += 1
	secondFrameVirusCounts = countTrackTotalPoints[1]

	outcountpath = imagepath + '_counts.csv'
	f = open(outcountpath, 'wb')
	writer = csv.writer(f)
	writer.writerow(['Frame', 'Counts', 'TotalCounts'])
	for index, count in enumerate(countTrackLastPoints):
		arow = [ index + 1, count, countTrackTotalPoints[index]]
		writer.writerow(arow)
	f.close()
	
	return tracks, cellarea, postcounts, postcounts2, precounts, imp, reallength, secondFrameVirusCounts


def main(parentpath, cellNo): 
	rootname = "cell" + str(cellNo)
	roizipname = 'RoiSet_' + rootname + '.zip'
	imagename = rootname + '_virus_median.tif'

	
	
	#pp = '/Users/miura/Desktop/161122 ctrl croped and 16 frames/RoiSet_cell5.zip'
	# unzipping http://stackoverflow.com/questions/3451111/unzipping-files-in-python
	#pp = '/Users/miura/Desktop/161122 ctrl croped and 16 frames/RoiSet_cell5/0005-0419-0327.roi'
	zippp = os.path.join(parentpath, roizipname)
	rm = RoiManager(False)
	rm.runCommand("Open", zippp)
	cellroi = rm.getRoi(2)
	if cellroi.getType() != 3:
		print "ROI type mismatch! ... ABORT"
		sys.exit()
	
	#imagepath = '/Users/miura/Desktop/161122 ctrl croped and 16 frames/cell1_virus_median.tif'
	imagepath = os.path.join(parentpath, imagename)
	
	tracks, cellarea, postcounts, postcounts2, precounts, imp, reallength, secondFrameVirusCounts = core(cellroi, imagepath)
 
	return tracks, cellarea, postcounts, postcounts2, precounts, imp, reallength, secondFrameVirusCounts


def batchProcess(parentpath, theExp):
	resarray = []
	#for ind in range(8):
	for ind in theExp:
		cellNo = ind + 1
		tracks, cellarea, postcounts, postcounts2, precounts, imp, reallength, secondFrameVirusCounts = main(parentpath, cellNo)
		
		scaledCellArea = cellarea * reallength * reallength 
		#density = len(tracks) / float(scaledCellArea)
		density = secondFrameVirusCounts / float(scaledCellArea)
		#ripoffRatio = postcounts / float(len(tracks))
		ripoffRatio = postcounts / float(secondFrameVirusCounts)
		ripoffRatio2 = postcounts2 / float(secondFrameVirusCounts)		
		print "Total number of detected dots: ", len(tracks)
		print "Second Frame Counts", secondFrameVirusCounts
		print "cell area [um2]", scaledCellArea
		print "Dot Density:[count / um2]:" , density
		print "Number of Ripped off (", preDrug_Starts, " to ", preDrug_Ends, " frame):", precounts	
		print "Number of Ripped off (", postDrug_Starts, " to ", postDrug_Ends, " frame):", postcounts
		print "Number of Ripped off (", postDrug2_Starts, " to ", postDrug2_Ends, " frame):", postcounts2
	
		#resarray.append([cellNo, scaledCellArea, len(tracks), density, precounts, postcounts, ripoffRatio])
		resarray.append([cellNo, scaledCellArea, secondFrameVirusCounts, density, precounts, postcounts, postcounts2, ripoffRatio])		
		outname = "cell" + str(cellNo) + '_dots.tif'
		savefilepath = os.path.join(parentpath, outname)	
		IJ.saveAsTiff(imp, savefilepath)
		
	# exporting results
	outcsvpath = os.path.join(parentpath, 'results.csv')
	f = open(outcsvpath, 'wb')
	writer = csv.writer(f)
	writer.writerow(['CellID', 'Area[um2]', 'Dots Total', 'Density', 'RipOff counts2_4', 'RipOff counts5_7', 'RipOff counts8_10', 'RipOff Ratio'])
	for arow in resarray:
		#arow = [cellNo, cellarea, len(tracks), postcounts, precounts]
		writer.writerow(arow)
	f.close()
#imp.show()  # this should be saved.

###############

datalistpath = '/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/data_lists/test.csv'

#grandparentpath, folderDict = fetch_datalist(datalistfile.getPath())
grandparentpath, folderDict = fetch_datalist(datalistpath)

print "Folder containing  data: ", grandparentpath
folders = folderDict.keys()
folders.sort()
for afolder in folders:
	print afolder
	print folderDict[afolder]


for key, val in folderDict.iteritems():
	parentpath = os.path.join(grandparentpath, key)
	batchProcess(parentpath, val)


