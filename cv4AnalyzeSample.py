import numpy as np
from math import *
from collections.abc import Iterable

class CV4_ANALYZE_SAMPLE(object):

  #__INIT__#
  def __init__(self):
    #calib constants
    #self.sarWeights = [3584*4,2048*4,1024*4,640*4,384*4,256*4,128*4,224*4,128*4,64*4,32*4,24*4,16*4,10*4,6*4,4*4,2*4,1*4,2,1]
    #self.mdacWeights = [4288*4, 4288*4, 4288*4, 4288*4, 4288*4, 4288*4, 4288*4, 4288*4]
    self.mdacWeights = [4288, 4288, 4288, 4288, 4288, 4288, 4288, 4288]
    self.sarWeights = [3584,2048,1024,640,384,256,128,224,128,64,32,24,16,10,6,4,2,1,0.5,0.25]
    
    self.dropOverFlowSamples = False
    self.applyDdpuCorr = False
    self.sampScaleFactor = 1.0 
    self.sampOffsetVal = 0.
    self.maxCode = 32767
    self.limitSamples = False
    self.maxNumSamples = 100

  def getColutaSampleValue(self,sarBitList=[],mdacBitList=[]):
    val = 0
    #for bitNum in range(len(sarBits)):
    #    val += self.sarWeights[bitNum]*int(sarBits[bitNum])
    #for bitNum in range(len(mdacBits)):
    #    val += self.mdacWeights[bitNum]*int(mdacBits[bitNum])
    val = np.dot(sarBitList,self.sarWeights) + np.dot(mdacBitList,self.mdacWeights)
    return val

  def convertIntTo16BitWord(self,sampInt):
    sampBin = bin(sampInt)[2:].zfill(16) #this is recovering the actual COLUTA data word which was stored as an int in the dict, this is silly
    return sampBin

  def convertIntTo32BitWord(self,sampInt):
    sampBin = bin(sampInt)[2:].zfill(32) #this is recovering the actual COLUTA data word which was stored as an int in the dict, this is silly
    return sampBin

  def getWaveformVals(self,colutaWf,maxNum=None):
    if isinstance(colutaWf, Iterable) == False: 
      return []
    vals = []
    #print ("SAR BITS USED",self.sarWeights)
    #print(len(colutaWf), colutaWf)
    if len(colutaWf) == 0 :
      return []
    for samp in colutaWf :
      sampBin = self.convertIntTo32BitWord(samp)
      header = sampBin[0:2]
      clockCheck = sampBin[2:4]
      mdacBits = sampBin[4:12]
      sarBits = sampBin[12:32]
      sarBitList = [int(x) for x in sarBits]
      mdacBitList = [int(x) for x in mdacBits]
      val = self.getColutaSampleValue(sarBitList,mdacBitList)
      if self.applyDdpuCorr == True :
        val = int(val*self.sampScaleFactor/4. - self.sampOffsetVal)
        #val = int(val*self.sampScaleFactor - self.sampOffsetVal)
      if self.dropOverFlowSamples == True and val > self.maxCode :
        val = 32767
      if self.dropOverFlowSamples == True and val < 0 :
        val = 0
      vals.append(val)
      #if maxNum != None:
      if self.limitSamples :
        #if len(vals) > maxNum :
        if len(vals) >= self.maxNumSamples :
          break
      if False and len(vals) < 10 :
        print("\t",samp,sampBin,header,clockCheck,mdacBits,sarBits,val)
    return vals
    #return [np.mean(vals)]

def main():
  return
#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
