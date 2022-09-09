import numpy as np
from math import *
import matplotlib.pyplot as plt
import sys
import glob
import pickle

from cv4AnalyzeSample import CV4_ANALYZE_SAMPLE

from scipy.stats import norm

def model_pdf(x, mu, sigma):
    return norm.pdf(x, mu, sigma)

def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

def gaussian(x, amp, cen, wid):
    return amp * exp(-1.*(x-cen)*(x-cen) / (2*wid*wid))

def close_event():
    plt.close() #timer calls this function after 3 seconds and closes the window 

#BEGIN SLICE_ANALYZE_WAVEFORM CLASS
class CV4_ANALYZE_WAVEFORM(object):

  #__INIT__#
  def __init__(self, fileName = None):
    self.runResultsDict = None
    self.fileNameDict = None
    self.fileName = fileName
    self.dropInitialSamples = False
    self.numDroppedInitialSamples = 1
    self.applyCuts = False
    self.applyDdpuCorr = False
    self.dropOverFlowSamples = False
    
    #sample handling object
    self.analyzeSample = CV4_ANALYZE_SAMPLE()
    self.analyzeSample.dropOverFlowSamples = False
    self.analyzeSample.applyDdpuCorr = False
    self.analyzeSample.sampScaleFactor = 1.0
    self.analyzeSample.sampOffsetVal = 0.

  def SINAD(self,fourier):
    sum2 = 0
    for normBin in fourier:
      if normBin==1: continue
      sum2 += normBin**2
    return -10*np.log10(sum2)

  def ENOB(self,fourier):
    return (self.SINAD(fourier)-1.76)/6.02

  def getFftWaveform(self,vals):
    #print("num Samples",len(vals))
    fft_wf = np.fft.fft(vals)
    fftWf_x = []
    fftWf_y = []
    psd = []
    psd_x = []
    for sampNum,samp in enumerate(fft_wf) :
      if sampNum > float( len(fft_wf) ) / 2. :
        continue
      freqVal = 40. * sampNum/float(len(fft_wf))
      sampVal = np.abs(samp)
      if sampNum == 0 :
        sampVal = 0
      fftWf_x.append(freqVal)
      fftWf_y.append(sampVal)
    if np.max(fftWf_y) <= 0 :
      return psd_x,psd,0,0

    #fourier_fftWf_y = fftWf_y/np.max(fftWf_y)
    fourier_fftWf_y = fftWf_y
    for sampNum,samp in enumerate(fourier_fftWf_y) :
      if sampNum == 0 :
        continue
      else:
        psd_x.append( fftWf_x[sampNum] )
        psd.append( 20*np.log10(samp) )
    sinad = self.SINAD(fourier_fftWf_y)
    enob = self.ENOB(fourier_fftWf_y)
    #print ("ENOB\t",enob)
    return psd_x,psd,sinad,enob

    coarseVal = 1
    freqBinDict = {}
    
    if len(fftWf_x) != len(fftWf_y) :
      return psd_x,psd,sinad,enob
    
    for freqNum, freq in enumerate( fftWf_x ):
      freqBin = int(freqNum/coarseVal)
      freqVal = fftWf_x[freqNum]
      fftVal = fftWf_y[freqNum]
      if freqBin not in freqBinDict :
        freqBinDict[freqBin] = {"freqVal":[],"fftVal":[]}
      freqBinDict[freqBin]["freqVal"].append(freqVal)
      freqBinDict[freqBin]["fftVal"].append(fftVal)
      #print(freqNum,"\t",freq,"\t",freqBinDict[freqBin]["freqVal"])
    
    new_fftWf_x = []
    new_fftWf_y = []
    for freqBin in freqBinDict :
      #print(freqBin)
      #print("\t", freqBinDict[freqBin]["freqVal"],"\t",freqBinDict[freqBin]["fftVal"] )
      new_fftWf_x.append( np.mean(freqBinDict[freqBin]["freqVal"]) )
      new_fftWf_y.append( np.mean(freqBinDict[freqBin]["fftVal"]) )
    fftWf_x = new_fftWf_x[:-1]
    fftWf_y = new_fftWf_y[:-1]

    #fourier_fftWf_y = fftWf_y/np.max(fftWf_y)
    fourier_fftWf_y = fftWf_y
    #fourier_fftWf_y = fftWf_y/np.sum(fftWf_y)
    for sampNum,samp in enumerate(fourier_fftWf_y) :
      if sampNum == 0 :
        continue
      else:
        psd_x.append( fftWf_x[sampNum] )
        psd.append( 20*np.log10(samp) )
        #psd.append( samp )
    #sinad = self.SINAD(fourier_fftWf_y)
    #enob = self.ENOB(fourier_fftWf_y)
    return psd_x,psd,sinad,enob

  def printEnob(self,chId=None,measNum=None):
    vals = self.getMeasChData(chId=chId,measNum=measNum)
    if len(vals) == 0 :
      return None
    print(len(vals))
    psd_x,psd,sinad,enob = self.getFftWaveform(vals)
    print(measNum,"\t",enob)

  def plotVals(self,vals_x=[],vals_y=[],label=""):
    
    print("MEAN\t",round(np.mean(vals_y),2),"\tRMS\t",round(np.std(vals_y),2))
    
    fig, axes = plt.subplots(1,1,figsize=(13, 6))
    axes.plot(vals_x,vals_y,".")
    axes.set_xlabel('Sample #', fontsize=20)
    axes.set_ylabel('ADC CODE [ADC]', fontsize=20)
    axes.set_title("COLUTA WAVEFORM", fontsize=20)
    axes.tick_params(axis="x", labelsize=12)
    axes.tick_params(axis="y", labelsize=12)
    #axes[0].set_xlim(0,4000)
    #axes[0].set_ylim(17087-50,17087+50)    
    fig.suptitle(label, fontsize=16)
    fig.tight_layout()
    
    plt.show()
    #plt.draw()
    return None

  def plotValsFft(self,measNum="",vals_x=[],vals_y=[],psd_x=[],psd=[],label=""):
    #return
    
    fig, axes = plt.subplots(1,2,figsize=(13, 6))
    #fig, axes = plt.subplots(1,3,figsize=(13, 6))
    #timer = fig.canvas.new_timer(interval = 1000) #creating a timer object and setting an interval of 3000 milliseconds
    #timer.add_callback(close_event)
    
    #axes[0].plot(vals[0:100],".")
    axes[0].plot(vals_x,vals_y,".")
    axes[0].set_xlabel('Sample #', fontsize=20)
    axes[0].set_ylabel('ADC CODE [ADC]', fontsize=20)
    axes[0].set_title("COLUTA WAVEFORM", fontsize=20)
    axes[0].tick_params(axis="x", labelsize=12)
    axes[0].tick_params(axis="y", labelsize=12)
    #axes[0].set_xlim(0,4000)
    #axes[0].set_xlim(3600,6400)
    #axes[0].set_xlim(16077,16094)
    #axes[0].set_ylim(17087-50,17087+50)
    axes[1].plot(psd_x,psd,"-")
    axes[1].set_xlabel('Frequency [MHz]', horizontalalignment='right', fontsize=20)
    axes[1].set_ylabel('PSD [dB]', horizontalalignment='left', fontsize=20)
    axes[1].set_title("PSD", fontsize=20)
    axes[1].tick_params(axis="x", labelsize=12)
    axes[1].tick_params(axis="y", labelsize=12)
    #axes[1].set_xlim(0.,2.5)
    #axes[1].set_xlim(0.756-0.02+0.2544*2,0.756+0.02+0.2544*2)
    
    if False :
      bins = np.arange(int(np.min(vals))-0.5, int(np.max(vals))-0.5+1, 1)
      centers = (0.5*(bins[1:]+bins[:-1]))
      histTotal, edges = np.histogram(vals,bins=bins)

      #axes[2].plot(psd_x,psd,"-")
      axes[2].bar(centers, histTotal, width=0.8)
      axes[2].set_xlabel('COLUTA CODE [ADC]', fontsize=20)
      axes[2].set_ylabel('Number Entries', fontsize=20)
      axes[2].tick_params(axis="x", labelsize=12)
      axes[2].tick_params(axis="y", labelsize=12)
      axes[2].set_title("Sample Distribution", fontsize=20)
    
    
    fig.suptitle(label, fontsize=20)
    fig.tight_layout()
    
    #timer.start()
    plt.show()
    #plt.draw()

  def getMeasChData(self,chId=None,measNum=None,get32Bit=False):
    if chId == None or self.runResultsDict == None or measNum == None :
      print("ERROR, INCORRECT INPUT",chId,measNum, self.runResultsDict)
      return None
    if "results" not in self.runResultsDict:
      print("ERROR, MISSING RESULTS")
      return None
    runData = self.runResultsDict["results"]
    if measNum not in runData :
      print("ERROR, MISSING MEASNUM",measNum)
      return None
    measInfo = runData[measNum]
    if "data" not in measInfo:
      print("ERROR, MISSING DATA")
      return None
    measData = measInfo["data"]
    if chId not in measData:
      print("ERROR, MISSING CHID",chId)
      return None
    chanData = measData[chId]
    if "attrs" not in measInfo:
      print("ERROR, MISSING ATTRS")
      return None
    measAttrs = measInfo["attrs"]
    if '32BitMode' not in measAttrs:
      print("ERROR, MISSING ATTRS 32BITMODE")
      return None
    chWf = chanData["wf"]
    if measAttrs['32BitMode'] == chId and get32Bit == False:
      #print("DO 32BIT MODE CORRECTION")
      chWf = self.analyzeSample.getWaveformVals(chWf)
    #account for bad samples
    if self.dropInitialSamples :
      chWf = chWf[self.numDroppedInitialSamples:]
    return chWf


  def viewWaveform(self,chId=None,measNum=None,doPrint=False,doPlot=False):
    #chWf = self.getMeasChData(chId=chId,measNum=measNum)
    vals = self.getMeasChData(chId=chId,measNum=measNum)
    #vals == None :
    #  return None
    if len(vals) == 0 :
      return
    #print( self.getFftWaveform(vals) )
    psd_x,psd,sinad,enob = self.getFftWaveform(vals)
    #print("Number of samples ",len(chWf))
    if doPrint:
      print("MEAS#\t",measNum,"\tMEAN\t",np.mean(vals),"\tRMS\t",np.std(vals),"\tMIN\t",np.min(vals),"\tMAX\t",np.max(vals),"\tRANGE\t",np.max(vals)-np.min(vals) )
      print("\t", len(vals),"\t",enob )
      print("\t", self.runResultsDict["results"][measNum]["attrs"] )
      if False and self.runResultsDict["results"][measNum]["attrs"]["SARConstsDdpu"][0] > 0 :
        #print("\tSAR", [ round(x /self.runResultsDict["results"][measNum]["attrs"]["SARConstsDdpu"][0] *3476.50,2) for x in self.runResultsDict["results"][measNum]["attrs"]["SARConstsDdpu"] ] )
        #print("\tSAR", [ round(x /4.,2) for x in self.runResultsDict["results"][measNum]["attrs"]["SARConstsDdpu"] ] )
        MDACConstsDdpu = self.runResultsDict["results"][measNum]["attrs"]["MDACConstsDdpu"]
        convMdac = []
        prevVal = 0
        for corrNum,corr in enumerate(MDACConstsDdpu) :
          convMdac.append( corr/4. - prevVal )
          prevVal = corr/4.
        convMdac = convMdac[::-1]
        print("\tMDAC", convMdac )
      print("\n")
    #return
    if doPlot:
      vals_x = [num for num,x in enumerate(vals) ]      
      #self.plotVals(vals_x=vals_x[0:75],vals_y=vals[0:75],label="")
      #self.plotVals(vals_x=vals_x,vals_y=vals,label="")
      self.plotValsFft(measNum="",vals_x=vals_x[0:200],vals_y=vals[0:200],psd_x=psd_x,psd=psd,label="")
    #return np.mean(vals)
    return np.mean(vals), np.std(vals), np.max(vals),np.min(vals),enob


  def dumpFile(self):
    if self.runResultsDict == None:
      return
    if "results" not in self.runResultsDict:
      return
    runData = self.runResultsDict["results"]
    for measNum in runData:
      measData = runData[measNum]
      print("Measurement ", measNum)
      print("\t",measData)
    return

  #open file
  def openFile(self):
    if self.fileName == None :
      print("ERROR no input file specified")
      return None
    self.runResultsDict = pickle.load( open( self.fileName, "rb" ) )
    return

def main():
  if len(sys.argv) != 2 :
    print("ERROR, program requires filename argument")
    return
  fileName = sys.argv[1]
  cv3tbAnalyzeFile = CV3TB_ANALYZE_WAVEFORM(fileName)
  
#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
