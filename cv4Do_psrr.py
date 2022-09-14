import os
import json
from pathlib import Path
import sys
import numpy as np
from math import *
import matplotlib.pyplot as plt

from cv4ProcessFile import CV4_PROCESS_FILE
from cv4AnalyzeWaveform import CV4_ANALYZE_WAVEFORM
from cv4AnalyzeDacScan import CV4_ANALYZE_DACSCAN
from cv4AnalyzePulses import CV4_ANALYZE_PULSES

def simpleTest(resultsDict):
  cv4AnalyzeWaveform = CV4_ANALYZE_WAVEFORM("")
  cv4AnalyzeWaveform.runResultsDict = resultsDict
  for measNum in resultsDict["results"] :
    #cv4AnalyzeWaveform.printEnob(chId="channel2",measNum=measNum)
    #cv4AnalyzeWaveform.viewWaveform(chId="channel1",measNum=measNum,doPrint=True,doPlot=True)
    #print( measNum ,resultsDict[measNum])
    #print("\n")
    continue
  return None

def plotQuick(vals_x,vals_y,label=""):
    fig, axes = plt.subplots(1,1,figsize=(13, 6))
    axes.plot(vals_x,vals_y,".")
    #axes.set_xlabel('AWG AMP [V]', fontsize=20)
    #axes.set_xlabel('Input Signal Amplitude [mV diff]', fontsize=20)
    #axes.set_ylabel('ENOB', fontsize=20)
    #axes.set_title("COLUTA ENOB VS AWG AMP", fontsize=20)
    axes.tick_params(axis="x", labelsize=12)
    axes.tick_params(axis="y", labelsize=12)
    #axes.set_xlim(-2.5,2.5)
    #axes[0].set_ylim(17087-50,17087+50)    
    axes.grid()
    #axes.set_ylim(6,12)
    fig.suptitle(label, fontsize=16)
    fig.tight_layout()
    plt.show()

def findMaxPsd(chFFt_x,chFft_y,lowFreq,highFreq):
  maxPsd = None
  maxPsd_x = None
  for psdNum,freq in enumerate(chFFt_x) :
    if freq < lowFreq : continue
    if freq > highFreq : continue
    psdVal = chFft_y[psdNum]
    #print( psdNum , freq , psdVal )
    if maxPsd == None :
      maxPsd = psdVal
      maxPsd_x = freq
    if psdVal > maxPsd :
      maxPsd = psdVal
      maxPsd_x = freq
  #print("MAX",maxPsd,maxPsd_x)
  return maxPsd,maxPsd_x

def main():
  print("HELLO")
  if len(sys.argv) != 3 :
    print("ERROR, program requires filename argument")
    return
  fileName = sys.argv[1]
  #chanName = "channel1"
  chanName = str(sys.argv[2])
  testNum = None
  
  print("FILE",fileName)

  print("PROCESS FILE")
  cv4ProcessFile = CV4_PROCESS_FILE(fileName)
  cv4ProcessFile.limitNumSamples = False
  cv4ProcessFile.maxNumSamples = 128
  cv4ProcessFile.openFile()
  cv4ProcessFile.processFile()
  #cv4ProcessFile.dumpFile()
  #cv4ProcessFile.outputTextFiles()
  #cv4ProcessFile.outputFile()
  #return

  if cv4ProcessFile.runResultsDict["results"] == None :
    print("NO RESULTS")
    return
  
  #get average FFT for specified channel
  chFftList = [] 
  chFFt_x = None
  psdDiffs = []
  if True :
    cv4AnalyzeWaveform = CV4_ANALYZE_WAVEFORM("")
    cv4AnalyzeWaveform.runResultsDict = cv4ProcessFile.runResultsDict
    for measNum in cv4AnalyzeWaveform.runResultsDict["results"] :
      #if measNum != "Measurement_891" : continue
      #doPlot = True
      #meanvals, stdvals, maxval,minval,enob = cv4AnalyzeWaveform.viewWaveform(chId=chanName,measNum=measNum,doPrint=True,doPlot=doPlot)
      #return None
      vals = cv4AnalyzeWaveform.getMeasChData(chId=chanName,measNum=measNum)
      if len(vals) == 0 : continue
      sampFreq = 40.E+6 #Hz
      sampPeriod = 1/float(sampFreq) #s
      times = []
      for sampNum, samp in enumerate(vals) :
        sampTime = sampNum * sampPeriod #s
        times.append(sampTime)
      #plotQuick(times,vals,label="WF "+str(measNum)+" "+str(chanName))
      psd_x,psd,sinad,enob = cv4AnalyzeWaveform.getFftWaveform(vals)
      print(len(psd))
      
      if chFFt_x == None : chFFt_x = psd_x
      chFftList.append(psd)

      #carrier analysis
      if False :
        expFreq = 8
        distFreq = 0.1
        #find max tone in individual FFT
        maxPsd,maxPsd_x = findMaxPsd(psd_x,psd,expFreq-0.05,expFreq+0.05)
        #find sub tone
        maxSubPsd,maxSuPsd_x = findMaxPsd(psd_x,psd,expFreq+distFreq-0.05,expFreq+distFreq+0.05)
        diffPsd = maxPsd - maxSubPsd 
        psdDiffs.append(diffPsd)
      #plotQuick(psd_x,psd,label="FFT "+str(measNum)+" "+str(chanName))
  
  #print( "PSD DIFFS", np.mean(psdDiffs), "STD", np.std(psdDiffs) )
  #return None

  if len(chFftList) == 0 :
    return None
  chAvgFft = np.mean(chFftList,axis=0)

  #find tone
  """
  maxPsd = None
  maxPsd_x = None
  for psdNum,freq in enumerate(chFFt_x) :
    if freq < 0.09 : continue
    psdVal = chAvgFft[psdNum]
    #print( psdNum , freq , psdVal )
    if maxPsd == None :
      maxPsd = psdVal
      maxPsd_x = freq
    if psdVal > maxPsd :
      maxPsd = psdVal
      maxPsd_x = freq
  """
  #find disturbance impact
  
  expFreq = 0.1
  maxPsd,maxPsd_x = findMaxPsd(chFFt_x,chAvgFft,expFreq-0.045,expFreq+0.045)
  print("MAX",maxPsd,maxPsd_x)
  #find noise floor
  noiseFloorPsd = []
  for psdNum,freq in enumerate(chFFt_x) :
    if freq < 10 : continue
    if freq == maxPsd_x : continue
    psdVal = chAvgFft[psdNum]
    noiseFloorPsd.append(psdVal)
  print("NOISE FLOOR",np.mean(noiseFloorPsd),"STD",np.std(noiseFloorPsd))
  print("DIFFERENCE", maxPsd - np.mean(noiseFloorPsd) )
  plotQuick(chFFt_x,chAvgFft,label="REAL DATA AVG FFT "+str(chanName)+" FREQ "+str(expFreq) + "MHz")
  
  #carrier analysis
  """
  expFreq = 5
  distFreq = 0.1
  #find max tone in individual FFT
  maxPsd,maxPsd_x = findMaxPsd(chFFt_x,chAvgFft,expFreq-0.045,expFreq+0.045)
  print("MAX",maxPsd,maxPsd_x)
  #find sub tone
  maxSubPsd,maxSuPsd_x = findMaxPsd(chFFt_x,chAvgFft,expFreq-distFreq-0.04,expFreq-distFreq+0.04)
  print("LEFT",maxSubPsd,maxSuPsd_x)  
  #find sub tone
  maxSubPsd,maxSuPsd_x = findMaxPsd(chFFt_x,chAvgFft,expFreq+distFreq-0.04,expFreq+distFreq+0.04)
  print("RIGHT",maxSubPsd,maxSuPsd_x)  
  """
  #plotQuick(chFFt_x,chAvgFft,label="REAL DATA AVG FFT "+str(chanName))
  #simpleTest(cv4ProcessFile.runResultsDict)
  
  return

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
