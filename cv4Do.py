import os
import json
from pathlib import Path
import sys
import numpy as np
from math import *
import pickle
import matplotlib.pyplot as plt

from cv4ProcessFile import CV4_PROCESS_FILE
from cv4AnalyzeWaveform import CV4_ANALYZE_WAVEFORM
from cv4AnalyzeDacScan import CV4_ANALYZE_DACSCAN


def plotQuick(vals_x,vals_y,label=""):
    fig, axes = plt.subplots(1,1,figsize=(13, 6))
    axes.plot(vals_x,vals_y,".")
    axes.set_xlabel('AWG AMP [V]', fontsize=20)
    axes.set_ylabel('ENOB', fontsize=20)
    axes.set_title("COLUTA ENOB VS AWG AMP", fontsize=20)
    axes.tick_params(axis="x", labelsize=12)
    axes.tick_params(axis="y", labelsize=12)
    #axes[0].set_xlim(0,4000)
    #axes[0].set_ylim(17087-50,17087+50)    
    axes.grid()
    axes.set_ylim(6,12)    
    fig.suptitle(label, fontsize=16)
    fig.tight_layout()
    
    plt.show()

def main():
  print("HELLO")
  #chanName = "channel3"
  #runResultsDict  = pickle.load( open( "output_cv4ProcessFile_Run_0920_Output.hdf5.pickle", "rb" ) )

  if len(sys.argv) != 3 :
    print("ERROR, program requires filename + channel arguments")
    return
  fileName = sys.argv[1]
  #chanName = "channel1"
  chanName = str(sys.argv[2])
  runResultsDict  = pickle.load( open( fileName, "rb" ) )

  if runResultsDict["results"] == None :
    print("NO RESULTS")
    return
  
  print("HERE")
  
  if True :
    cv4AnalyzeWaveform = CV4_ANALYZE_WAVEFORM("")
    cv4AnalyzeWaveform.runResultsDict = runResultsDict
    mdacWeights = [4196.59, 4197.03, 4197.74, 4200.05, 4200.08, 4202.27, 4205.72, 4204.92]
    sarWeights = [3476.50, 1986.45, 994.19, 621.88, 373.28, 249.05, 124.10, 214.03, 123.09, 61.18, 30.77, 22.86, 15.26, 9.55, 5.75, 4.00, 2.00, 1.00, 0.50, 0.25]
    cv4AnalyzeWaveform.analyzeSample.mdacWeights = [x*1. for x in mdacWeights]
    cv4AnalyzeWaveform.analyzeSample.sarWeights =  [x*1 for x in sarWeights]
    
    plot_x = []
    plot_y = []
    for measNum in runResultsDict["results"] :
      awgAmp = float(cv4AnalyzeWaveform.runResultsDict["results"][measNum]["attrs"]['awgAmp'])
      #cv4AnalyzeWaveform.printEnob(chId="channel2",measNum=measNum)
      meanvals, stdvals, maxval,minval,enob = cv4AnalyzeWaveform.viewWaveform(chId=chanName,measNum=measNum,doPrint=True,doPlot=True)
  
  if False :
    cv4AnalyzeWaveform = CV4_ANALYZE_WAVEFORM("")
    cv4AnalyzeWaveform.runResultsDict = runResultsDict
    #run 938
    #FMINSEARCH
    #mdacWeights = [4151.21, 4152.23, 4152.10, 4153.28, 4151.33, 4151.85, 4152.88, 4151.38]
    #sarWeights = [3476.50, 1986.84, 994.48, 621.99, 373.48, 249.33, 124.25, 212.32, 121.43, 60.96, 30.76, 22.80, 15.21, 9.52, 5.75, 4.00, 2.00, 1.00, 0.50, 0.25]
    #run 1029
    #FMINSEARCH
    #mdacWeights = [4153.52, 4152.21, 4150.56, 4150.31, 4147.08, 4146.32, 4145.69, 4145.56]
    #sarWeights = [3476.50, 1986.84, 994.48, 621.99, 373.48, 249.33, 124.25, 212.32, 121.43, 60.96, 30.76, 22.80, 15.21, 9.52, 5.75, 4.00, 2.00, 1.00, 0.50, 0.25]
    #run 1031
    mdacWeights = [4196.59, 4197.03, 4197.74, 4200.05, 4200.08, 4202.27, 4205.72, 4204.92]
    sarWeights = [3476.50, 1986.45, 994.19, 621.88, 373.28, 249.05, 124.10, 214.03, 123.09, 61.18, 30.77, 22.86, 15.26, 9.55, 5.75, 4.00, 2.00, 1.00, 0.50, 0.25]
    cv4AnalyzeWaveform.analyzeSample.mdacWeights = [x*1. for x in mdacWeights]
    cv4AnalyzeWaveform.analyzeSample.sarWeights =  [x*1 for x in sarWeights]
    
    plot_x = []
    plot_y = []
    for measNum in runResultsDict["results"] :
      awgAmp = float(cv4AnalyzeWaveform.runResultsDict["results"][measNum]["attrs"]['awgAmp'])
      #cv4AnalyzeWaveform.printEnob(chId="channel2",measNum=measNum)
      meanvals, stdvals, maxval,minval,enob = cv4AnalyzeWaveform.viewWaveform(chId=chanName,measNum=measNum,doPrint=True,doPlot=True)
      plot_x.append(awgAmp)
      plot_y.append(enob)
      #print( measNum ,resultsDict[measNum])
      #print("\n")
      continue
    plotQuick(plot_x,plot_y,label="COLUTA Ch1,5.005MHz Sine Wave,8000 Samples")
  #return None   
  
  if False : #DAC scan
    cv4AnalyzeDacScan = CV4_ANALYZE_DACSCAN("")
    cv4AnalyzeDacScan.runResultsDict = runResultsDict
    
    #run 920 weights fmin
    mdacWeights = [4054.67, 4054.36, 4054.83, 4055.49, 4055.18, 4053.4, 4054.35, 4054.09]
    sarWeights = [3310.33, 1891.79, 946.251, 591.556, 355.171, 237.098, 118.073, 202.343, 115.931, 57.8463, 29.0422, 21.4246, 14.283, 9.04592, 5.47516, 3.80881, 1.9044, 0.952202, 0.476101, 0.238051]
    cv4AnalyzeDacScan.analyzeSample.mdacWeights = [x*1. for x in mdacWeights]
    cv4AnalyzeDacScan.analyzeSample.sarWeights =  [x*1 for x in sarWeights]
    
    cv4AnalyzeDacScan.analyzeSample.applyDdpuCorr = False
    cv4AnalyzeDacScan.analyzeSample.sampScaleFactor = 1.
    cv4AnalyzeDacScan.analyzeSample.sampOffsetVal = 3780
    cv4AnalyzeDacScan.analyzeSample.limitSamples = True
    cv4AnalyzeDacScan.analyzeSample.maxNumSamples = 1000
    #cv4AnalyzeDacScan.plotDacLinearityData(chId=chanName)
    cv4AnalyzeDacScan.measureDnl(chId=chanName)
    
  #simpleTest(cv4ProcessFile.runResultsDict)
  
  return

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
