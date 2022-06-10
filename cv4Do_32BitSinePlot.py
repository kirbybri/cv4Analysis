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

def getPlotVals(cv4AnalyzeWaveform,runResultsDict,mdacWeights,sarWeights,chanName):
    cv4AnalyzeWaveform.analyzeSample.mdacWeights = [x*1. for x in mdacWeights]
    cv4AnalyzeWaveform.analyzeSample.sarWeights =  [x*1 for x in sarWeights]
    
    plot_x = []
    plot_y = []
    for measNum in runResultsDict["results"] :
      awgAmp = float(cv4AnalyzeWaveform.runResultsDict["results"][measNum]["attrs"]['awgAmp'])
      #cv4AnalyzeWaveform.printEnob(chId="channel2",measNum=measNum)
      doPlot = True
      #if awgAmp > 0.4 : doPlot = True
      meanvals, stdvals, maxval,minval,enob = cv4AnalyzeWaveform.viewWaveform(chId=chanName,measNum=measNum,doPrint=False,doPlot=doPlot)
      plot_x.append(awgAmp)
      plot_y.append(enob)
      #print( measNum ,resultsDict[measNum])
      #print("\n")
      continue
    return plot_x,plot_y

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
    #run 938
    #FMINSEARCH
    #mdacWeights = [4151.21, 4152.23, 4152.10, 4153.28, 4151.33, 4151.85, 4152.88, 4151.38]
    #sarWeights = [3476.50, 1986.84, 994.48, 621.99, 373.48, 249.33, 124.25, 212.32, 121.43, 60.96, 30.76, 22.80, 15.21, 9.52, 5.75, 4.00, 2.00, 1.00, 0.50, 0.25]
    #run 1029
    #FMINSEARCH
    #mdacWeights = [4153.52, 4152.21, 4150.56, 4150.31, 4147.08, 4146.32, 4145.69, 4145.56]
    #sarWeights = [3476.50, 1986.84, 994.48, 621.99, 373.48, 249.33, 124.25, 212.32, 121.43, 60.96, 30.76, 22.80, 15.21, 9.52, 5.75, 4.00, 2.00, 1.00, 0.50, 0.25]
    #run 1031
    #mdacWeights = [4196.59, 4197.03, 4197.74, 4200.05, 4200.08, 4202.27, 4205.72, 4204.92]
    #sarWeights = [3476.50, 1986.45, 994.19, 621.88, 373.28, 249.05, 124.10, 214.03, 123.09, 61.18, 30.77, 22.86, 15.26, 9.55, 5.75, 4.00, 2.00, 1.00, 0.50, 0.25]
    #run 1042
    #plotTitle = "Run1097,Ch1,18.005MHz Sine Wave,8000 Samples"
    #mdacWeights = [4216.61, 4201.68, 4202.03, 4202.95, 4201.57, 4202.17, 4202.25, 4218.72]
    #sarWeights = [3476.50, 1986.45, 994.19, 621.88, 373.28, 249.05, 124.10, 214.03, 123.09, 61.18, 30.77, 22.86, 15.26, 9.55, 5.75, 4.00, 2.00, 1.00, 0.50, 0.25 ]
    
    #run1097,1099, ch1
    #plotTitle = "Run1097,Ch1,5.005MHz Sine Wave,8000 Samples"
    #plotTitle = "Run1099,Ch1,5.005MHz Sine Wave,8000 Samples"
    #plotTitle = "Run1100,Ch1,5.005MHz Sine Wave,8000 Samples"
    #plotTitle = "Run1135,Ch1,18.005MHz Sine Wave,8000 Samples"
    #mdacWeights = [4200.49, 4201.84, 4201.23, 4203.01, 4201.35, 4201.15, 4202.49, 4201.03]
    #sarWeights = [3476.50, 1986.45, 994.19, 621.88, 373.28, 249.05, 124.10, 214.03, 123.09, 61.18, 30.77, 22.86, 15.26, 9.55, 5.75, 4.00, 2.00, 1.00, 0.50, 0.25]
    
    #run 1106,1107,1113 ch 2
    #plotTitle = "Run1106,Ch1,18.005MHz Sine Wave,8000 Samples"
    #plotTitle = "Run1107,Ch1,8.005MHz Sine Wave,8000 Samples"
    #plotTitle = "Run1113,Ch1,5.005MHz Sine Wave,8000 Samples"
    #plotTitle = "Run1131,Ch1,18.005MHz Sine Wave,8000 Samples"
    #mdacWeights = [4225.09, 4228.00, 4227.52, 4228.67, 4227.00, 4227.46, 4227.55, 4225.67]
    #sarWeights = [3476.50, 1986.19, 994.47, 621.63, 373.65, 249.45, 124.00, 214.50, 123.78, 61.23, 30.83, 22.76, 15.19, 9.29, 5.51, 3.61, 2.00, 1.00, 0.50, 0.25]
    
    #run 1119, board 153 ch4
    #plotTitle = "Run1119,Ch4,18.005MHz Sine Wave,8000 Samples"
    #plotTitle = "Run1120, Ch4, 8.005MHz Sine Wave, 8000 Samples"
    #plotTitle = "Run1121, Ch4, 5.005MHz Sine Wave, 8000 Samples"
    #mdacWeights = [4203.30, 4204.95, 4203.56, 4205.91, 4205.78, 4205.70, 4207.39, 4205.66]
    #sarWeights = [3476.50, 1986.81, 995.10, 622.51, 373.38, 249.32, 124.34, 213.77, 122.59, 61.27, 30.94, 22.87, 15.23, 9.57, 5.75, 4.00, 2.00, 1.00, 0.50, 0.25]
    plotTitle = "Run1150, Ch4, 5.005MHz Sine Wave, 8000 Samples"
    mdacWeights = [4216.61, 4218.00, 4181.83, 4186.88, 4193.65, 4188.97, 4218.72, 4218.72]
    sarWeights = [3476.50, 1986.45, 994.19, 621.88, 373.28, 249.05, 124.10, 214.03, 123.09, 61.18, 30.77, 22.86, 15.26, 9.55, 5.75, 4.00, 2.00, 1.00, 0.50, 0.25]
    print( mdacWeights )
    plot_x_fmin,plot_y_fmin = getPlotVals(cv4AnalyzeWaveform,runResultsDict,mdacWeights,sarWeights,chanName)
    
    for measNum in runResultsDict["results"] :
      mdacWeights = runResultsDict["results"][measNum]["attrs"]['MDACConsts']
      mdacWeights = mdacWeights[::-1]
      sarWeights = runResultsDict["results"][measNum]["attrs"]['SARConstsDdpu']
      sarWeights = [x/4. for x in sarWeights]
      print( mdacWeights )
      #print( sarWeights )
      break
    plot_x_const,plot_y_const = getPlotVals(cv4AnalyzeWaveform,runResultsDict,mdacWeights,sarWeights,chanName)

    fig, axes = plt.subplots(1,1,figsize=(13, 6))
    axes.plot(plot_x_fmin,plot_y_fmin,".",label="FMIN weights")
    axes.plot(plot_x_const,plot_y_const,".",label="On-chip weights")
    axes.set_xlabel('AWG AMP [V]', fontsize=20)
    axes.set_ylabel('ENOB', fontsize=20)
    axes.set_title("COLUTA ENOB VS AWG AMP", fontsize=20)
    axes.tick_params(axis="x", labelsize=12)
    axes.tick_params(axis="y", labelsize=12)
    #axes[0].set_xlim(0,4000)
    #axes[0].set_ylim(17087-50,17087+50)    
    axes.grid()
    axes.legend()
    axes.set_ylim(6,12)    
    fig.suptitle(plotTitle, fontsize=16)
    fig.tight_layout()
    plt.show()
  return None   
  

    
  #simpleTest(cv4ProcessFile.runResultsDict)
  
  return

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
