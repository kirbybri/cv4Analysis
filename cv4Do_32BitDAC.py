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
    #axes.set_xlabel('AWG AMP [V]', fontsize=20)
    #axes.set_ylabel('ENOB', fontsize=20)
    #axes.set_title("COLUTA ENOB VS AWG AMP", fontsize=20)
    axes.tick_params(axis="x", labelsize=12)
    axes.tick_params(axis="y", labelsize=12)
    #axes[0].set_xlim(0,4000)
    #axes[0].set_ylim(17087-50,17087+50)    
    axes.grid()
    #axes.set_ylim(6,12)    
    fig.suptitle(label, fontsize=16)
    fig.tight_layout()
    
    plt.show()


def main():
  print("HELLO")
  if len(sys.argv) != 3 :
    print("ERROR, program requires filename + channel arguments")
    return
  fileName = sys.argv[1]
  chanName = str(sys.argv[2])
  runResultsDict  = pickle.load( open( fileName, "rb" ) )

  if runResultsDict["results"] == None :
    print("NO RESULTS")
    return
  print("HERE")
  
  if True :
    cv4AnalyzeWaveform = CV4_ANALYZE_WAVEFORM("")
    cv4AnalyzeWaveform.runResultsDict = runResultsDict
    cv4AnalyzeWaveform.applyDdpuCorr = True

    for measNum in runResultsDict["results"] :
      mdacWeights = runResultsDict["results"][measNum]["attrs"]['MDACConstsDdpu']
      mdacWeights[0] = 4096*4 + mdacWeights[0]
      #mdacWeights = [x/4. for x in mdacWeights]
      mdacWeights = mdacWeights[::-1]
      sarWeights = runResultsDict["results"][measNum]["attrs"]['SARConstsDdpu']
      #ssarWeights = [x/4. for x in sarWeights]
      print( mdacWeights )
      print( sarWeights )
      break
    
    cv4AnalyzeWaveform.analyzeSample.mdacWeights = mdacWeights
    cv4AnalyzeWaveform.analyzeSample.sarWeights = sarWeights

    for measNum in runResultsDict["results"] :
        chWf = cv4AnalyzeWaveform.getMeasChData(chId=chanName,measNum=measNum,get32Bit=True)
        chWfNormal = cv4AnalyzeWaveform.getMeasChData(chId=chanName,measNum=measNum,get32Bit=False)
        print(measNum)
        dacAVal = cv4AnalyzeWaveform.runResultsDict["results"][measNum]["attrs"]["dacAVal"]
        dacBVal = cv4AnalyzeWaveform.runResultsDict["results"][measNum]["attrs"]["dacBVal"]
        print(dacAVal,dacBVal,chWf[0],chWfNormal[0])


  return

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
