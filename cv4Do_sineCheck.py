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
  if len(sys.argv) != 2 :
    print("ERROR, program requires filename argument")
    return
  fileName = sys.argv[1]
  chanName = "channel1"
  #chanName = str(sys.argv[2])
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
  
  if True :
    cv4AnalyzeWaveform = CV4_ANALYZE_WAVEFORM("")
    cv4AnalyzeWaveform.runResultsDict = cv4ProcessFile.runResultsDict
    
    plot_x = []
    plot_y = []
    for measNum in cv4AnalyzeWaveform.runResultsDict["results"] :
      awgAmp = float(cv4AnalyzeWaveform.runResultsDict["results"][measNum]["attrs"]['awgAmp'])
      meanvals, stdvals, maxval,minval,enob,sinad = cv4AnalyzeWaveform.viewWaveform(chId=chanName,measNum=measNum,doPrint=False,doPlot=False)
      print( measNum , enob )
      plot_x.append(awgAmp)
      plot_y.append(enob)
    #plotQuick(plot_x,plot_y,label="COLUTA Ch1,5.005MHz Sine Wave,8000 Samples")
    print("AVG",np.mean(plot_y),"RMS",np.std(plot_y))
  return

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
