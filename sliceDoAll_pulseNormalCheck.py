import json
import sys
import subprocess
import numpy as np
from math import *
import matplotlib.pyplot as plt

from sliceProcessFile import SLICE_PROCESS_FILE
from sliceProcessNormalData import SLICE_PROCESS_NORMALDATA
from sliceAnalyzeWaveform import SLICE_ANALYZE_WAVEFORM
from sliceAnalyzePulses import SLICE_ANALYZE_PULSES
from collections.abc import Iterable

def main():
  print("HELLO")
  if len(sys.argv) != 2:
    print("ERROR, program requires filename and channel name as arguments")
    return
  fileName = sys.argv[1]
  
  print("FILE",fileName)

  print("PROCESS FILE")
  processFile = SLICE_PROCESS_FILE(fileName)
  processFile.limitNumSamples = True
  processFile.maxNumSamples = 100000
  processFile.processFile()
  #processFile.dumpFile()

  if processFile.runResultsDict == None :
    return
  outputResults = processFile.runResultsDict
  print("DONE PROCESS FILE")
  
  print("PROCESS NORMAL DATA")
  processNormal = SLICE_PROCESS_NORMALDATA(fileName)
  processNormal.runResultsDict = outputResults
  processNormal.processFile()
  outputResults = processNormal.runResultsDict
  print("DONE PROCESS NORMAL DATA")
  if outputResults == None :
    return

  #define waveform processing class
  print("ANALYZE DATA")
  analyzeFile = SLICE_ANALYZE_WAVEFORM(fileName)
  analyzeFile.runResultsDict = outputResults
  analyzeFile.dropInitialSamples = False
  analyzeFile.applyDdpuCorr = False
  analyzeFile.numDroppedInitialSamples = 1

  #get channel data in each measurement
  fitResults = {}
  for measNum in analyzeFile.runResultsDict["results"] :
    #if measNum != "Measurement_070" : continue #optional
    measInfo = analyzeFile.runResultsDict["results"][measNum]
    measData = measInfo["data"]
    for chan in measData:
      chanInfo = measData[chan]
      for gain in chanInfo :
        print(measNum,chan,gain)
        #if chan != "channel079" : continue
        #if gain != "hi" : continue #optional
        #if chan == "channel064" or chan == "channel065" or chan == "channel068" or chan == "channel069" or chan == "channel072" or chan == "channel073" or chan == "channel076" or chan == "channel077" : continue
        #if chan == "channel048" or chan == "channel049" or chan == "channel052" or chan == "channel053" or chan == "channel056" or chan == "channel057" or chan == "channel060" or chan == "channel061" : continue
        chWf = analyzeFile.getMeasChData(measNum,chan,gain)
        if isinstance(chWf, Iterable) == False : continue
        chWfArray = np.array( chWf ,  dtype='u2')
    
        measAttrs = measInfo["attrs"]
        print(measAttrs)
        awg_amp = 0
        if "awg_amp" in measAttrs :
          awg_amp = measAttrs["awg_amp"]
        if 'att_val' in measAttrs :
          att_val = int(measAttrs['att_val'])
          if att_val > 0 :
            #convert AWG amp to signal current here
            awg_amp = float(awg_amp)/250.*1000./(  pow(10,att_val/20.)   )
            awg_amp = round(float(awg_amp),2)
        
        print("#Samples ",len(chWfArray) )
        if False :
          print("PLOT WAVEFORM")
          vals_x = [cnt for cnt in range(0,1000,1)]
          analyzeFile.plotVals(vals_x=vals_x,vals_y=chWfArray[0:1000],label=chan +" " + gain)

        #average pulse analysis
        if True :
          print("ANALYZE PULSE")
          analyzePulse = SLICE_ANALYZE_PULSES(fileName)   
          fitPulseSetDict = analyzePulse.findPulses(chWfArray)
          if len(fitPulseSetDict) == 0 :
            print("NO PULSES FOUND!")
            break
          result = analyzePulse.getPulses(chWfArray,fitPulseSetDict)
          if result == None : continue
          fitPulseSetDict, binsDict, pedVal, pedRms, avgPulseHeight = result
          print( "PED",pedVal,"RMS",pedRms,"Pulse Height",avgPulseHeight)
          if True : #control drawing pulse
            plotTitle = str(chan) + ", " + str(gain) + ", AMP=" + str(awg_amp)  + "mA"
            analyzePulse.plotPulses(fitPulseSetDict,binsDict,chWfBin=[],plotTitle=plotTitle)
            #minDiff = analyzePulse.plotAvgPulseDiffs(fitPulseSetDict,binsDict,plotTitle=plotTitle)
            #maxResid = analyzePulse.plotPulseResiduals(fitPulseSetDict,binsDict,plotTitle=plotTitle)
          if False : #control fit analysis
            fitAmps = analyzePulse.fitPulses(fitPulseSetDict,binsDict,pedVal,pedRms)
            if len(fitAmps) != 0 :
              print( "FIT AMPS",np.mean(fitAmps),np.std(fitAmps),len(fitAmps))
              fitResults[measNum] = {"awg_amp":awg_amp,"fit_mean":np.mean(fitAmps),"fit_rms":np.std(fitAmps)}  
        #
  #end measurement loop
  print("FIT RESULTS",fitResults)
  
  return

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
