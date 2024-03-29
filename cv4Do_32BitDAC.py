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
    return

def getHist(vals):
    bins = np.arange(int(np.min(vals))-0.5, int(np.max(vals))-0.5+1, 1)
    centers = (0.5*(bins[1:]+bins[:-1]))
    #centers = [int(x) for x in centers]
    histTotal, edges = np.histogram(vals,bins=bins)
    return histTotal,edges,centers

def plotHist(centers, histTotal):
    fig, axes = plt.subplots(1,1,figsize=(13, 6))
    axes.bar(centers, histTotal, width=0.8)
    axes.set_xlabel('COLUTA CODE [ADC]', fontsize=20)
    axes.set_ylabel('Number Entries', fontsize=20)
    axes.tick_params(axis="x", labelsize=12)
    axes.tick_params(axis="y", labelsize=12)
    axes.set_title("Sample Distribution", fontsize=20)
    #fig.suptitle(label, fontsize=20)
    fig.tight_layout()
    plt.show()

def getFlippedSamples(cv4AnalyzeSample,colutaWf,maxNum=None,flipSarBitNum=0):
    #if isinstance(colutaWf, Iterable) == False: 
    #  return []
    vals = []
    #print ("SAR BITS USED",self.sarWeights)
    #print(len(colutaWf), colutaWf)
    if len(colutaWf) == 0 :
      return []
    for samp in colutaWf :
      sampBin = cv4AnalyzeSample.convertIntTo32BitWord(samp)
      header = sampBin[0:2]
      clockCheck = sampBin[2:4]
      mdacBits = sampBin[4:12]
      sarBits = sampBin[12:32]
      sarBitList = [int(x) for x in sarBits]
      mdacBitList = [int(x) for x in mdacBits]
      if sarBitList[19-flipSarBitNum] == 1 :
            sarBitList[19-flipSarBitNum] = 0
      else :
            sarBitList[19-flipSarBitNum] = 1
      val = cv4AnalyzeSample.getColutaSampleValue(sarBitList,mdacBitList)
      if cv4AnalyzeSample.applyDdpuCorr == True :
        val = int(val*cv4AnalyzeSample.sampScaleFactor/4. - cv4AnalyzeSample.sampOffsetVal)
        #val = int(val*self.sampScaleFactor - self.sampOffsetVal)
      if cv4AnalyzeSample.dropOverFlowSamples == True and val > cv4AnalyzeSample.maxCode :
        val = 32767
      if cv4AnalyzeSample.dropOverFlowSamples == True and val < 0 :
        val = 0
      vals.append(val)
      #if maxNum != None:
      if cv4AnalyzeSample.limitSamples :
        #if len(vals) > maxNum :
        if len(vals) >= cv4AnalyzeSample.maxNumSamples :
          break
      if False and len(vals) < 10 :
        print("\t",samp,sampBin,header,clockCheck,mdacBits,sarBits,val)
    return vals
    
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
        cv4AnalyzeWaveform.analyzeSample.applyDdpuCorr = True
        cv4AnalyzeWaveform.analyzeSample.sampOffsetVal = 4096
        for measNum in runResultsDict["results"] :
            mdacWeights = runResultsDict["results"][measNum]["attrs"]['MDACConstsDdpu']
            mdacSum = 0
            for mdacNum in range(1,8,1):
                mdacSum = mdacSum + mdacWeights[mdacNum-1]
                mdacWeights[mdacNum] = mdacWeights[mdacNum] - mdacSum
            mdacWeights[0] = 4096*4 + mdacWeights[0]
            #mdacWeights = [x/4. for x in mdacWeights]
            mdacWeights = mdacWeights[::-1]
            #mdacWeights = runResultsDict["results"][measNum]["attrs"]['MDACConsts']
            #mdacWeights[0] = 4096 + mdacWeights[0]
            #mdacWeights = mdacWeights[::-1]
            sarWeights = runResultsDict["results"][measNum]["attrs"]['SARConstsDdpu']
            #sarWeights = [x/4. for x in sarWeights]
            #print( mdacWeights )
            #print( sarWeights )
            cv4AnalyzeWaveform.analyzeSample.mdacWeights = mdacWeights
            cv4AnalyzeWaveform.analyzeSample.sarWeights = sarWeights
            break

        for measNum in runResultsDict["results"] :
            if measNum != "Measurement_385" : continue
            chWf32Bit = cv4AnalyzeWaveform.getMeasChData(chId=chanName,measNum=measNum,get32Bit=True)
            chWfNormal = cv4AnalyzeWaveform.getMeasChData(chId=chanName,measNum=measNum,get32Bit=False)
            print(measNum)
            dacAVal = cv4AnalyzeWaveform.runResultsDict["results"][measNum]["attrs"]["dacAVal"]
            dacBVal = cv4AnalyzeWaveform.runResultsDict["results"][measNum]["attrs"]["dacBVal"]
            modWfNormal = getFlippedSamples(cv4AnalyzeWaveform.analyzeSample,chWf32Bit,maxNum=None,flipSarBitNum=19)
            histTotal,edges,centers = getHist(chWfNormal)
            modHistTotal,modEdges,modCenters = getHist(modWfNormal)
            #plotHist(centers,histTotal)
            #plotHist(modCenters,modHistTotal)
            print(dacAVal,dacBVal,hex(chWf32Bit[0]),chWfNormal[0],modWfNormal[0])
            #break

    return

#-------------------------------------------------------------------------
if __name__ == "__main__":
    main()
