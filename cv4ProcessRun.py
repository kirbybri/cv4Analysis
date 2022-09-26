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
    axes.set_xlim(-2.5,2.5)
    #axes[0].set_ylim(17087-50,17087+50)    
    axes.grid()
    #axes.set_ylim(6,12)
    fig.suptitle(label, fontsize=16)
    fig.tight_layout()
    plt.show()
  
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
    
  if False :
    cv4AnalyzeWaveform = CV4_ANALYZE_WAVEFORM("")
    cv4AnalyzeWaveform.runResultsDict = cv4ProcessFile.runResultsDict
    cv4AnalyzePulse = CV4_ANALYZE_PULSES()

    if True :
        measNum = list(cv4AnalyzeWaveform.runResultsDict["results"].keys())[0]
        SARConstsDdpu  = cv4AnalyzeWaveform.runResultsDict["results"][measNum]["attrs"]["SARConstsDdpu"]
        MDACConstsDdpu = cv4AnalyzeWaveform.runResultsDict["results"][measNum]["attrs"]["MDACConstsDdpu"]
        convSar  = [ x/4. for x in SARConstsDdpu]
        convMdac = []
        prevVal = 0
        for corrNum,corr in enumerate(MDACConstsDdpu) :
          convMdac.append( corr/4. - prevVal )
          prevVal = corr/4.
        convMdac = convMdac[::-1]
        convMdac[7] = 4096 + convMdac[7]
        cv4AnalyzeWaveform.analyzeSample.mdacWeights = convMdac
        cv4AnalyzeWaveform.analyzeSample.sarWeights = convSar
        print("MDAC",convMdac)
        print("SAR",convSar)
    
    measResults = {}
    for measNum in cv4AnalyzeWaveform.runResultsDict["results"] :
      vals = cv4AnalyzeWaveform.getMeasChData(chId=chanName,measNum=measNum)
      vals32bit = cv4AnalyzeWaveform.getMeasChData(chId=chanName,measNum=measNum,get32Bit=True)
      fitPulseSetDict = cv4AnalyzePulse.findPulses(vals)
      if len(fitPulseSetDict) == 0 :
        print("NO PULSES FOUND!")
        continue
      result = cv4AnalyzePulse.getPulses(vals,fitPulseSetDict,chWfBin=vals32bit)
      if result == None : 
        print("NO PULSES SAVED")
        continue
      measResults[measNum] = result
      
    if False :
      print(cv4AnalyzePulse.slopes)
      print("SLOPE", np.mean(cv4AnalyzePulse.slopes), np.std(cv4AnalyzePulse.slopes) )
      print(cv4AnalyzePulse.intercepts)
      print("INTERCEPT", np.mean(cv4AnalyzePulse.intercepts), np.std(cv4AnalyzePulse.intercepts) )
      plot_x = []
      plot_y = []
      for intNum,intVal in enumerate(cv4AnalyzePulse.intercepts) :
        plot_x.append(intNum)
        plot_y.append(intVal)
      fig, axes = plt.subplots(1,1,figsize=(13, 8))
      axes.plot(plot_x,plot_y,".",label="Start")
      axes.set_xlabel('Measurement #', fontsize=20)
      axes.set_ylabel('Start Time [ns]', fontsize=20)
      axes.set_title("Pulse Train Start Time vs Measurement", fontsize=20)
      #fig.suptitle("Pulse Shape", fontsize=20)
      #axes.set_xlim(-75,75)
      axes.legend(fontsize=12)
      fig.tight_layout()
      plt.show()
      return
    
    plot_x = []
    plot_y = []
    plot_x_same = []
    plot_y_same = []
    mdacSamples = {}
    sampA = []
    sampB = []
    for measNum in measResults :
      #if measNum != "Measurement_146" : continue
      fitPulseSetDict = measResults[measNum]
      for pulseSet in fitPulseSetDict :
        #print("PULSE SET ", pulseSet)
        pulseSetInfo = fitPulseSetDict[pulseSet]
        for pulseNum in pulseSetInfo :
          #if pulseNum < 4 : continue
          #if pulseNum > 29 + 4 : continue
          pulseTimeInfo = pulseSetInfo[pulseNum]
          print(measNum,pulseSet,pulseNum,pulseTimeInfo.keys())
          prevMdacBits = None
          prevBits = None
          #isGood = False
          #for sampNum,sampTime in enumerate(pulseTimeInfo["pulse_fit_x"]):
          #  if sampTime*25 > -40 and sampTime*25 < -30 :  isGood = True
          #if isGood == False : continue
          for sampNum,sampTime in enumerate(pulseTimeInfo["pulse_fit_x"]):
            sampVal = pulseTimeInfo["pulse_y"][sampNum]
            plot_x.append(sampTime)
            plot_y.append( sampVal )
            continue #normal mode case
            bin32 = "{0:b}".format(pulseTimeInfo["pulse_y_bin"][sampNum])
            if len(bin32) != 32 : continue
            mdacBits = bin32[4:12]
            #print(sampNum,sampTime,pulseTimeInfo["pulse_y"][sampNum],bin32,mdacBits)
            #break
            if mdacBits not in mdacSamples:
              mdacSamples[mdacBits] = {"x":[],"y":[]}
            #mdacSamples[mdacBits]["x"].append(sampTime)
            #mdacSamples[mdacBits]["y"].append(sampVal)
            #continue
            
            if prevMdacBits == None :
              prevMdacBits = mdacBits
              prevBits = bin32
              continue

            #case 1
            #if sampTime*25 > -46 and sampTime*25 < -45.5 and sampVal > 28585 : sampA.append( (bin32,prevBits) )
            #if sampTime*25 > -46 and sampTime*25 < -45.5 and sampVal < 28585 : sampB.append( (bin32,prevBits) )
            #case 2
            #if sampTime*25 > -45.5 and sampTime*25 < -45.2 : sampA.append( (bin32,prevBits,round(sampTime*25,2)) )
            #if sampTime*25 > -45.2 and sampTime*25 < -45.  : sampB.append( (bin32,prevBits,round(sampTime*25,2)) )
            #case 3
            #if sampTime*25 > -44.7 and sampTime*25 < -44.5 : sampA.append( (bin32,prevBits,round(sampTime*25,2)) )
            #if sampTime*25 > -44.5 and sampTime*25 < -44.3 : sampB.append( (bin32,prevBits,round(sampTime*25,2)) )
            #case 4
            #if sampTime*25 > -21 and sampTime*25 < -20.5 and sampVal > 28585 : sampA.append( (bin32,prevBits,round(sampTime*25,2)) )
            #if sampTime*25 > -21 and sampTime*25 < -20.5 and sampVal < 28585 : sampB.append( (bin32,prevBits,round(sampTime*25,2)) )
            #case 5
            #if sampTime*25 > -20.5 and sampTime*25 < -20.0 and sampVal > 28950 : sampA.append( (bin32,prevBits) )
            #if sampTime*25 > -20.5 and sampTime*25 < -20.0 and sampVal < 28950 : sampB.append( (bin32,prevBits) )
            #case 6
            if sampTime*25 > -20.0 and sampTime*25 < -19.0 and sampVal > 29150 : sampA.append( (bin32,prevBits) )
            if sampTime*25 > -20.0 and sampTime*25 < -19.0 and sampVal < 29150 : sampB.append( (bin32,prevBits) )
              
            if mdacBits != prevMdacBits :
              mdacSamples[prevMdacBits]["x"].append(sampTime)
              mdacSamples[prevMdacBits]["y"].append(sampVal)
            else :
              plot_x_same.append(sampTime)
              plot_y_same.append(sampVal)
            #mdacSamples[prevMdacBits]["x"].append(sampTime)
            #mdacSamples[prevMdacBits]["y"].append(sampVal)
            prevMdacBits = mdacBits
    
    for samp in sampA :
      print( samp )
    print("SAMPB")
    for samp in sampB :
      print( samp )        
            
    mdacCodeDict = { "00000000": ("MDAC SR 0"), "00000001": ("MDAC SR 1"), "00000011": ("MDAC SR 2"), "00000111": ("MDAC SR 3"),\
                     "00001111": ("MDAC SR 4"), "00011111": ("MDAC SR 5"), "00111111": ("MDAC SR 6"), "01111111": ("MDAC SR 7"),\
                     "11111111": ("MDAC SR 8") }
    if True :    
      fig, axes = plt.subplots(1,1,figsize=(13, 8))
      axes.plot([x*25 for x in plot_x],plot_y,".",label="Normal Mode")
      #axes.plot([x*25 for x in plot_x_same],plot_y_same,".",label="No Transition")
      #axes.plot(plot_avg_x,plot_avg_y,"-",label="Average")
      #axes.plot(plot_samp_x,plot_samp_y,".",label="Samples")
      for mdacCode in mdacCodeDict :
        if True and mdacCode in mdacSamples:
          axes.plot( [x*25 for x in mdacSamples[mdacCode]["x"]],mdacSamples[mdacCode]["y"],".",label=mdacCodeDict[mdacCode] )
      axes.set_xlabel('Sample Time [ns]', fontsize=20)
      axes.set_ylabel('Sample Value [ADC]', fontsize=20)
      axes.set_title("Sample Value vs Sample Time", fontsize=20)
      fig.suptitle("Run1813 Channel 5 Pulse Shape", fontsize=20)
      #axes.set_xlim(-46.8,-43.8) #case 1,2,3
      #axes.set_ylim(13000,14300) #case 1,2,3
      #axes.set_xlim(-22.0,-18.0) #case 4,5,6
      #axes.set_ylim(28400,29400) #case 4,5,6
      #axes.set_xlim(-16,-11)
      axes.set_xlim(-170,400)
      #axes.set_ylim(14000,15000)
      #axes.set_xlim(-75,75)
      axes.legend(fontsize=12)
      fig.tight_layout()
      plt.show()
  
    #plotQuick(plot_x,plot_y,label="CV4 Board Pulse")
    return None
      
  
  if True :
    cv4AnalyzeWaveform = CV4_ANALYZE_WAVEFORM("")
    cv4AnalyzeWaveform.runResultsDict = cv4ProcessFile.runResultsDict
    for measNum in cv4AnalyzeWaveform.runResultsDict["results"] :
      #if measNum != "Measurement_891" : continue
      doPlot = True
      meanvals, stdvals, maxval,minval,enob,sinad = cv4AnalyzeWaveform.viewWaveform(chId=chanName,measNum=measNum,doPrint=False,doPlot=doPlot)
    return None
    
  if False :
    cv4AnalyzeWaveform = CV4_ANALYZE_WAVEFORM("")
    cv4AnalyzeWaveform.runResultsDict = cv4ProcessFile.runResultsDict

    plot_x = []
    plot_y = []
    for measNum in cv4AnalyzeWaveform.runResultsDict["results"] :
      awgAmp = float(cv4AnalyzeWaveform.runResultsDict["results"][measNum]["attrs"]['awgAmp'])
      #cv4AnalyzeWaveform.printEnob(chId="channel2",measNum=measNum)
      print( cv4AnalyzeWaveform.runResultsDict["results"][measNum]["attrs"]["SARConstsDdpu"] )
      if False :
      #if True and cv4AnalyzeWaveform.runResultsDict["results"][measNum]["attrs"]["SARConstsDdpu"][0] > 0 :
        #print("\tSAR", [ round(x /self.runResultsDict["results"][measNum]["attrs"]["SARConstsDdpu"][0] *3476.50,2) for x in self.runResultsDict["results"][measNum]["attrs"]["SARConstsDdpu"] ] )
        #print("\tSAR", [ round(x /4.,2) for x in self.runResultsDict["results"][measNum]["attrs"]["SARConstsDdpu"] ] )
        SARConstsDdpu  = cv4AnalyzeWaveform.runResultsDict["results"][measNum]["attrs"]["SARConstsDdpu"]
        MDACConstsDdpu = cv4AnalyzeWaveform.runResultsDict["results"][measNum]["attrs"]["MDACConstsDdpu"]
        convSar  = [ x/4. for x in SARConstsDdpu]
        convMdac = []
        prevVal = 0
        for corrNum,corr in enumerate(MDACConstsDdpu) :
          convMdac.append( corr/4. - prevVal )
          prevVal = corr/4.
        convMdac = convMdac[::-1]
        convMdac[7] = 4096 + convMdac[7]
        cv4AnalyzeWaveform.analyzeSample.mdacWeights = convMdac
        cv4AnalyzeWaveform.analyzeSample.sarWeights = convSar
        print("MDAC",convMdac)
        print("SAR",convSar)
        
      doPlot = False
      #if awgAmp > 2. : doPlot = True
      meanvals, stdvals, maxval,minval,enob = cv4AnalyzeWaveform.viewWaveform(chId=chanName,measNum=measNum,doPrint=True,doPlot=doPlot)
      
      #plot_x.append(awgAmp)
      plot_x.append( (maxval - minval)/16.377 )
      #plot_x.append( (maxval - minval) )
      plot_y.append(enob)
      #print( measNum ,resultsDict[measNum])
      #print("\n")
      continue
    plotQuick(plot_x,plot_y,label="CV4 Board E170685 Ch1, 4.985MHz Sine Wave, 8000 Samples")
    return None
  
  if False: #DAC scan
    print("HERE")
    cv4AnalyzeDacScan = CV4_ANALYZE_DACSCAN("")
    cv4AnalyzeDacScan.runResultsDict = cv4ProcessFile.runResultsDict
    # = [4735.869653510426, 4734.5631849608735+2, 4735.1601059361, 4737.514001857277, 4737.187384719889-2, 4739.642644580255-4, 4740.091395244682-3, 4739.316027442867]
    #sarWeights = [3913.774318704171, 2236.3024890577885, 1119.2178184223656, 700.094823863811, 420.2289877617106, 280.3996861818317, 139.7223063107344, 239.79441840852482, 137.6567344832514, 68.68093902340458, 34.540325411804176, 25.59833657630727, 17.05786156659972, 10.617298234440584, 6.369022916411194, 4.2264482831239265, 2.118157186602517, 1.0942642691259024, 0.5651523904188553, 0.3296175099882698]
    #cv4AnalyzeDacScan.analyzeSample.mdacWeights = [x*1. for x in mdacWeights]
    #cv4AnalyzeDacScan.analyzeSample.sarWeights =  [x*1 for x in sarWeights]
    #cv4AnalyzeDacScan.analyzeSample.mdacWeights = [x/0.8878897240938903 for x in mdacWeights]
    #cv4AnalyzeDacScan.analyzeSample.sarWeights =  [x/0.8878897240938903 for x in sarWeights]
    #cv4AnalyzeDacScan.viewDacScanWaveforms(chId=chanName)
    #cv4AnalyzeDacScan.analyzeSample.applyDdpuCorr = False
    #cv4AnalyzeDacScan.analyzeSample.dropOverFlowSamples = True
    #cv4AnalyzeDacScan.analyzeSample.sampScaleFactor = 1.
    #cv4AnalyzeDacScan.analyzeSample.sampOffsetVal = 7275.0-41
    cv4AnalyzeDacScan.plotDacLinearityData(chId=chanName,plotTitle="Board E170685 Run 1764 Channel 1")
    #cv4AnalyzeDacScan.viewDacScanDist(chId=chanName)
    #cv4AnalyzeDacScan.printMeasSamples(chId=chanName)
    #cv4AnalyzeDacScan.plotMdacBits(chId=chanName,plotTitle="Run 911: Board 153, Ch3 Connected, MDAC RESETEN OFF")
    
  #simpleTest(cv4ProcessFile.runResultsDict)
  
  return

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
