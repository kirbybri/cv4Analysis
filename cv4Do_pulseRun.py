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

import statsmodels.api as sm
import scipy.stats
import json

from scipy.optimize import curve_fit

from cv4AnalyzeSample import CV4_ANALYZE_SAMPLE

def simpleTest(resultsDict):
  cv4AnalyzeWaveform = CV4_ANALYZE_WAVEFORM("")
  cv4AnalyzeWaveform.runResultsDict = resultsDict
  for measNum in resultsDict["results"] :
    if measNum != "Measurement_6000" : continue
    #cv4AnalyzeWaveform.printEnob(chId="channel2",measNum=measNum)
    cv4AnalyzeWaveform.viewWaveform(chId="channel6",measNum=measNum,doPrint=True,doPlot=True)
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

def measureLinearity(xs,ys,ysErr,lowLim,upLim):
    if len(xs) < 3 or len(ys) < 3 :
      print("ERROR TOO FEW POINTS")
      return None
    if len(xs) != len(ys) :
      print("ERROR MISMATCHED LENGTHS")
      return None
    xsFit = []
    ysFit = []
    for num in range(0,len(xs),1):
      if xs[num] <= lowLim or xs[num] > upLim :
        continue
      xsFit.append(xs[num])
      ysFit.append(ys[num])
    if len(ysFit ) < 3 :
      print("ERROR TOO FEW POINTS")
      return None   

    xsFit = sm.add_constant(xsFit)
    #model = sm.OLS(ysFit,xsFit)
    model = sm.GLM(ysFit,xsFit)
    results = model.fit()
    if len(results.params) < 2 :
      print("ERROR FIT FAILED")
      return None
    slope = results.params[1]
    intercept = results.params[0]
    slopeErr = results.bse[1]
    interceptErr = results.bse[0]

    #calculate reduced chi-sq
    chiSq = 0
    resid_x = []
    resid_y = []
    resid_yRms = []
    for num in range(0,len(xs),1):
      if xs[num] <= lowLim or xs[num] > upLim :
        continue
      predY = xs[num]*slope + intercept
      resid = ys[num] - predY
      if ysErr[num] > 0 :
        chiSq = chiSq + resid*resid/ysErr[num]/ysErr[num]
      resid_x.append(xs[num])
      resid_y.append(resid)
      #resid_y.append(resid / ys[num] * 100. )
      resid_yRms.append(ysErr[num])
    chiSq = chiSq / float( len(ysFit ) - 2 )

    print( "SLOPE ", slope , "\tERR ", slopeErr,"\tCHI2 ",chiSq )
    print( results.summary() )
    return slope, intercept, slopeErr, interceptErr, chiSq,resid_x,resid_y,resid_yRms

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

  if True :
    simpleTest(cv4ProcessFile.runResultsDict)
    
  if True :
    cv4AnalyzeWaveform = CV4_ANALYZE_WAVEFORM("")
    cv4AnalyzeWaveform.runResultsDict = cv4ProcessFile.runResultsDict
    cv4AnalyzePulse = CV4_ANALYZE_PULSES()
    
    measResults = {}
    for measNum in cv4AnalyzeWaveform.runResultsDict["results"] :
      vals = cv4AnalyzeWaveform.getMeasChData(chId=chanName,measNum=measNum)
      vals32bit = cv4AnalyzeWaveform.getMeasChData(chId=chanName,measNum=measNum,get32Bit=True)
      fitPulseSetDict = cv4AnalyzePulse.findPulses(vals)
      if fitPulseSetDict == None :
        print("NO PULSES FOUND!")
        continue
      if len(fitPulseSetDict) == 0 :
        print("NO PULSES FOUND!")
        continue
      result = cv4AnalyzePulse.getPulses(vals,fitPulseSetDict,chWfBin=vals32bit)
      if result == None : 
        print("NO PULSES SAVED")
        continue
      measResults[measNum] = result
    
    plot_x = []
    plot_y = []
    ampSamples = {}
    for measNum in measResults :
      #if measNum != "Measurement_146" : continue
      measAttrs = cv4AnalyzeWaveform.runResultsDict["results"][measNum]["attrs"]
      awgAmp = 0.0
      if 'awgAmp' in measAttrs :
        awgAmp = measAttrs['awgAmp']
        print(measNum,"\t",awgAmp)
        awgAmp = float(awgAmp)/512.*1000. #convert to mA
        awgAmp = round(float(awgAmp),2)
      print("AMP VAL", awgAmp )
      attVal = 0.0
      if 'attVal' in measAttrs :
        attVal = measAttrs['attVal']
      #attVal = 0.0
      print("ATT VAL", attVal)
      if attVal > 0 :
        #convert AWG amp to signal current here
        awgAmp = awgAmp/(  pow(10,attVal/20.)   )
        awgAmp = round(float(awgAmp),5)
      print("SIGNAL VAL", awgAmp )
      if awgAmp < 5. : continue
      #if awgAmp < 0.3 or awgAmp > 11 : continue #LG
      #if awgAmp < 0.025 or awgAmp > 0.5 : continue #HG
      #if awgAmp < 0.01 or awgAmp > 0.5 : continue #HG
      #if awgAmp < 9.6 or awgAmp > 10 : continue
      if awgAmp not in ampSamples:
        ampSamples[awgAmp] = {"x":[],"y":[],"count":0}
      ampSamples[awgAmp]["count"] = ampSamples[awgAmp]["count"] + 1
      if ampSamples[awgAmp]["count"] > 5 : continue
      fitPulseSetDict = measResults[measNum]
      for pulseSet in fitPulseSetDict :
        #print("PULSE SET ", pulseSet)
        pulseSetInfo = fitPulseSetDict[pulseSet]
        for pulseNum in pulseSetInfo :
          #if pulseNum < 4 : continue
          #if pulseNum > 29 + 4 : continue
          pulseTimeInfo = pulseSetInfo[pulseNum]
          for sampNum,sampTime in enumerate(pulseTimeInfo["pulse_fit_x"]):
            sampVal = pulseTimeInfo["pulse_y"][sampNum]
            plot_x.append(sampTime)
            plot_y.append( sampVal )
            ampSamples[awgAmp]["x"].append(sampTime)
            ampSamples[awgAmp]["y"].append(sampVal )
    
    #measure pedestal and pulse height
    for awgAmp in ampSamples :
      result = cv4AnalyzePulse.getAvgPulses(samp_x=ampSamples[awgAmp]["x"],samp_y=ampSamples[awgAmp]["y"],pedTime=-5.)
      if result == None : continue
      binsDict, pedVal, pedRms, avgPulseHeight, riseTime, avg_pulse_x, avg_pulse_y = result
      if avgPulseHeight > 25000 : continue
      ampSamples[awgAmp]["avg_x"] = avg_pulse_x
      ampSamples[awgAmp]["avg_y"] = avg_pulse_y
      ampSamples[awgAmp]["pulseHeight"] = avgPulseHeight
      ampSamples[awgAmp]["pedVal"] = pedVal
      ampSamples[awgAmp]["pedRms"] = pedRms
      ampSamples[awgAmp]["riseTime"] = riseTime
      print( awgAmp ,"\t", avgPulseHeight )
      
    #plot shapes   
    if True :    
      fig, axes = plt.subplots(1,1,figsize=(13, 8))
      #axes.plot([x*25 for x in plot_x],plot_y,".",label="Normal Mode")
      for awgAmp in ampSamples :
        if "pulseHeight" in ampSamples[awgAmp] :
          #axes.plot( [x*25 for x in ampSamples[awgAmp]["x"]],ampSamples[awgAmp]["y"],".",label=str(awgAmp) +" mA" )
          axes.plot( [x*25 for x in ampSamples[awgAmp]["avg_x"]],ampSamples[awgAmp]["avg_y"],".-",label="AVG "+str(awgAmp) +" mA" )
      axes.set_xlabel('Sample Time [ns]', fontsize=20)
      axes.set_ylabel('Sample Value [ADC]', fontsize=20)
      axes.set_title("Sample Value vs Sample Time", fontsize=20)
      #fig.suptitle("Run1813 Channel 5 Pulse Shape", fontsize=20)
      #axes.set_xlim(-46.8,-43.8) #case 1,2,3
      #axes.set_ylim(13000,14300) #case 1,2,3
      #axes.set_xlim(-22.0,-18.0) #case 4,5,6
      #axes.set_ylim(28400,29400) #case 4,5,6
      #axes.set_xlim(-16,-11)
      axes.set_xlim(-170,400)
      #axes.set_ylim(14000,15000)
      #axes.set_xlim(-75,75)
      axes.legend(loc=1,ncol=2,fontsize=12)
      fig.tight_layout()
      plt.show()
    
    #plot lineartiy
    if True :
      amps = []
      heights = []
      rms = []
      for amp in ampSamples :
        if "pulseHeight" not in ampSamples[amp]: continue
        amps.append(amp)
        heights.append( ampSamples[amp]["pulseHeight"] )
        rms.append( ampSamples[amp]["pedRms"] )
      #lineResult = measureLinearity(amps,heights,rms,0,10)
      lineResult = None
      if len(amps) > 3:
        lineResult = measureLinearity(amps,heights,rms,np.min(amps)-0.1,np.max(amps)+0.1)
      X_plotFit = []
      Y_plotFit = []
      resid_x = []
      resid_y = []
      resid_yRms = []
      textstr = ""
      #slope, intercept, slopeErr, interceptErr, chiSq,resid_x,resid_y,resid_yRms = result
      if lineResult != None :
        #slope, intercept, slopeErr, interceptErr, chiSq,resid_x,resid_y,resid_yRms = linResult
        slope = lineResult[0]
        intercept = lineResult[1]
        slopeErr = lineResult[2]
        interceptErr = lineResult[3]
        chiSq = lineResult[4]
        resid_x = lineResult[5]
        resid_y = lineResult[6]
        resid_yRms = lineResult[7]
        X_plotFit = np.linspace(np.min(amps)-0.1,np.max(amps)+0.1,1000)
        Y_plotFit = X_plotFit*slope + intercept
        textstr = '\n'.join((
          r'$m=%.3f\pm%.3f$' % (slope, slopeErr, ),
          r'$b=%.2f\pm%.2f$' % (intercept,interceptErr, ),
          r'$\chi^2=%.2f$' % (chiSq, )
        ))
      
      fig, axes = plt.subplots(1,1,figsize=(13, 8))
      axes.plot( amps , heights ,".",label="Linearity" )
      axes.plot(X_plotFit,Y_plotFit,"-")
      axes.set_xlabel('Signal Amplitude [mA]', fontsize=20)
      axes.set_ylabel('Pulse Height [ADC]', fontsize=20)
      axes.set_title("Pulse Height vs Signal Amplitude: "+str(fileName), fontsize=20)
      axes.text(0.05, 0.95, textstr, transform=axes.transAxes, fontsize=14, verticalalignment='top', bbox=dict(facecolor='white', alpha=1.0, edgecolor='black'))
      axes.legend(fontsize=12)
      fig.tight_layout()
      plt.show()
      
      fig, axes = plt.subplots(1,1,figsize=(13, 8))
      #axes.plot( amps , resid_y / np.max(heights)*100 ,".",label="Non-linearity" )
      axes.plot( amps , resid_y  ,".",label="Non-linearity" )
      #axes.plot(X_plotFit,Y_plotFit,"-")
      axes.set_xlabel('Signal Amplitude [mA]', fontsize=20)
      #axes.set_ylabel('Pulse Height Residual [%]', fontsize=20)
      axes.set_ylabel('Pulse Height Residual [ADC]', fontsize=20)
      axes.set_title("Pulse Height Residual vs Signal Amplitude: "+str(fileName), fontsize=20)
      #axes.text(0.05, 0.95, textstr, transform=axes.transAxes, fontsize=14, verticalalignment='top', bbox=dict(facecolor='white', alpha=1.0, edgecolor='black'))
      axes.legend(fontsize=12)
      fig.tight_layout()
      plt.show()
      
    #plot peaking time
    if True :
      fig, axes = plt.subplots(1,1,figsize=(13, 8))
      axes.plot( [x for x in ampSamples],[ ampSamples[x]["riseTime"] for x in ampSamples],".",label="Rise time" )
      axes.set_xlabel('Signal Amplitude [mA]', fontsize=20)
      axes.set_ylabel('Rise Time [ns]', fontsize=20)
      axes.set_title("Pulse Rise Time vs Signal Amplitude", fontsize=20)
      axes.legend(fontsize=12)
      axes.grid()
      fig.tight_layout()
      plt.show()
  
    #plotQuick(plot_x,plot_y,label="CV4 Board Pulse")
    return None
      

  return

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
