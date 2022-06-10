import json
import sys
import subprocess
import numpy as np
from math import *
#import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import statsmodels.api as sm
import scipy.stats

from templateFit import TEMPLATE_FIT

def measureLinearity(xs,ys):
    if len(xs) < 3 or len(ys) < 3 :
      print("MEASURE LINEARITY ERROR TOO FEW POINTS")
      return None
    if len(xs) != len(ys) :
      print("MEASURE LINEARITY ERROR MISMATCHED LENGTHS")
      return None
    xsFit = xs
    ysFit = ys

    xsFit = sm.add_constant(xsFit)
    #model = sm.OLS(ysFit,xsFit)
    model = sm.GLM(ysFit,xsFit)
    results = model.fit()
    if len(results.params) < 2 :
      print("MEASURE LINEARITY ERROR FIT FAILED")
      return None
    slope = results.params[1]
    intercept = results.params[0]
    slopeErr = results.bse[1]
    interceptErr = results.bse[0]
    
    #print( "MEASURE LINEARITY  SLOPE ", slope , "\tERR ", slopeErr )
    #print( results.summary() )
    return slope, intercept, slopeErr, interceptErr

#BEGIN SLICE_ANALYZE_WAVEFORM CLASS
class CV4_ANALYZE_PULSES(object):

  #__INIT__#
  def __init__(self, fileName = None):
    self.fileName = fileName
    self.slopes = []
    self.intercepts = []
    return

  def findPulses(self,chWf):
    #print( len(chWf),chWf )
    #find maximum sample
    maxSamp = np.max(chWf)
    meanSamp = np.mean(chWf)
    threshold = meanSamp + (maxSamp-meanSamp)*0.5
    threshold_fall = meanSamp + (maxSamp-meanSamp)*0.3
    threshold_2 = meanSamp + (maxSamp-meanSamp)*0.7
    threshold_fall_2 = meanSamp + (maxSamp-meanSamp)*0.8
    #threshold = 16400      #fixed threshold
    #threshold_fall = 12700 #fixed threshold
    #print(maxSamp,meanSamp,threshold,threshold_fall)
    setDiffThreshold = 1300
  
    prevSamp = chWf[0]
    pulseTimeList = []
    pulse_x = []
    pulse_y = []
    pulse_x_lo = []
    pulse_y_lo = []
  
    #get pulse times
    doFallTime = True
    for idx, samp in enumerate(chWf):
      #if idx > 50000 : break
      if idx <= 10 : continue
      if idx > len(chWf) - 25 : continue #require some offset from waveform ends
      if samp > threshold and prevSamp <= threshold :
        #have a pulse rising edge 
        sampNum_rise_low = idx-1
      
        #find falling edge
        sampNum_fall_low = -1
        for sampNum in range(idx,idx+10,1):
          if chWf[sampNum-1] > threshold_fall and chWf[sampNum] <= threshold_fall :
            sampNum_fall_low = sampNum-1
            break
        #if sampNum_fall_low < 0 :
        #  print("ERROR could not find falling edge, should not happen, exiting")
        #  return None

        #rising edge pulse time
        #print(idx, (idx-1) + (threshold - prevSamp)/float(samp-prevSamp) , (idx) - (samp - threshold)/float(samp-prevSamp) )
        pulseTime_rise = sampNum_rise_low + (threshold - float(chWf[sampNum_rise_low]) )/float( float(chWf[sampNum_rise_low+1])-float(chWf[sampNum_rise_low]) )
        pulseTime_fall = 0
        if doFallTime and chWf[sampNum_fall_low+1] !=  chWf[sampNum_fall_low] :
          pulseTime_fall = sampNum_fall_low + (threshold_fall - float(chWf[sampNum_fall_low]) )/float( float(chWf[sampNum_fall_low+1])-float(chWf[sampNum_fall_low]) )
      
        pulseTime = pulseTime_rise
        if doFallTime :      
          pulseTime = (pulseTime_rise+pulseTime_fall)/2.
        pulseTimeList.append( (idx,pulseTime) )
        #if len(pulseTimeList) > 30 : break
      prevSamp = samp
    if len(pulseTimeList) == 0 :
      print("NO PULSES FOUND")
      return None
    pulseTimeList.pop(0)
  
    if False :
      fig, axes = plt.subplots(1,1,figsize=(13, 8))
      axes.plot(chWf,".")
      axes.set_xlabel('Sample #', fontsize=20)
      axes.set_ylabel('Sample Value', fontsize=20)
      axes.set_title("Sample Vs Sample #", fontsize=20)
      #fig.suptitle(label, fontsize=16)
      #axes.set_xlim(0,10000)
      #axes.set_ylim(6700,6850)
      fig.tight_layout()
      plt.show()    
  
    if len(pulseTimeList) == 0 :
      print("NO PULSES FOUND")
      return None

    if len(pulseTimeList) < 3 :
      print("< 3 PULSES FOUND")
      return None
  
    #divide pulses into groups
    setCount = 0
    prevPulseIdx = pulseTimeList[0][0]
    pulseSetDict = {}
    for pulseTimeInfo in pulseTimeList:
      idx = pulseTimeInfo[0]
      pulseTime = pulseTimeInfo[1]
      #print( idx - prevPulseIdx )
      if idx > prevPulseIdx + setDiffThreshold :
        setCount += 1
      if setCount not in pulseSetDict :
        pulseSetDict[setCount] = []
      pulseSetDict[setCount].append(pulseTimeInfo)
      prevPulseIdx = idx
    #print(setCount,len(pulseSetDict),pulseSetDict.keys())
  
    #for each group, perform linear fit
    fitPulseSetDict = {}
    for pulseSet in pulseSetDict :
      #print("PULSE SET ", pulseSet)
      pulseSetInfo = pulseSetDict[pulseSet]
      xs = []
      ys = []
      for pulseNum, pulseTimeInfo in enumerate(pulseSetInfo) :
        idx = pulseTimeInfo[0]
        pulseTime = pulseTimeInfo[1]
        xs.append(pulseNum)
        ys.append(pulseTime)
      fitResult = measureLinearity(xs,ys)
      if fitResult == None : continue
      slope, intercept, slopeErr, interceptErr = fitResult
      #if intercept < 23.83 or intercept > 23.8445: continue
      #if intercept > 216.85 :continue
      self.slopes.append(slope)
      self.intercepts.append(intercept)
      #slope = 191.97032546677994
      #intercept = 215.54061395436648
      #print("INTERCEPT",intercept)
      #print("SLOPE",slope)
      resid_x = []
      resid_y = []
      for pulseNum, pulseTimeInfo in enumerate(pulseSetInfo) :
        idx = pulseTimeInfo[0]
        pulseTime = pulseTimeInfo[1]
        #print( pulseTime , slope*pulseNum + intercept )
        if setCount not in fitPulseSetDict :
          fitPulseSetDict[setCount] = {}
        fitPulseSetDict[setCount][pulseNum] = {"idx":idx,"pulseTime":pulseTime,"fitPulseTime":slope*pulseNum + intercept }
        resid_x.append(pulseNum)
        resid_y.append(pulseTime - (slope*pulseNum + intercept) )
      if False :
        fig, axes = plt.subplots(1,1,figsize=(13, 8))
        axes.plot(resid_x,resid_y,".")
        axes.set_xlabel('Pulse #', fontsize=20)
        axes.set_ylabel('Residual Time', fontsize=20)
        #axes.set_title("Pulse Sample Vs Pulse #", fontsize=20)
        #fig.suptitle(label, fontsize=16)
        fig.tight_layout()
        plt.show()
    #print("Number groups",setCount,len(fitPulseSetDict),fitPulseSetDict.keys())
  
    if False :
      pulseIdxList = [ x[0] for x in pulseTimeList]
      fig, axes = plt.subplots(1,1,figsize=(13, 8))
      axes.plot(pulseIdxList,".")
      axes.set_xlabel('Pulse #', fontsize=20)
      axes.set_ylabel('Pulse Sample', fontsize=20)
      axes.set_title("Pulse Sample Vs Pulse #", fontsize=20)
      #fig.suptitle(label, fontsize=16)
      fig.tight_layout()
      plt.show()
    return fitPulseSetDict
    
  def getPulses(self,chWf,fitPulseSetDict,chWfBin=[]):
    if len(chWf) == 0 :
      print("NO SAMPLES")
      return None
    #don't break code
    if len(chWfBin)==0 :
      chWfBin = chWf #just use sample array
    #define pulse shapes
    pulse_x = []
    pulse_y = []
    pulse_y_bin = []
    pulse_ind = []
    pedVals = []
    for pulseSet in fitPulseSetDict :
      #print("PULSE SET ", pulseSet)
      pulseSetInfo = fitPulseSetDict[pulseSet]
      #for pulseNum, pulseTimeInfo in enumerate(pulseSetInfo) :
      for pulseNum in pulseSetInfo :
        pulseTimeInfo = pulseSetInfo[pulseNum]
        idx = pulseTimeInfo["idx"]
        pulseTime = pulseTimeInfo["pulseTime"]
        fitPulseTime = pulseTimeInfo["fitPulseTime"]
        pulseTimeInfo["pulse_x"] = []
        pulseTimeInfo["pulse_y"] = []
        pulseTimeInfo["pulse_y_bin"] = []
        pulseTimeInfo["pulse_fit_x"] = []
        pulseTimeInfo["pulse_ind_x"] = []
        for sampNum in range(-10,40,1):
          if idx+sampNum < 0 or idx+sampNum >= len(chWf) : continue
          pulse_x.append(sampNum + idx - pulseTime)
          pulse_y.append(chWf[idx+sampNum])
          pulse_y_bin.append(chWfBin[idx+sampNum])
          pulse_ind.append(idx+sampNum)
          pulseTimeInfo["pulse_x"].append(sampNum)
          pulseTimeInfo["pulse_ind_x"].append(sampNum + idx)
          pulseTimeInfo["pulse_y"].append(chWf[idx+sampNum])
          pulseTimeInfo["pulse_y_bin"].append(chWfBin[idx+sampNum])
          pulseTimeInfo["pulse_fit_x"].append(sampNum + idx - pulseTime)
          if sampNum + idx - pulseTime < -5 :
            pedVals.append( chWf[idx+sampNum]) 
        #print( pulseTimeInfo )
    return fitPulseSetDict

  def getAvgPulses(self,samp_x=[],samp_y=[],pedTime=0.):
    if len(samp_x) < 3 :
      return None
    if len(samp_x) != len(samp_y) :
      return None
    
    if False :
      fig, axes = plt.subplots(1,1,figsize=(13, 8))
      axes.plot(samp_x,samp_y,".",label="Samples")
      axes.set_xlabel('Sample Number [#]', fontsize=20)
      axes.set_ylabel('Sample Value [ADC]', fontsize=20)
      axes.set_title("Sample Value vs Sample Number", fontsize=20)
      axes.legend(fontsize=12)
      fig.tight_layout()
      plt.show()

    #get pedestal
    pedVals = []
    for sampNum,sampTime in enumerate(samp_x) :
      if sampTime < pedTime :
        pedVals.append( samp_y[sampNum] )
    if len(pedVals) < 3:
      return None
    pedVal = round(np.mean(pedVals),2)
    pedRms = round(np.std(pedVals),2)
    pedErr = round(pedRms/np.sqrt(len(pedVals)),2)

    #get avg pulse
    #bins = np.linspace(-10,25,301) #correct
    bins = np.linspace(-10,25,701)
    #print( bins)

    binsDict = {}
    for binNum,binEdge in enumerate(bins):
      if binNum == 0 : continue
      binsDict[binNum] = {"low":bins[binNum-1],"high":bins[binNum],"samples":[]}
      #print( binsDict[binNum] )
  
    avg_pulse_x = []
    avg_pulse_y = []
    for binNum in binsDict :
      low = binsDict[binNum]["low"]
      high = binsDict[binNum]["high"]
      for sampNum, samp in enumerate(samp_x):
        if samp_x[sampNum] >= low and samp_x[sampNum] < high :
          binsDict[binNum]["samples"].append( samp_y[sampNum] )
      #if len(binsDict[binNum]["samples"]) == 0 : continue # no samples in bin, skip
      if len(binsDict[binNum]["samples"]) == 0 :
        binsDict[binNum]["samples"].append(0)
      binsDict[binNum]["center"] = (high+low)/2.
      binsDict[binNum]["sampMean"] = np.mean( binsDict[binNum]["samples"] )
      #print( binsDict[binNum]["center"] , binsDict[binNum]["sampMean"] )   
      avg_pulse_x.append(binsDict[binNum]["center"])
      avg_pulse_y.append(binsDict[binNum]["sampMean"])

    avgPulseHeight = np.max(avg_pulse_y) - pedVal
    #print( avgPulseHeight )

    #rise time measurement
    perRisePrev = 100
    sampTimePrev = 0
    sampTime_rise5Per = None
    sampTime_rise100Per = None
    for sampNum, sampVal in enumerate(avg_pulse_y) :
      sampTime = avg_pulse_x[sampNum]*25
      pedCorrSampVal = sampVal - pedVal
      perRise = pedCorrSampVal/avgPulseHeight*100
      if perRise > 5 and perRisePrev < 5 :
        #find 5% sample
        #sampTime_rise5Per = sampTime
        slope = (perRise - perRisePrev)/(sampTime - sampTimePrev)
        sampTime_rise5Per = sampTimePrev + (5 - perRisePrev )/slope
      if perRise > 99.999 :
        sampTime_rise100Per = sampTime
      #print( sampNum,"\t",sampTime,"\t",pedCorrSampVal,"\t",pedCorrSampVal/avgPulseHeight*100)
      perRisePrev = perRise
      sampTimePrev = sampTime
    riseTime = 0.
    if sampTime_rise5Per != None and sampTime_rise100Per != None :
      riseTime = sampTime_rise100Per - sampTime_rise5Per
    

    #simple function for avg pulse display
    if False :
      print("Average pulse samples")
      #print_x = [round(x,2) for x in avg_pulse_x]
      #print_y = [round(x,2) for x in avg_pulse_y]
      #if avgPulseHeight == 0 : avgPulseHeight = 1 # shouldn't happen
      #print_y = [round( (x - pedVal)/avgPulseHeight ,2) for x in avg_pulse_y]
    
      plot_x = [x*25 for x in samp_x]
      plot_y = [x for x in samp_y]
      plot_avg_x = [x*25 for x in avg_pulse_x]
      plot_avg_y = [x for x in avg_pulse_y]
      #print("X values")
      #print(print_x)
      #print("Y values")
      #print(print_y)
      print( "PED",pedVal,"RMS",pedRms,"Pulse Height",avgPulseHeight)
      print("RISETIME ", riseTime)
    
      #for sampNum, samp in enumerate(plot_y) :
      #  print( plot_y[sampNum],"\t",plot_y_bin[sampNum])
  
      fig, axes = plt.subplots(1,1,figsize=(13, 8))
      #axes.plot(plot_x,plot_y,".",label="Sample Data")
      #axes.plot(plot_avg_x,plot_avg_y,"-",label="Average")
      axes.plot(samp_x,samp_y,".",label="Sample Data")
      axes.plot(avg_pulse_x,avg_pulse_y,"-",label="Average")
      #axes.set_xlabel('Sample #', fontsize=20)
      axes.set_xlabel('Sample Time [ns]', fontsize=20)
      axes.set_ylabel('Sample Value [ADC]', fontsize=20)
      axes.set_title("Sample Value vs Sample Time", fontsize=20)
      #axes.set_xlim(40-4,44+4)
      #axes.set_xlim(-25,10)
      #axes.set_ylim(18000,21000)
      #fig.suptitle(label, fontsize=16)
      axes.legend(fontsize=12)
      fig.tight_layout()
      plt.show()

    return binsDict, pedVal, pedRms, avgPulseHeight, riseTime, avg_pulse_x, avg_pulse_y
    
  def plotPulses_32BitInfo(self,fitPulseSetDict,binsDict,chWfBin=[],plotTitle=""):

    allPulse_fit_x = []
    allPulse_y = []
    mdacSamples = {}
    for pulseSet in fitPulseSetDict :
      #print("PULSE SET ", pulseSet)
      pulseSetInfo = fitPulseSetDict[pulseSet]
      for pulseNum, pulseTimeInfo in enumerate(pulseSetInfo) :
        #print(pulseTimeInfo.keys())
        pulse_x = pulseTimeInfo["pulse_x"]
        pulse_ind_x = pulseTimeInfo["pulse_ind_x"]
        pulse_y = pulseTimeInfo["pulse_y"]
        pulse_fit_x = pulseTimeInfo["pulse_fit_x"]
        for sampNum, samp in enumerate(pulse_x):
          allPulse_fit_x.append(pulse_fit_x[sampNum])
          allPulse_y.append(pulse_y[sampNum])
          samp32bit = chWfBin[ pulse_ind_x[sampNum] ]
          #access 32-bit info for this sample
          #print( samp32bit ,"\t", samp32bit[20:28], "\t", pulse_y[sampNum] )
          mdacBits = samp32bit[20:28] #get MDAC info
          if mdacBits not in mdacSamples:
            mdacSamples[mdacBits] = {"x":[],"y":[]}
          mdacSamples[mdacBits]["x"].append(pulse_fit_x[sampNum])
          mdacSamples[mdacBits]["y"].append(pulse_y[sampNum])
    #end loop

    avg_pulse_x = []
    avg_pulse_y = []
    for binNum in binsDict :
      avg_pulse_x.append(binsDict[binNum]["center"])
      avg_pulse_y.append(binsDict[binNum]["sampMean"])

    #note: bits here 0 to the left, MSB to the right, confusing!
    mdacCodeDict = { "00000000": ("MDAC SR 0"), "10000000": ("MDAC SR 1"), "11000000": ("MDAC SR 2"), "11100000": ("MDAC SR 3"),\
                     "11110000": ("MDAC SR 4"), "11111000": ("MDAC SR 5"), "11111100": ("MDAC SR 6"), "11111110": ("MDAC SR 7"),\
                     "11111111": ("MDAC SR 8") }

    #simple function for avg pulse display
    if True :    
      plot_avg_x = [x*25 for x in avg_pulse_x]
      plot_avg_y = [x for x in avg_pulse_y]
      plot_samp_x = [x*25 for x in allPulse_fit_x]
      plot_samp_y = [x for x in allPulse_y]
      fig, axes = plt.subplots(1,1,figsize=(13, 8))
      axes.plot(plot_avg_x,plot_avg_y,"-",label="Average")
      #axes.plot(plot_samp_x,plot_samp_y,".",label="Samples")
      for mdacCode in mdacCodeDict :
        if True and mdacCode in mdacSamples:
          axes.plot( [x*25 for x in mdacSamples[mdacCode]["x"]],mdacSamples[mdacCode]["y"],".",label=mdacCodeDict[mdacCode] )
      axes.set_xlabel('Sample Time [ns]', fontsize=20)
      axes.set_ylabel('Sample Value [ADC]', fontsize=20)
      axes.set_title("Sample Value vs Sample Time", fontsize=20)
      fig.suptitle(plotTitle, fontsize=20)
      axes.set_xlim(-75,75)
      #axes.set_ylim(22000,30000)
      #axes.set_xlim(-50,-30)
      axes.legend(fontsize=12)
      fig.tight_layout()
      plt.show()
  
    return None
    
  def plotPulses(self,fitPulseSetDict,binsDict,chWfBin=[],plotTitle=""):

    allPulse_fit_x = []
    allPulse_y = []
    for pulseSet in fitPulseSetDict :
      #print("PULSE SET ", pulseSet)
      pulseSetInfo = fitPulseSetDict[pulseSet]
      for pulseNum, pulseTimeInfo in enumerate(pulseSetInfo) :
        pulse_x = pulseTimeInfo["pulse_x"]
        pulse_y = pulseTimeInfo["pulse_y"]
        pulse_fit_x = pulseTimeInfo["pulse_fit_x"]
        for sampNum, samp in enumerate(pulse_x):
          allPulse_fit_x.append(pulse_fit_x[sampNum])
          allPulse_y.append(pulse_y[sampNum])
    #end loop

    avg_pulse_x = []
    avg_pulse_y = []
    for binNum in binsDict :
      if len( binsDict[binNum]["samples"] ) == 0 : continue
      if binsDict[binNum]["sampMean"] == 0 : continue
      avg_pulse_x.append(binsDict[binNum]["center"])
      avg_pulse_y.append(binsDict[binNum]["sampMean"])
    
    #simple function for avg pulse display
    if True :    
      plot_avg_x = [x*25 for x in avg_pulse_x]
      plot_avg_y = [x for x in avg_pulse_y]
      plot_samp_x = [x*25 for x in allPulse_fit_x]
      plot_samp_y = [x for x in allPulse_y]
      fig, axes = plt.subplots(1,1,figsize=(13, 8))
      axes.plot(plot_avg_x,plot_avg_y,".",label="Average")
      axes.plot(plot_samp_x,plot_samp_y,"o",label="Samples")
      axes.set_xlabel('Sample Time [ns]', fontsize=20)
      axes.set_ylabel('Sample Value [ADC]', fontsize=20)
      axes.set_title("Sample Value vs Sample Time", fontsize=20)
      fig.suptitle(plotTitle, fontsize=20)
      #fig.suptitle(label, fontsize=16)
      #axes.set_xlim(-75,475)
      #axes.set_xlim(-75,75)
      #axes.set_xlim(20,60)
      #axes.set_xlim(-25,15)
      #axes.set_ylim(22000,30000)
      #axes.set_ylim(np.max( allPulse_y ) - 2000,np.max( allPulse_y ) + 2000)
      axes.legend(fontsize=12)
      fig.tight_layout()
      plt.show()
  
    return None    

  def getAvgVal(self,avg_pulse_x,avg_pulse_y,time ):
    lowNum = None
    for avgSampNum in range(1,len(avg_pulse_x),1) :
      if avg_pulse_x[avgSampNum-1] <= time and avg_pulse_x[avgSampNum] > time :
        lowNum = avgSampNum-1
        break
    if lowNum == None :
      return None
    lowTime = avg_pulse_x[lowNum]
    highTime = avg_pulse_x[lowNum+1]
    if highTime <= lowTime :
      return None
    lowVal = avg_pulse_y[lowNum]
    highVal = avg_pulse_y[lowNum+1]
    avgVal = lowVal + (highVal - lowVal)*( time - lowTime )/(highTime - lowTime )
    return avgVal

  def plotPulseResiduals(self,fitPulseSetDict,binsDict,chWfBin=[],plotTitle=""):
    avg_pulse_x = []
    avg_pulse_y = []
    for binNum in binsDict :
      avg_pulse_x.append(binsDict[binNum]["center"])
      avg_pulse_y.append(binsDict[binNum]["sampMean"])

    allPulse_fit_x = []
    allPulse_y = []
    allDiff_x = []
    allDiff_y = []
    for pulseSet in fitPulseSetDict :
      pulseSetInfo = fitPulseSetDict[pulseSet]
      for pulseNum, pulseTimeInfo in enumerate(pulseSetInfo) :
        pulse_x = pulseTimeInfo["pulse_x"]
        pulse_y = pulseTimeInfo["pulse_y"]
        pulse_fit_x = pulseTimeInfo["pulse_fit_x"]
        for sampNum, samp in enumerate(pulse_x):
          allPulse_fit_x.append(pulse_fit_x[sampNum])
          allPulse_y.append(pulse_y[sampNum])
          avgVal = self.getAvgVal(avg_pulse_x,avg_pulse_y,pulse_fit_x[sampNum])
          if avgVal == None : continue
          allDiff_x.append(pulse_fit_x[sampNum])
          allDiff_y.append( pulse_y[sampNum] - avgVal )
    if len(allDiff_y) == 0 :
      return None
    maxResid = round(np.max(allDiff_y),2)
    print("MAX RESIDUAL",maxResid)

    #simple function for avg pulse display
    if False :    
      fig, axes = plt.subplots(1,2,figsize=(13, 8))
      axes[0].plot( [x*25 for x in avg_pulse_x   ] ,avg_pulse_y,"-",label="Average")
      axes[0].plot( [x*25 for x in allPulse_fit_x] ,allPulse_y,".",label="Samples")
      axes[0].set_xlabel('Sample Time [ns]', fontsize=20)
      axes[0].set_ylabel('Sample Value [ADC]', fontsize=20)
      axes[0].set_title("Sample Value vs Sample Time", fontsize=20)
      #axes[0].set_title(plotTitle, fontsize=20)
      axes[0].set_xlim(-75,75)
      axes[0].legend(fontsize=12)
      axes[1].plot( [x*25 for x in allDiff_x]  , allDiff_y ,".",label="Residual")
      axes[1].set_xlabel('Sample Time [ns]', fontsize=20)
      axes[1].set_ylabel('Residual [ADC]', fontsize=20)
      axes[1].set_title("Sample Value - Average Value vs Sample Time", fontsize=20)
      textstr = '\n'.join((
        r'$\mu=%.2f$' % (np.mean(allDiff_y), ),
        r'$Max=%.2f$' % (maxResid, ),
        r'$\sigma=%.2f$' % (np.std(allDiff_y), )))
      axes[1].text(0.05, 0.95, textstr, transform=axes[1].transAxes, fontsize=14,
        verticalalignment='top')
      axes[1].set_xlim(-75,75)
      axes[1].legend(fontsize=12)
      fig.suptitle(plotTitle, fontsize=20)
      fig.tight_layout()
      plt.show()

    return maxResid


  def plotAvgPulseDiffs(self,fitPulseSetDict,binsDict,chWfBin=[],plotTitle=""):
    avg_pulse_x = []
    avg_pulse_y = []
    avg_pulseDiff_x = []
    avg_pulseDiff_y = []
    prev_avg_center = None
    prev_avg_mean = None
    slope = None
    minDiff = 10000.
    for binNum in binsDict :
      if len( binsDict[binNum]["samples"] ) == 0 : continue
      if binsDict[binNum]["sampMean"] == 0 : continue
      avg_pulse_x.append(binsDict[binNum]["center"])
      avg_pulse_y.append(binsDict[binNum]["sampMean"])
      if prev_avg_center != None :
        #do something here
        #avg_pulseDiff_x.append(binsDict[binNum]["center"])
        #avg_pulseDiff_y.append( binsDict[binNum]["sampMean"] - prev_avg_mean )
        if slope != None :
          predMean = prev_avg_mean + slope*(binsDict[binNum]["center"]-prev_avg_center)
          diff = binsDict[binNum]["sampMean"] - predMean
          avg_pulseDiff_x.append(binsDict[binNum]["center"])
          avg_pulseDiff_y.append( diff )
          if binsDict[binNum]["center"] > -20 and binsDict[binNum]["center"] < 10 :
            if diff < minDiff :
              minDiff = diff
        slope = (binsDict[binNum]["sampMean"] - prev_avg_mean)/( binsDict[binNum]["center"] - prev_avg_center )
      prev_avg_center = binsDict[binNum]["center"]
      prev_avg_mean = binsDict[binNum]["sampMean"]
      #if binNum - 1 in binsDict :
      #  avg_pulseDiff_x.append(binsDict[binNum]["center"])
      #  avg_pulseDiff_y.append( binsDict[binNum]["sampMean"] - binsDict[binNum-1]["sampMean"] )

    allPulse_fit_x = []
    allPulse_y = []
    for pulseSet in fitPulseSetDict :
      pulseSetInfo = fitPulseSetDict[pulseSet]
      for pulseNum, pulseTimeInfo in enumerate(pulseSetInfo) :
        pulse_x = pulseTimeInfo["pulse_x"]
        pulse_y = pulseTimeInfo["pulse_y"]
        pulse_fit_x = pulseTimeInfo["pulse_fit_x"]
        for sampNum, samp in enumerate(pulse_x):
          allPulse_fit_x.append(pulse_fit_x[sampNum])
          allPulse_y.append(pulse_y[sampNum])

    if False :    
      fig, axes = plt.subplots(1,2,figsize=(13, 8))
      axes[0].plot( [x*25 for x in avg_pulse_x   ] ,avg_pulse_y,"o",label="Average")
      axes[0].plot( [x*25 for x in allPulse_fit_x] ,allPulse_y,".",label="Samples")
      axes[0].set_xlabel('Sample Time [ns]', fontsize=20)
      axes[0].set_ylabel('Sample Value [ADC]', fontsize=20)
      axes[0].set_title("Sample Value vs Sample Time", fontsize=20)
      axes[0].set_xlim(-75,75)
      axes[0].legend(fontsize=12)
      axes[1].plot( [x*25 for x in avg_pulseDiff_x]  , avg_pulseDiff_y ,".",label="Avg Diffs")
      axes[1].set_xlabel('Sample Time [ns]', fontsize=20)
      axes[1].set_ylabel('Avg - Prev Avg [ADC]', fontsize=20)
      axes[1].set_title("Average - Previous Average Value vs Average Sample Time", fontsize=20)
      #axes[1].set_xlim(-75,75)
      axes[1].set_xlim(-25,10)
      #axes[0].set_xlim(-50,-30)
      axes[1].legend(fontsize=12)
      fig.suptitle(plotTitle, fontsize=20)
      fig.tight_layout()
      plt.show()
      
    return minDiff


  def fitPulses(self,fitPulseSetDict,binsDict,pedVal=0,pedRms=1):
    # get average pulse
    avg_pulse_x = []
    avg_pulse_y = []
    for binNum in binsDict :
      avg_pulse_x.append(binsDict[binNum]["center"])
      avg_pulse_y.append(binsDict[binNum]["sampMean"] )
  
    #loop through pulse sets
    allFitAmps = []
    maxNum = 200
    for pulseSetNum,pulseSet in enumerate(fitPulseSetDict) :
      #print("PULSE SET ", pulseSet)
      pulseSetInfo = fitPulseSetDict[pulseSet]
      #loop through pulses
      for pulseNum, pulseTimeInfo in enumerate(pulseSetInfo) :
        rangeLow = 7
        rangeHigh = 14
        pulse_x = pulseTimeInfo["pulse_x"][rangeLow:rangeHigh]
        pulse_y = pulseTimeInfo["pulse_y"][rangeLow:rangeHigh]
        pulse_fit_x = pulseTimeInfo["pulse_fit_x"][rangeLow:rangeHigh]
        pulse_ind_x = pulseTimeInfo["pulse_ind_x"][rangeLow:rangeHigh]
        pulse_time = pulseTimeInfo["pulseTime"]
        #print(   pulse_time )
        plot_avg_x = []
        plot_avg_y = []
        for sampNum,samp in enumerate(avg_pulse_x) :
          plot_avg_x.append( avg_pulse_x[sampNum] + pulse_time )
          plot_avg_y.append( avg_pulse_y[sampNum])
        #do pedestal correction
        if True:
          plot_avg_y = [x - pedVal for x in plot_avg_y]
          pulse_y = [x - pedVal for x in pulse_y]
        amp = self.fitPulse(pulse_x,pulse_fit_x,pulse_ind_x,pulse_y,plot_avg_x,plot_avg_y,0,1)
        #print( pulseSet, pulseNum, amp )
        allFitAmps.append(amp)
        if len( allFitAmps ) > maxNum : break
      if len( allFitAmps ) > maxNum : break
    #end loop
    print(allFitAmps)
    return allFitAmps
    
  def fitPulse(self,plot_x=[],plot_fit_x=[],plot_ind_x=[],plot_y=[],avg_pulse_x=[],avg_pulse_y=[],pedVal=0,pedRms=1):
  
    #things to fit: plot_ind_x,plot_y vs avg_pulse_x,avg_pulse_y
    tempFit = TEMPLATE_FIT()
    tempFit.doTemplateFit(plot_ind_x,plot_y,avg_pulse_x,avg_pulse_y,pedVal,pedRms)
    return tempFit.amp
  
    if False :
      fig, axes = plt.subplots(1,1,figsize=(13, 8))
      axes.plot(plot_ind_x,plot_y,"o",label="Samples")
      axes.plot(avg_pulse_x,avg_pulse_y,"-",label="Average")
      axes.set_xlabel('Sample [Number]', fontsize=20)
      axes.set_ylabel('Sample Value [ADC]', fontsize=20)
      axes.set_title("Sample Value vs Sample Number", fontsize=20)
      #axes.set_xlim(-75,75)
      #axes.set_xlim(-50,-30)
      axes.legend(fontsize=12)
      fig.tight_layout()
      #plt.draw()
      #plt.pause(1)
      plt.show()
    return None
    

#END CLASS

def main():
  print("HELLO")
  return

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
