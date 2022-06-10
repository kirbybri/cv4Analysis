import numpy as np
from math import *
import matplotlib.pyplot as plt
import sys
import glob
import pickle
import statsmodels.api as sm
import scipy.stats
import json

from scipy.optimize import curve_fit

from cv4AnalyzeSample import CV4_ANALYZE_SAMPLE

def gaussian(x, amp, cen, wid):
    #print(x,amp,cen,wid)
    return amp * exp(-1.*(x-cen)*(x-cen) / (2*wid*wid))

#BEGIN CV3TB_ANALYZE_DACSCAN CLASS
class CV4_ANALYZE_DACSCAN(object):

  #__INIT__#
  def __init__(self, fileName = None):
    self.runResultsDict = None
    self.fileNameDict = None
    self.dacDataDict = None
    self.runAvgDataDict = None
    self.fileName = fileName
    
    #sample handling object
    self.analyzeSample = CV4_ANALYZE_SAMPLE()
    self.analyzeSample.dropOverFlowSamples = False
    self.analyzeSample.applyDdpuCorr = False
    self.analyzeSample.sampScaleFactor = 1.0
    self.analyzeSample.sampOffsetVal = 0.
    
    self.lowRun = -1000
    self.highRun = 1000
    self.rmsCut = 1.5
    self.applyCuts = False

  def measureLinearity(self,xs,ys,ysErr,lowLim,upLim):
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
      resid_yRms.append(ysErr[num])
    chiSq = chiSq / float( len(ysFit ) - 2 )

    print( "SLOPE ", slope , "\tERR ", slopeErr,"\tCHI2 ",chiSq )
    print( results.summary() )
    return slope, intercept, slopeErr, interceptErr, chiSq,resid_x,resid_y,resid_yRms

  def getSarVsDacVals(self,chId=None,runData=None):
    if chId == None or runData == None:
      print("ERROR,getSarVsDacVals invalid inputs")
      return None
    x = []
    y = []
    yRms = []
    for measNum in runData:
      measInfo = runData[measNum]
      if "data" not in measInfo or "attrs" not in measInfo:
        continue
      measData = measInfo["data"]
      measAttrs = measInfo["attrs"]
      if chId not in measData:
        continue
      measNumVal = measNum.split("_")
      measNumVal = int(measNumVal[1])
      dacVal = int(measNumVal)
      if "dacAVal" not in measAttrs :
        continue
      dacVal = int(measAttrs["dacAVal"])
      #dacVolt = (65535-dacVal)*0.036622 #mV
      #dacVolt = (dacVal)*0.036622 #mV
      dacVolt = 0.0364008*dacVal -1193.0563 #board 153 calib
      vals = self.getMeasChData(chId=chId,measNum=measNum)
      if len(vals) == 0 :
        continue
      avgVal = int(np.mean(vals))
      stdDev = float(np.std(vals))
      #apply cuts
      if self.applyCuts == True :
        #if measNumVal < self.lowRun or measNumVal > self.highRun :
        #  continue
        #if stdDev == 0 or stdDev > self.rmsCut :
        if stdDev == 0 :
          continue
      #x.append(dacVal)
      #if avgVal < 0 : avgVal = 0
      #if avgVal > 32767 : avgVal = 32767
      x.append(dacVolt)
      y.append(avgVal)
      yRms.append(stdDev)
      #if stdDev > 5 :
      #if True :
      #  print(measNumVal,"\t",dacVal,"\t",dacVolt,"\t",avgVal,"\t",chWf[0][16:],"\t",vals[0])
    #end for l
    orderX = np.argsort(x)
    xOrd = np.array(x)[orderX]
    yOrd = np.array(y)[orderX]
    yRmsOrd = np.array(yRms)[orderX]
    return xOrd,yOrd,yRmsOrd

  def plotMdacBits(self,chId=None,plotTitle=""):
    if chId == None:
      return None
    if self.runResultsDict == None:
      print("ERROR, results dict not defined")
      return
    if "results" not in self.runResultsDict:
      print("ERROR, no results in dict")
      return
    runData = self.runResultsDict["results"]

    x_dacVal = []
    y_mdacRange = []
    for measNum in runData:
      measInfo = runData[measNum]
      if "data" not in measInfo or "attrs" not in measInfo: continue
      measData = measInfo["data"]
      measAttrs = measInfo["attrs"]
      if chId not in measData: continue
      measNumVal = measNum.split("_")
      measNumVal = int(measNumVal[1])
      if "dacAVal" not in measAttrs : continue
      dacVal = measAttrs["dacAVal"]
      if "wf" not in measData[chId] : continue
      colutaWf = measData[chId]["wf"]
      
      #print("DAC VALUE ", dacVal)
      mdacVals = []
      for sampNum,sampInt in enumerate(colutaWf) :
        samp = self.analyzeSample.convertIntTo32BitWord(sampInt)
        header = samp[0:2]
        clockCheck = samp[2:4]
        mdacBits = samp[4:12]
        sarBits = samp[12:32]
        #print("\t",mdacBits)
        #convert thermometer curve to integer
        mdacCount = 0
        for bit in mdacBits :
          mdacCount = mdacCount + int(bit)
        #print("\t\t",mdacCount)
        mdacVals.append(mdacCount)

      #dacVolt = (65535-dacVal)*0.036622 #mV
      #dacVolt = (dacVal)*0.036622 #mVde
      dacVolt = 0.0364008*dacVal -1193.0563 #board 153 calib

      x_dacVal.append(dacVolt)
      y_mdacRange.append(np.mean(mdacVals))
      #end measNum loop
    if len(x_dacVal) == 0 : return None

    fig, axes = plt.subplots(figsize=(8, 6))
    axes.plot(x_dacVal,y_mdacRange,".")
    for num in range(0,8,1):    
      plt.plot([-1000+num*250+125, -1000+num*250+125], [-0.25, 8.25], lw=2, color='red', marker='', linestyle='dashed')
    axes.set_xlabel('DAC Voltage,diff [mV]', horizontalalignment='right', x=1.0, fontsize=16)
    axes.set_ylabel('Avg. MDAC Range', horizontalalignment='center', x=1.0, fontsize=16)
    axes.set_title("MDAC Range vs DAC Code", fontsize=20)
    axes.tick_params(axis="x", labelsize=16)
    axes.tick_params(axis="y", labelsize=16)
    axes.grid()

    fig.suptitle(plotTitle, fontsize=16)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

    return None

  def plotDacLinearityData(self,chId=None,plotTitle=""):
    if chId == None:
      return None
    if self.runResultsDict == None:
      print("ERROR, results dict not defined")
      return
    if "results" not in self.runResultsDict:
      print("ERROR, no results in dict")
      return
    runData = self.runResultsDict["results"]

    print("plotDacLinearityData: Using MDAC",self.analyzeSample.mdacWeights)
    print("plotDacLinearityData: Using SAR",self.analyzeSample.sarWeights)

    x,y,yRms = self.getSarVsDacVals(chId,runData)
    if len(x) != len(y) :
      print("ERROR, plotDacLinearityData couldn't get sample vs DAC values")
      return None
    #lineResult = self.measureLinearity(x,y,yRms,self.lowRun,self.highRun)
    lineResult = self.measureLinearity(x,y,yRms,self.lowRun+100,self.highRun-100)
    #lineResult = None
    X_plotFit = []
    Y_plotFit = []
    resid_x = []
    resid_y = []
    resid_yRms = []
    textstr = ""
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
      X_plotFit = np.linspace(self.lowRun,self.highRun,1000)
      Y_plotFit = X_plotFit*slope + intercept
      textstr = '\n'.join((
        r'$m=%.3f\pm%.3f$' % (slope, slopeErr, ),
        r'$b=%.2f\pm%.2f$' % (intercept,interceptErr, ),
        r'$\chi^2=%.2f$' % (chiSq, )
      ))
      #print("GAIN ", slope/-0.0366 , "ADCs/mV" )
      print("GAIN ", slope , "ADCs/mV" )
      cutResid = []
      for residNum,resid_xVal in enumerate(resid_x) :
        #if resid_xVal < self.lowRun or resid_xVal > self.highRun : continue
        cutResid.append( resid_y[residNum] )
      print("RESID RMS ",np.std(cutResid))

      #calculate optimal scale and offsets
      #optimal slope = 32767 / 2000mV,
      #optimal offset = -100mV*slope
      scaleFactor = 1.
      if slope != 0.:
        scaleFactor = (32767 / 2000.)/slope
      print("REC SCALE FACTOR",scaleFactor)
      offset = - scaleFactor*intercept
      print("REC OFFSET",offset)

    print("START PLOT")

    fig, axes = plt.subplots(2,2,figsize=(10, 6))
    axes[0][0].plot(x,y,".")
    if lineResult != None :
      axes[0][0].plot(X_plotFit,Y_plotFit,"-")
    for num in range(0,8,1):    
      axes[0][0].plot([-1000+num*250+125, -1000+num*250+125], [0, 33000], lw=2, color='red', marker='', linestyle='dashed')
    #axes[0][0].set_xlabel('DAC Code', horizontalalignment='right', x=1.0)
    axes[0][0].set_xlabel('DAC Diff. Voltage [mV DE]', horizontalalignment='right', x=1.0, fontsize=20)
    axes[0][0].set_ylabel('Sample Value [ADC]', horizontalalignment='center', x=1.0, fontsize=20)
    #axes[0][0].set_title("COLUTA SAmple Value vs DAC Code")
    axes[0][0].set_title("Mean Waveform Sample Value vs DAC Voltage", fontsize=20)
    #axes[0][0].set_xlim(self.lowRun,self.highRun)
    axes[0][0].set_ylim(-2000,32767+2000)
    axes[0][0].tick_params(axis="x", labelsize=16)
    axes[0][0].tick_params(axis="y", labelsize=16)
    axes[0][0].set_xlim(-1200,1200)
    axes[0][0].grid()
    axes[0][0].text(0.05, 0.95, textstr, transform=axes[0][0].transAxes, fontsize=14, verticalalignment='top', bbox=dict(facecolor='white', alpha=1.0, edgecolor='black'))

    axes[0][1].plot(x,yRms,".")
    #axes[0][1].set_xlabel('DAC Code', horizontalalignment='right', x=1.0)
    for num in range(0,8,1):    
      axes[0][1].plot([-1000+num*250+125, -1000+num*250+125], [0, 10], lw=2, color='red', marker='', linestyle='dashed')
    axes[0][1].set_xlabel('DAC Diff Voltage [mV DE]', horizontalalignment='right', x=1.0, fontsize=20)
    axes[0][1].set_ylabel('Waveform RMS [ADC]', horizontalalignment='center', x=1.0, fontsize=20)
    axes[0][1].set_xlim(-1200,1200)
    axes[0][1].tick_params(axis="x", labelsize=12)
    axes[0][1].tick_params(axis="y", labelsize=12)
    #axes[0][1].set_ylim(0,5.)
    axes[0][1].set_ylim(0,20.)
    axes[0][1].grid()    
    axes[0][1].set_title("Waveform RMS vs DAC Voltage", fontsize=20)

    if lineResult != None :
      axes[1][0].plot(resid_x,resid_y,".")
      #axes[1][0].errorbar(x=resid_x, y=resid_y, yerr=resid_yRms, fmt='.', ecolor='g', capthick=2)
    for num in range(0,8,1):    
      axes[1][0].plot([-1000+num*250+125, -1000+num*250+125], [-10, 10], lw=2, color='red', marker='', linestyle='dashed')
    axes[1][0].set_xlabel('DAC Diff. Voltage [mV DE]', horizontalalignment='right', x=1.0, fontsize=20)
    axes[1][0].set_ylabel('Residual [ADC]', horizontalalignment='center', x=1.0, fontsize=20)
    axes[1][0].set_xlim(-1200,1200)
    #axes[1][0].set_xlim(self.lowRun,self.highRun)
    #axes[1][0].set_ylim(-5,5)
    axes[1][0].tick_params(axis="x", labelsize=12)
    axes[1][0].tick_params(axis="y", labelsize=12)
    axes[1][0].grid()
    axes[1][0].set_title("Fit Residuals vs DAC Voltage", fontsize=20)

    if len(yRms) > 1:
      axes[1][1].hist(yRms,bins = np.arange(min(yRms), max(yRms)+1, 0.5))
    axes[1][1].set_xlabel('Waveform RMS', horizontalalignment='right', x=1.0, fontsize=20)
    axes[1][1].set_ylabel('Number of Waveforms', horizontalalignment='center', x=1.0, fontsize=20)
    axes[1][1].set_title("Waveform RMS Distribution", fontsize=20)
    axes[1][1].set_xlim(0,5)
    axes[1][1].tick_params(axis="x", labelsize=12)
    axes[1][1].tick_params(axis="y", labelsize=12)
    axes[1][1].grid()
    fig.suptitle(plotTitle, fontsize=16)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

    return

  def getFftWaveform(self,vals):
    fft_wf = np.fft.fft(vals)
    fftWf_x = []
    fftWf_y = []
    psd = []
    psd_x = []
    for sampNum,samp in enumerate(fft_wf) :
      if sampNum > float( len(fft_wf) ) / 2. :
        continue
      freqVal = 40. * sampNum/float(len(fft_wf))
      sampVal = np.abs(samp)
      if sampNum == 0 :
        sampVal = 0
      fftWf_x.append(freqVal)
      fftWf_y.append(sampVal)
    if np.max(fftWf_y) <= 0 :
      return psd_x,psd

    #fourier_fftWf_y = fftWf_y/np.max(fftWf_y)
    fourier_fftWf_y = fftWf_y
    for sampNum,samp in enumerate(fourier_fftWf_y) :
      if sampNum == 0 :
        continue
      else:
        psd_x.append( fftWf_x[sampNum] )
        psd.append( 20*np.log10(samp) )
    return psd_x,psd

  def getMeasChData(self,chId=None,measNum=None):
    if chId == None or self.runResultsDict == None or measNum == None :
      return None
    if "results" not in self.runResultsDict:
      print("ERROR, MISSING RESULTS")
      return None
    runData = self.runResultsDict["results"]
    if measNum not in runData :
      print("ERROR, MISSING MEASNUM",measNum)
      return None
    measInfo = runData[measNum]
    if "data" not in measInfo:
      print("ERROR, MISSING DATA")
      return None
    measData = measInfo["data"]
    if chId not in measData:
      print("ERROR, MISSING CHID",chId)
      return None
    chanData = measData[chId]
    if "attrs" not in measInfo:
      print("ERROR, MISSING ATTRS")
      return None
    measAttrs = measInfo["attrs"]
    if '32BitMode' not in measAttrs:
      print("ERROR, MISSING ATTRS 32BITMODE")
      return None
    chWf = chanData["wf"]
    if measAttrs['32BitMode'] == chId :
      #print("DO 32BIT MODE CORRECTION")
      chWf = self.analyzeSample.getWaveformVals(chWf)
    #account for bad samples
    #if self.dropInitialSamples :
    #  chWf = chWf[self.numDroppedInitialSamples:]
    return chWf
  
  def viewDist(self,measNum=None,chId=None, measInfo=None):
    if measNum == None or chId == None or measInfo == None:
      return None
    if "data" not in measInfo or "attrs" not in measInfo:
      return None
    measData = measInfo["data"]
    measAttrs = measInfo["attrs"]
    if chId not in measData:
      return None
    measNumVal = measNum.split("_")
    measNumVal = int(measNumVal[1])
    dacVal = int(measNumVal)
    if "dacAVal" in measAttrs :
      dacVal = int(measAttrs["dacAVal"])
    #dacVolt = (65535-dacVal)*0.0366 #mV
    dacVolt = 0.0364008*dacVal -1193.0563 #board 153 calib
    vals = self.getMeasChData(chId=chId,measNum=measNum)
    #print( measNum , dacVolt , vals[0:3] )  

    allMean = np.mean(vals)
    allRms = np.std(vals)

    if allMean < 100 or allMean > 32767 : return None
    
    fig, axes = plt.subplots(1,1,figsize=(12, 6))
    histTotal, bins, patches = axes.hist( vals, bins = np.arange(int(min(vals)), int(max(vals))+1, 1), density=False )
    centers = (0.5*(bins[1:]+bins[:-1]))
    
    #print("centers", centers)
    #print("histtotal",histTotal)
    
    #do fit
    v_gaussian = np.vectorize(gaussian)
    pars, cov = curve_fit(v_gaussian, centers, histTotal, p0=[max(histTotal),allMean,allRms])
    print("FIT PARS  ", pars )
    
    #get fit plot
    fit_x = []
    fit_y = []
    #plotRange = np.arange(int(-10)-0.5, int(10)+1, 1)
    plotRange = bins
    #for point in centers:
    for point in plotRange:
      for num in range(0,10,1):
        fit_x.append(point+num/10.)
        fit_y.append( gaussian(point+num/10.,pars[0],pars[1],pars[2]) )
    
    #get chisq
    chiSq = 0
    num = 0
    for pointNum,point in enumerate(centers):
      if histTotal[pointNum] == 0 :
        continue
      num = num + 1
      fitVal = gaussian(point,pars[0],pars[1],pars[2])
      err = np.sqrt( histTotal[pointNum] )
      resid = histTotal[pointNum] - fitVal
      chiSq = chiSq + resid*resid/err/err
    if num > 3:
      chiSq = chiSq / float( num - 3 )
    
    errMu = 0.
    if np.sqrt(cov[1,1]) != inf :
      errMu = np.sqrt(cov[1,1])
    errSigma = 0.
    if np.sqrt(cov[2,2]) != inf :
      errSigma = np.sqrt(cov[2,2])
    textstr = '\n'.join((
      r'$\mu=%.2f\pm%.2f$' % (pars[1], errMu, ),
      r'$\sigma=%.2f\pm%.2f$' % (np.abs(pars[2]),errSigma, ),
      r'$\chi^2 / ndf=%.2f / %d$' % (float(chiSq),float(num), )
    ))      
    
    #axes.plot(centers,histTotal,marker=".")
    axes.plot( fit_x,fit_y)
    axes.text(0.05, 0.95, textstr, transform=axes.transAxes, fontsize=14, verticalalignment='top')
    axes.set_xlabel('Sample Value [ADC counts]', horizontalalignment='right', x=1.0)
    axes.set_ylabel('# Samples / ADC code', horizontalalignment='left', x=1.0)
    axes.set_title("Waveform Sample Distribution, "+str(chId) + ", DAC " + str(round(dacVolt,2)) +"mVdiff" )
    axes.set_yscale('log')
    fig.tight_layout()
    plt.show()
    
    return None
  
  def viewDacScanDist(self,chId=None):
    if chId == None:
      print("NO CHANNEL ID")
      return None
    if self.runResultsDict == None:
      print("NO RUNS DICT")
      return
    if "results" not in self.runResultsDict:
      print("NO RUNS RESULT")
      return
    runData = self.runResultsDict["results"]

    for measNum in runData:
      measInfo = runData[measNum]
      self.viewDist(measNum,chId,measInfo)
    return

  def viewDacScanWaveforms(self,chId=None):
    if chId == None:
      print("NO CHANNEL ID")
      return None
    if self.runResultsDict == None:
      print("NO RUNS DICT")
      return
    if "results" not in self.runResultsDict:
      print("NO RUNS RESULT")
      return
    runData = self.runResultsDict["results"]

    for measNum in runData:
      measInfo = runData[measNum]
      if "data" not in measInfo or "attrs" not in measInfo:
        continue
      measData = measInfo["data"]
      measAttrs = measInfo["attrs"]
      if chId not in measData:
        continue
      measNumVal = measNum.split("_")
      measNumVal = int(measNumVal[1])
      dacVal = int(measNumVal)
      if "dacAVal" in measAttrs :
        dacVal = int(measAttrs["dacAVal"])
      dacVolt = (65535-dacVal)*0.0366 #mV
      #vals = self.cv3tbAnalyzeSample.getWaveformVals(chWf)
      vals = self.getMeasChData(chId=chId,measNum=measNum)

      #apply some cuts
      stdDev = float(np.std(vals))
      #if self.applyCuts == True :
      #  #if stdDev == 0 or stdDev > self.rmsCut :
      #  #if stdDev == 0 :
      #  #  continue
      #if stdDev ==0 or dacVolt > 250: continue
      #if stdDev ==0 or dacVolt < 2100 or dacVolt > 2200: continue
          
      psd_x,psd = self.getFftWaveform(vals)

      textstr = '\n'.join((
        r'$\mu=%.2f$' % ( np.mean(vals), ),
        r'$\sigma=%.2f$' % ( np.std(vals), )
      ))  

      fig, axes = plt.subplots(1,2,figsize=(10, 6))
      axes[0].plot(vals,".")
      axes[0].set_xlabel('Sample #', horizontalalignment='right', x=1.0)
      axes[0].set_ylabel('ADC CODE [ADC]', horizontalalignment='left', x=1.0)
      axes[0].set_title("COLUTA WAVEFORM" )
      axes[0].text(0.05, 0.95, textstr, transform=axes[0].transAxes, fontsize=14, verticalalignment='top')
      axes[1].plot(psd_x,psd,".")
      axes[1].set_xlabel('Frequency [MHz]', horizontalalignment='right', x=1.0)
      axes[1].set_ylabel('PSD [dB]', horizontalalignment='left', x=1.0)
      axes[1].set_title(str(measNum) + ", Voltage " + str(dacVolt)+"mV")
      #axes[1].set_xlim(0,0.5)
      fig.tight_layout()
      plt.show()

  def plotBitVsDac(self,chId=None,bitToPlot=None):
    if chId == None or bitToPlot == None:
      return None
    if self.runResultsDict == None:
      return
    if "results" not in self.runResultsDict:
      return
    runData = self.runResultsDict["results"]
    x = []
    y = []
    for measNum in runData:
      measInfo = runData[measNum]
      if "data" not in measInfo or "attrs" not in measInfo:
        continue
      measData = measInfo["data"]
      measAttrs = measInfo["attrs"]
      if chId not in measData:
        continue
      measNumVal = measNum.split("_")
      measNumVal = int(measNumVal[1])
      chWf = measData[chId]
      bitVals = []
      for samp in chWf :
        bitVal = samp[31-bitToPlot]
        bitVals.append( int(bitVal) )
      avgBitVal = np.mean(bitVals)
      x.append(measNumVal)
      y.append(avgBitVal)

    orderX = np.argsort(x)
    xOrd = np.array(x)[orderX]
    yOrd = np.array(y)[orderX]
    prevAvgBitVal = 0
    foundTransition = False
    transitionMeasNum = 0
    for cnt in range(0,len(xOrd),1):
      measNumVal = xOrd[cnt]
      avgBitVal = yOrd[cnt]
      #print(measNumVal ,"\t", avgBitVal)
      if cnt > 0 and foundTransition == False:
        if avgBitVal < 0.5 and prevAvgBitVal >= 0.5 :
          print(bitToPlot,"\t",measNumVal,"\t",prevAvgBitVal,"\t",avgBitVal)
          foundTransition = True
          transitionMeasNum = measNumVal
      prevAvgBitVal = avgBitVal

    fig, axes = plt.subplots(1,1,figsize=(10, 6))
    axes.plot(x,y,".")
    #axes.set_xlim(transitionMeasNum-500,transitionMeasNum+500)
    axes.set_xlabel('Measurement #', horizontalalignment='right', x=1.0)
    #axes.set_xlabel('DAC DC Differential Voltage [mV]', horizontalalignment='right', x=1.0)
    axes.set_ylabel('Avg Bit Value', horizontalalignment='left', x=1.0)
    #axes.set_title("Avg Bit Value vs Measurement #, bit # " + str(bitToPlot) )
    axes.set_title("Avg Bit Value vs DAC Voltage, bit # " + str(bitToPlot) )
    fig.tight_layout()
    plt.show()    
    return

  def printMeasSamples(self,chId=None):
    if chId == None:
      return None
    if self.runResultsDict == None:
      return
    if "results" not in self.runResultsDict:
      return
    runData = self.runResultsDict["results"]
    measVals = []
    chSampleStr = []
    for measNum in runData:
      measInfo = runData[measNum]
      if "data" not in measInfo or "attrs" not in measInfo:
        continue
      measData = measInfo["data"]
      measAttrs = measInfo["attrs"]
      if chId not in measData:
        continue
      measNumVal = measNum.split("_")
      measNumVal = int(measNumVal[1])

      if "dacAVal" not in measAttrs : continue

      chWf = measData[chId]
      #vals = self.getWaveformVals(chWf)
      #vals = self.cv3tbAnalyzeSample.getWaveformVals(chWf)
      #vals = self.getMeasChData(chId=chId,measNum=measNum)
      #stdDev = float(np.std(vals))
      measVals.append(measNumVal)

      binSamps = []
      for sampNum,samp in enumerate(chWf["wf"]):
        #binSamps.append(self.analyzeSample.convertIntTo32BitWord(int(samp)))
        binSamps.append( '{0:032b}'.format(samp) )
        if sampNum > 101 : break

      #chSampleStr.append(str(measNumVal)+","+",".join(chWf[0:100]))
      #chSampleStr.append(str(measAttrs["dac_aval"])+","+",".join(chWf[0:100]))
      chSampleStr.append(str(measAttrs["dacAVal"])+","+",".join(binSamps[0:100]))
      #chSampleStr.append(str(measNumVal)+","+",".join(chWf))
      #print("\t",str(measNumVal),"\t",",".join(chWf[0:10]) )
    orderMeasVals = np.argsort(measVals)
    measValsOrd = np.array(measVals)[orderMeasVals]
    chSampleStrOrd = np.array(chSampleStr)[orderMeasVals]
    f = open("output_cv4AnalyzeDacScan_data.txt", "w")
    for strVal in chSampleStrOrd:
      f.write(strVal+"\n")
    f.close()
    return

  def getMeasAvgSarValForMdacVal(self,chId=None,measNum=None,mdacVal=None):
    if chId == None or measNum == None or mdacVal == None:
      return None
    if self.runResultsDict == None:
      return
    if "results" not in self.runResultsDict:
      return
    runData = self.runResultsDict["results"]
    if measNum not in runData:
      return
    measData = runData[measNum]
    if chId not in measData:
      return
    chWf = measData[chId]
    zeroMdacWeights = [0,0,0,0,0,0,0,0]
    vals = []
    for samp in chWf:
      mdacBits = samp[4:12]
      sarBits = samp[12:32]
      print(samp,"\t",mdacBits,"\t",sarBits)
      sarVal = self.cv3tbAnalyzeSample.getColutaSampleValue(sarBits,zeroMdacWeights)
      if mdacBits == mdacVal :
        vals.append(sarVal)
    avgVal = float(np.mean(vals))
    print(measNum,"\t",mdacVal,"\t",len(vals),"\t",avgVal)
    return avgVal

  def dumpFile(self):
    if self.runResultsDict == None:
      return
    if "results" not in self.runResultsDict:
      return
    runData = self.runResultsDict["results"]
    for measNum in runData:
      measData = runData[measNum]
      print("Measurement ", measNum)
      print("\t",measData)
    return

  #open file
  def openFile(self):
    if self.fileName == None :
      print("ERROR no input file specified")
      return None
    self.runResultsDict = pickle.load( open( self.fileName, "rb" ) )
    return

  ##very specific analysis to look at samples around MDAC transitions
  def analyzeMdacTransition(self,chId=None,mdacBit=None):
    if chId == None or mdacBit==None :
      return None
    if self.runResultsDict == None:
      return
    if "results" not in self.runResultsDict:
      return
    bitToPlot = mdacBit + 20
    runData = self.runResultsDict["results"]
    x_measNum = []
    x_measDac = []
    y_mdacBit = []
    y_sarVal = []
    y_sarVal_low = []
    y_sarVal_high = []
    x_sarVal_low_static = []
    x_sarVal_high_static = []
    y_sarVal_low_static = []
    y_sarVal_high_static = []
    yRms_sarVal_low_static = []
    yRms_sarVal_high_static = []
    tempMdacWeights = [0,0,0,0,0,0,0,0]
    mdacLowList = ["00000000","00000001","00000011","00000111","00001111","00011111","00111111","01111111"]
    mdacHighList = ["00000001","00000011","00000111","00001111","00011111","00111111","01111111","11111111"]
    mdacLow = mdacLowList[mdacBit]
    mdacHigh = mdacHighList[mdacBit]
    #print("HERE")
    for measNum in runData:
      measInfo = runData[measNum]
      if "data" not in measInfo or "attrs" not in measInfo:
        continue
      measData = measInfo["data"]
      measAttrs = measInfo["attrs"]
      if chId not in measData:
        continue
      measNumVal = measNum.split("_")
      measNumVal = int(measNumVal[1])
      dacVal = int(measNumVal)
      if "dac_aval" in measAttrs :
        dacVal = int(measAttrs["dac_aval"])
      dacVolt = (65535-dacVal)*0.0366 #mV
      #dacVolt = (dacVal)*0.0366 #mV
      chWf = measData[chId]
      bitVals = []
      sarVals = []
      sarVals_low = []
      sarVals_high = []
      for samp in chWf :
        samp = self.cv3tbAnalyzeSample.convertIntTo32BitWord(samp)
        bitVal = samp[31-bitToPlot]
        mdacBits = samp[4:12]
        sarBits = samp[12:32]
        if mdacBits != mdacLow and mdacBits != mdacHigh :
          continue
        bitVals.append(int(bitVal))
        sarVal = self.cv3tbAnalyzeSample.getColutaSampleValue(sarBits,tempMdacWeights)
        sarVal = int(sarVal*self.cv3tbAnalyzeSample.sampScaleFactor/4.)
        sarVals.append(sarVal)
        if mdacBits == mdacLow :
          sarVals_low.append(sarVal)
        if mdacBits == mdacHigh :
          sarVals_high.append(sarVal)
      if len(bitVals) == 0 : continue
      #print(dacVal,np.mean(bitVals) , bitVals)
      measNumVal = dacVolt
      x_measNum.append(measNumVal)
      x_measDac.append(dacVolt)
      y_mdacBit.append( np.mean(bitVals) )
      y_sarVal.append( np.mean(sarVals) )
      y_sarVal_low.append( np.mean(sarVals_low) )
      y_sarVal_high.append( np.mean(sarVals_high) )
      if np.mean(bitVals) == 0 :
        x_sarVal_low_static.append(measNumVal)
        y_sarVal_low_static.append( np.mean(sarVals_low) )
        yRms_sarVal_low_static.append( np.std(sarVals_low) )
      if np.mean(bitVals) == 1 :
        x_sarVal_high_static.append(measNumVal)
        y_sarVal_high_static.append( np.mean(sarVals_high) )
        yRms_sarVal_high_static.append( np.std(sarVals_high) )
      #print(measNum,dacVal, np.mean(bitVals) )

    print( y_mdacBit )
    orderX = np.argsort(x_measNum)
    xOrd = np.array(x_measNum)[orderX]
    xOrd_measDac = np.array(x_measDac)[orderX]
    yOrd_mdacBit = np.array(y_mdacBit)[orderX]
    yOrd_sarVal = np.array(y_sarVal)[orderX]
    yOrd_sarVal_low = np.array(y_sarVal_low)[orderX]
    yOrd_sarVal_high = np.array(y_sarVal_high)[orderX]

    orderX = np.argsort(x_sarVal_low_static)
    xOrd_sarVal_low_static = np.array(x_sarVal_low_static)[orderX][1:-1]
    yOrd_sarVal_low_static = np.array(y_sarVal_low_static)[orderX][1:-1]
    yRmsOrd_sarVal_low_static = np.array(yRms_sarVal_low_static)[orderX][1:-1]

    orderX = np.argsort(x_sarVal_high_static)
    xOrd_sarVal_high_static = np.array(x_sarVal_high_static)[orderX][1:-1]
    yOrd_sarVal_high_static = np.array(y_sarVal_high_static)[orderX][1:-1]
    yRmsOrd_sarVal_high_static = np.array(yRms_sarVal_high_static)[orderX][1:-1]

    prevAvgBitVal = 0
    foundTransition = False
    transitionMeas = None
    for cnt in range(0,len(xOrd),1):
      measNumVal = xOrd[cnt]
      avgBitVal = yOrd_mdacBit[cnt]
      print(measNumVal, avgBitVal)
      if cnt > 0 and foundTransition == False:
        if (avgBitVal >= 0.5 and prevAvgBitVal < 0.5) or (avgBitVal < 0.5 and prevAvgBitVal >= 0.5) :
          print(bitToPlot,"\t",measNumVal,"\t",prevAvgBitVal,"\t",avgBitVal,"\t",yOrd_sarVal_low[cnt],"\t",yOrd_sarVal_high[cnt],"\t",yOrd_sarVal_low[cnt]-yOrd_sarVal_high[cnt])
          transitionMeas = measNumVal
          foundTransition = True
      prevAvgBitVal = avgBitVal

    plotLowLim = transitionMeas-400
    plotHighLim = transitionMeas+400

    slope_low, intercept_low, slopeErr, interceptErr, chiSq,resid_x_low,resid_y_low,resid_yRms = self.measureLinearity(xOrd_sarVal_low_static,yOrd_sarVal_low_static,yRmsOrd_sarVal_low_static,plotLowLim,plotHighLim)
    X_plotFit_low = np.linspace(plotLowLim,plotHighLim,1000)
    Y_plotFit_low = X_plotFit_low*slope_low + intercept_low

    slope_high, intercept_high, slopeErr, interceptErr, chiSq,resid_x_high,resid_y_high,resid_yRms = self.measureLinearity(xOrd_sarVal_high_static,yOrd_sarVal_high_static,yRmsOrd_sarVal_high_static,plotLowLim,plotHighLim)
    X_plotFit_high = np.linspace(plotLowLim,plotHighLim,1000)
    Y_plotFit_high = X_plotFit_high*slope_high + intercept_high

    val_low = slope_low*transitionMeas + intercept_low
    val_high = slope_high*transitionMeas + intercept_high
    print("\t",val_low,"\t",val_high,"\t",val_low-val_high)
    print("MDAC BIT         " , mdacBit )
    print("TRANSITION POINT " , transitionMeas)
    print("LOW VAL          ", val_low )
    print("HIGH VAL         ", val_high )
    print("DIFF             ", val_low-val_high )

    fig, axes = plt.subplots(1,3,figsize=(10, 6))
    axes[0].plot(xOrd,yOrd_mdacBit,".")
    if transitionMeas != None :
      axes[0].set_xlim(plotLowLim,plotHighLim)
    #axes[0].set_xlabel('Measurement #', horizontalalignment='right', x=1.0)
    axes[0].set_xlabel('DAC Voltage [mV]', horizontalalignment='right', x=1.0)
    axes[0].set_ylabel('Avg Bit Value', horizontalalignment='left', x=1.0)
    axes[0].set_title("Avg Bit Value vs DAC Voltage, MDAC bit # " + str(mdacBit) )
    axes[1].plot(xOrd,yOrd_sarVal,".",label="Avg SAR Bit Value")
    if transitionMeas != None :
      axes[1].set_xlim(plotLowLim,plotHighLim)
    #axes[1].set_xlabel('Measurement #', horizontalalignment='right', x=1.0)
    axes[1].set_xlabel('DAC Voltage [mV]', horizontalalignment='right', x=1.0)
    axes[1].set_ylabel('Avg SAR Bit Value [ADC]', horizontalalignment='left', x=1.0)
    axes[1].set_title("Avg SAR Bit Value vs DAC Voltage: MDAC bit # " + str(mdacBit) )
    axes[1].plot(xOrd,yOrd_sarVal_low,".",label="MDAC Code " + str(mdacLow))
    axes[1].plot(xOrd,yOrd_sarVal_high,".",label="MDAC Code " + str(mdacHigh))
    #axes[1].plot(xOrd_sarVal_low_static,yOrd_sarVal_low_static,".")
    axes[1].plot(X_plotFit_low,Y_plotFit_low,"-")
    #axes[1].plot(xOrd_sarVal_high_static,yOrd_sarVal_high_static,".")
    axes[1].plot(X_plotFit_high,Y_plotFit_high,"-")
    axes[1].legend()
    
    axes[2].set_title("Residual vs DAC Voltage: MDAC bit # " + str(mdacBit) )
    axes[2].set_xlabel('DAC Voltage [mV]', horizontalalignment='right', x=1.0)
    axes[2].set_ylabel('Fit Residual [ADC]', horizontalalignment='left', x=1.0)
    axes[2].plot(resid_x_low,resid_y_low,".")
    axes[2].plot(resid_x_high,resid_y_high,".")
    
    fig.tight_layout()
    #plt.show()    
    return val_low-val_high, transitionMeas

  def viewMeasWaveform(self,chId=None,measNum=None):
    if chId == None or measNum == None:
      return None
    if self.runResultsDict == None:
      return
    if "results" not in self.runResultsDict:
      return
    runData = self.runResultsDict["results"]
    if measNum not in runData:
      return
    measInfo = runData[measNum]
    if "data" not in measInfo or "attrs" not in measInfo:
      return
    measData = measInfo["data"]
    measAttrs = measInfo["attrs"]
    if chId not in measData:
      return
    chWf = measData[chId]

    #print(chWf)

    x = []
    y = []
    x_same = []
    y_same = []
    x_diff = []
    y_diff = []

    prevMdacBits = None
    for sampNum,samp in enumerate(chWf) :
      mdacBits = samp[4:12]
      sarBits = samp[12:32]
      #sarVal = self.cv3tbAnalyzeSample.getColutaSampleValue(sarBits,[0,0,0,0,0,0,0,0])
      sarVal = self.cv3tbAnalyzeSample.getColutaSampleValue(sarBits,mdacBits)
      x.append(sampNum)
      y.append(sarVal)
      if mdacBits == prevMdacBits or sampNum == 0:
        x_same.append(sampNum)
        y_same.append(sarVal)
      else :
        x_diff.append(sampNum)
        y_diff.append(sarVal)
      prevMdacBits = mdacBits

    fig, axes = plt.subplots(1,1,figsize=(10, 6))
    #axes.plot(x,y,".")
    axes.plot(x_same,y_same,".",label="Sample with no MDAC transition", markersize=16)
    axes.plot(x_diff,y_diff,".",label="Sample with MDAC transition", markersize=16)
    axes.set_xlabel('Sample #', horizontalalignment='right', x=1.0)
    axes.set_ylabel('Sample Value Using SAR Bits Only', horizontalalignment='left', x=1.0)
    axes.set_title("COLUTA Sample Waveform, SAR Bit Values Only, " + str(measNum) )
    #axes.set_ylim(2700,2900)
    #axes.set_ylim(6200,6600)
    #axes.legend(fontsize=12)
    fig.tight_layout()
    plt.show()    
    return

  def analyzeMdacTransitionEffect(self,chId=None):
    if chId == None:
      return None
    if self.runResultsDict == None:
      return
    if "results" not in self.runResultsDict:
      return
    runData = self.runResultsDict["results"]

    x_measNum = []
    y_offset = []
    #loop over measurements
    for measNum in runData:
      measData = runData[measNum]
      if chId not in measData:
        continue
      chWf = measData[chId]
      measNumVal = measNum.split("_")
      measNumVal = int(measNumVal[1])
      x_same = []
      samp_same = []
      x_diff = []
      samp_diff = []
      prevMdacBits = None
      for sampNum,samp in enumerate(chWf) :
        mdacBits = samp[4:12]
        sarBits = samp[12:32]
        sarVal = self.cv3tbAnalyzeSample.getColutaSampleValue(sarBits,[0,0,0,0,0,0,0,0])
        if mdacBits == prevMdacBits:
          samp_same.append(sarVal)
          x_same.append(sampNum)
        elif sampNum > 0:
          samp_diff.append(sarVal)
          x_diff.append(sampNum)
        prevMdacBits = mdacBits
      if len(samp_diff) == 0 : #require some transition samples
        continue
      if len(samp_same) / len(chWf ) < 0.9 : #require most samples to be in same MDAC range
        continue

      baseVal = np.mean(samp_same)
      x_offset = []
      samp_offset = []
      for sampNum,samp in enumerate(samp_diff) :
        if abs(samp - baseVal) > 1000 :
          continue
        x_offset.append( x_diff[sampNum] )
        samp_offset.append(samp)
      offsetVal = np.mean(samp_offset)
      x_measNum.append(measNumVal)
      y_offset.append( offsetVal - baseVal )

      print(measNum,"\t", baseVal ,"\t",offsetVal,"\t", offsetVal - baseVal)
      continue
      #continue

      fig, axes = plt.subplots(1,1,figsize=(10, 6))
      axes.plot(x_same,samp_same,".",label="No MDAC transition")
      axes.plot(x_offset,samp_offset,".",label="MDAC transition")
      axes.set_xlabel('Sample #', horizontalalignment='right', x=1.0)
      axes.set_ylabel('Sample Value Using SAR Bits Only', horizontalalignment='left', x=1.0)
      axes.set_title("COLUTA Sample Waveform, SAR Bit Values Only, " + str(measNum) )
      axes.legend()
      #axes.set_ylim(baseVal - 200, baseVal + 200 )
      fig.tight_layout()
      plt.show()    
      #end meas loop

    fig, axes = plt.subplots(1,2,figsize=(10, 6))
    axes[0].plot(x_measNum,y_offset,".")
    axes[0].set_xlabel('Measurement #', horizontalalignment='right', x=1.0)
    axes[0].set_ylabel('Sample Offset Due to MDAC Transition', horizontalalignment='left', x=1.0)
    axes[0].set_title("Sample Offsets due to MDAC Transition vs Measurement #" )
    #axes[1].hist(y_offset,bins = np.arange(min(y_offset), max(y_offset)+1))
    axes[1].hist(y_offset)
    axes[1].set_xlabel('Sample Offset Due to MDAC Transition', horizontalalignment='right', x=1.0)
    axes[1].set_ylabel('Number of Samples', horizontalalignment='left', x=1.0)
    axes[1].set_title("Sample Offset due to MDAC Transition Distribution" )
    #axes.set_ylim(baseVal - 200, baseVal + 200 )
    fig.tight_layout()
    plt.show()   

    return

  def analyzeDacScanDacDataFile(self):
    if self.fileName == None :
      print("ERROR no input file specified")
      return None

    measNum = []
    dacAVolts = []
    dacBVolts = []
    diffVolts = []
    diffVoltsRms = []
    with open(self.fileName) as fp:
     for cnt, line in enumerate(fp):
       if cnt == 0:
         continue
       #if cnt == 2:
       #  break
       line = line.split("\t")
       #print(line)
       dacAVolt = float(line[0])*1000.
       dacBVolt = float(line[1])*1000.
       measNum.append(cnt-1)
       dacAVolts.append(dacAVolt)
       dacBVolts.append(dacBVolt)
       diffVolts.append(dacBVolt-dacAVolt)
       diffVoltsRms.append(0.00005*1000.)

    slope, intercept, slopeErr, interceptErr, chiSq,resid_x,resid_y,resid_yRms = self.measureLinearity(measNum,diffVolts,diffVoltsRms,self.lowRun,self.highRun)
    X_plotFit = np.linspace(self.lowRun,self.highRun,1000)
    Y_plotFit = X_plotFit*slope + intercept

    textstr = '\n'.join((
      r'$m=%.4f\pm%.4f$' % (slope, slopeErr, ),
      r'$b=%.2f\pm%.2f$' % (intercept,interceptErr, ),
      r'$\chi^2=%.2f$' % (chiSq, )
    ))

    fig, axes = plt.subplots(1,2,figsize=(10, 6))
    axes[0].plot(measNum,dacAVolts,".",label="DAC A voltage")
    axes[0].plot(measNum,dacBVolts,".",label="DAC B voltage")
    axes[0].plot(measNum,diffVolts,".",label="Differential voltage")
    axes[0].plot(X_plotFit,Y_plotFit,".",label="Differential voltage fit")
    axes[0].set_xlabel('Measurement #', horizontalalignment='right', x=1.0)
    axes[0].set_ylabel('DAC voltage [mV]', horizontalalignment='left', x=1.0)
    axes[0].set_title("DAC voltage vs Measurement #" )
    axes[0].text(0.05, 0.45, textstr, transform=axes[0].transAxes, fontsize=14, verticalalignment='top')
    #axes.set_ylim(baseVal - 200, baseVal + 200 )
    axes[0].legend()
    
    axes[1].plot(resid_x,resid_y,".",label="Differential Voltage Fit Residual")
    axes[1].set_xlabel('Measurement #', horizontalalignment='right', x=1.0)
    axes[1].set_ylabel('Fit Residual [mV]', horizontalalignment='left', x=1.0)
    axes[1].set_title("Fit Residual vs Measurement #" )
    fig.tight_layout()
    plt.show()   

    return

  def measureDnl(self,chId=None):
    if chId == None:
      return None
    if self.runResultsDict == None:
      return
    if "results" not in self.runResultsDict:
      return
    runData = self.runResultsDict["results"]

    #allVals = []
    #allVals = np.ushort([])
    allValsList = []
    totMeas = 0
    for measNum in runData:
      measInfo = runData[measNum]
      if "data" not in measInfo or "attrs" not in measInfo:
        continue
      measData = measInfo["data"]
      measAttrs = measInfo["attrs"]
      if chId not in measData:
        continue
      measNumVal = measNum.split("_")
      measNumVal = int(measNumVal[1])
      if measNumVal % 1000 == 0 :
        print("MEAS NUM", measNumVal)

      #chWf = measData[chId]
      #vals = self.getWaveformVals(chWf)
      #vals = self.cv3tbAnalyzeSample.getWaveformVals(chWf)
      vals = self.getMeasChData(chId=chId,measNum=measNum)
      if len(vals) == 0 :
        continue
      #print("\t",measNumVal,"\t",len(vals),"\t",measAttrs["dac_aval"])
      #print(measAttrs)
      totMeas = totMeas + 1
      #allVals = np.append( allVals, np.ushort(vals) )
      allValsList.append(np.ushort(vals))
      #print( all
      #for samp in vals:
      # allVals.append(  np.short(samp)  )
      #  allVals = np.append( allVals , np.short(samp) )

    allVals = np.array( allValsList , dtype = np.ushort )

    binEdges = []
    for edge in range(-1000,32767+1000,1):
      binEdges.append(edge)

    codeBins, bin_edges = np.histogram(allVals,bins=binEdges, density=False)
    center = (bin_edges[:-1] + bin_edges[1:]) / 2
    if len(codeBins) != len(bin_edges) - 1: 
      print("ERROR DNL")
      return
      
    data = {}
    for binNum, binCount in enumerate(codeBins):
      data[ int(bin_edges[binNum]) ] = int(binCount)
    
    with open("outputData_CV4_dacScanHist.json", 'w') as fp:
          json.dump(data, fp, sort_keys=True, indent=4)
    return
  
    #print(len(codeBins),len(bin_edges))
    #for num, sampBin in enumerate(codeBins):
    #  print(bin_edges[num],"\t",bin_edges[num+1],"\t",codeBins[num])
    #  #print(bin_edges[num])
    #return

    if False :
    
      center = (bin_edges[:-1] + bin_edges[1:]) / 2
      plot_x = []
      plot_y = []
      for codeNum,code in enumerate(center):
        if code < 15000 or code > 16000 :
          continue
        plot_x.append( code )
        plot_y.append( codeBins[codeNum] )
        print(codeNum,code,center[codeNum],codeBins[codeNum])
    
      fig, axes = plt.subplots(1,1,figsize=(10, 6))
      #n, bins, patches = axes.hist(allVals,bins = np.arange(min(allVals), max(allVals)+1))
     
      width = 0.7 * (bin_edges[1] - bin_edges[0])
      plt.bar(plot_x, plot_y, align='center',width=width )
      #axes.plot(plot_x, plot_y,".")
      #axes.set_xlim(6950,7050)
      fig.tight_layout()
      plt.show()
    
    #find minimum code with non-zero # entries
    """
    minBin = 0
    for binNum, binVal in enumerate(codeBins):
      if binVal == 0 : continue
      minBin = binNum
      break

    #find max code with non-zero # entries
    maxBin = 32767+1000
    for binNum, binVal in enumerate(reversed(codeBins)):
      if binVal == 0 : continue
      maxBin = binNum
      break
    maxBin = 32767+1000 - maxBin - 2
    """
    minBin = 2000
    maxBin = 31000
    print(minBin,binEdges[minBin],binEdges[minBin+1],center[minBin],codeBins[minBin])
    print(maxBin,binEdges[maxBin],binEdges[maxBin+1],center[maxBin],codeBins[maxBin])
    avgVal = 0
    avgValSimple = 0
    numBins = 0
    for binNum, binVal in enumerate(codeBins):
      if binNum < minBin : continue
      if binNum > maxBin : continue
      avgVal = avgVal + (binVal - avgVal)/float(numBins+1)
      avgValSimple = avgValSimple  + binVal
      numBins = numBins + 1
    avgValSimple = avgValSimple / float( numBins) 

    dnl_x = []
    dnl_y = []
    inl_y = []
    inlSum = 0
    for binNum, binVal in enumerate(codeBins):
      if binNum < minBin : continue
      if binNum > maxBin : continue
      
      binDnl = (binVal - avgVal)/avgVal
      
      dnl_x.append(bin_edges[binNum])
      dnl_y.append(binDnl)
      inlSum = inlSum + binDnl
      inl_y.append(inlSum)
      
      if binDnl <= -1 : 
        print( binNum, binVal , binDnl )
    

      
    #fig_dnl, axes_dnl = plt.subplots(1,1,figsize=(10, 6))
    #axes_dnl.plot(dnl_x,dnl_y,marker='o',linestyle='')
    #axes_dnl.set_xlabel('COLUTA Code [ADC]', fontsize=20)
    #axes_dnl.set_ylabel('DNL from Normal Mode DAC Scan', fontsize=20)
    #axes_dnl.tick_params(axis="x", labelsize=12)
    #axes_dnl.tick_params(axis="y", labelsize=12)
    #axes_dnl.set_title("DNL vs COLUTA Code Value" )
    #fig_dnl.tight_layout()
    #plt.show()   
      
    return
    
    print(minBin,"\t",codeBins[minBin])
    print(maxBin,"\t",codeBins[maxBin])
    if minBin >= maxBin :
      print("ERROR")
      return
    #count # of bins in buffer region
    offset = 500
    #numSamples = 8184
    #numSamples = 128
    numSamples = 512
    minRangeSum = 0
    for binNum in range(minBin-1,minBin+offset,1):
      minRangeSum = minRangeSum + codeBins[binNum]
    maxRangeSum = 0
    for binNum in range(maxBin-offset,maxBin+1,1):
      maxRangeSum = maxRangeSum + codeBins[binNum]
    minBin = minBin + offset
    maxBin = maxBin - offset
    numBins = maxBin - minBin + 1
    expNumBins = (totMeas*numSamples - minRangeSum - maxRangeSum )/float(numBins)


    inlSum = 0
    minDnl = 10
    minBinVal = 0
    minBinNum = -1
    for binNum, binVal in enumerate(codeBins):
      if binNum < minBin : continue
      if binNum > maxBin : continue
      binDnl = (binVal - expNumBins)/expNumBins
      if binDnl < minDnl :
        minDnl = binDnl
        minBinVal = binVal
        minBinNum = binNum
      #binDnl = (binVal - expNumBins)/expNumBins
      #print(binNum,"\t",bin_edges[binNum],"\t",binVal,"\t",binDnl)   
      dnl_x.append(bin_edges[binNum])
      dnl_y.append(binDnl)
      inlSum = inlSum + binDnl
      inl_y.append(inlSum)

    print(minBinNum,"\t",minDnl,"\t",minBinVal,"\t",expNumBins)

    #print(len(hist),"\t",len(bin_edges))
    #for binNum, binVal in enumerate(hist):
    #  print(binNum,"\t",bin_edges[binNum],"\t",binVal)

    """
    fig, axes = plt.subplots(1,2,figsize=(10, 6))
    #axes.hist(allVals,bins = np.arange(min(allVals), max(allVals)+1))
    axes[0].plot(dnl_x,dnl_y,marker='o',linestyle='')
    axes[0].set_xlabel('COLUTA Code Value [ADC]', horizontalalignment='right', x=1.0)
    axes[0].set_ylabel('DNL [LSB]', horizontalalignment='left', x=1.0)
    axes[0].set_title("DAC Scan DNL vs COLUTA Code Value" )
    #axes[0].set_xlim(minBin+1000,maxBin-1000)
    axes[1].plot(dnl_x,inl_y)
    axes[1].set_xlabel('COLUTA Code Value [ADC]', horizontalalignment='right', x=1.0)
    axes[1].set_ylabel('INL [LSB]', horizontalalignment='left', x=1.0)
    axes[1].set_title("DAC Scan INL vs COLUTA Code Value" )
    #axes[1].set_xlim(minBin+1000,maxBin-1000)
    #axes[1].set_ylim(-10,10)
    fig.tight_layout()
    """
        
    hist, bin_edges = np.histogram(dnl_y, bins = np.linspace(-1.1, 3.1, num=85) )
    print("HIST\t")
    print(*bin_edges, sep=", ")
    print(*hist, sep=", ")

    #textstr_dnl = '\n'.join((
    #  #r'$\mu=%.3f$' % ( np.mean(dnl_y) ),
    #  r'$\sigma=%.2f$' % ( np.std(dnl_y) )
    #))
    textstr_dnl = r'$\sigma=%.2f$' % ( np.std(dnl_y) )

    textstr_inl = '\n'.join((
      #r'$\mu=%.3f$' % ( np.mean(inl_y) ),
      r'$\sigma=%.2f$' % ( np.std(inl_y) )
    ))

    fig_inl, axes_inl = plt.subplots(1,1,figsize=(10, 6))
    #axes_inl.plot(dnl_x,inl_y)
    axes_inl.set_xlabel('COLUTA Code [ADC]', fontsize=20)
    axes_inl.set_ylabel('INL', fontsize=20)
    axes_inl.tick_params(axis="x", labelsize=12)
    axes_inl.tick_params(axis="y", labelsize=12)
    #axes_inl.set_title("DAC Scan INL vs COLUTA Code Value" )
    axes_inl.set_ylim(-15,15)
    fig_inl.tight_layout()

    fig_dnlHist, axes_dnlHist = plt.subplots(1,1,figsize=(10, 6))
    axes_dnlHist.hist(dnl_y,bins = np.arange(min(dnl_y), max(dnl_y)+1, 0.1))
    axes_dnlHist.set_xlabel('DNL from Normal Mode DAC Scan', fontsize=20)
    axes_dnlHist.set_ylabel('Number of Codes', fontsize=20)
    #axes_dnlHist.set_title("DAC Scan DNL Distribution" )
    axes_dnlHist.text(0.85, 0.75, textstr_dnl, transform=axes_dnl.transAxes, fontsize=14, verticalalignment='top')
    axes_dnlHist.tick_params(axis="x", labelsize=12)
    axes_dnlHist.tick_params(axis="y", labelsize=12)
    axes_dnlHist.set_yscale('log')
    fig_dnlHist.tight_layout()

    fig_inlHist, axes_inlHist = plt.subplots(1,1,figsize=(10, 6))
    #axes_inlHist.hist(inl_y,bins = np.arange(min(inl_y), max(inl_y)+1, 0.1))
    axes_inlHist.set_xlabel('INL', fontsize=20)
    axes_inlHist.set_ylabel('Number of Codes', fontsize=20)
    #axes_inlHist.set_title("DAC Scan INL Distribution" )
    #axes_inlHist.text(0.05, 0.45, textstr_inl, transform=axes_inl.transAxes, fontsize=14, verticalalignment='top')
    axes_inlHist.tick_params(axis="x", labelsize=12)
    axes_inlHist.tick_params(axis="y", labelsize=12)
    axes_inlHist.set_yscale('log')
    fig_inlHist.tight_layout()

    
    
    return
    
  """
  def scaleWeights(self):
    totWeight = np.sum(self.mdacWeights) + np.sum(self.sarWeights)
    if totWeight <= 0 :
      return
    scaleFactor = 32768.0 / totWeight
    for weightNum,weight in enumerate(self.mdacWeights):
      self.mdacWeights[weightNum] = self.mdacWeights[weightNum]*scaleFactor
    for weightNum,weight in enumerate(self.sarWeights):
      self.sarWeights[weightNum] = self.sarWeights[weightNum]*scaleFactor
  """

  def analyzeMdacConstants(self,chId=None):
    if chId == None:
      return None
    if self.runResultsDict == None:
      return
    if "results" not in self.runResultsDict:
      return
    runData = self.runResultsDict["results"]
    
    print(self.cv3tbAnalyzeSample.applyDdpuCorr)
    print(self.cv3tbAnalyzeSample.sarWeights)
    print(self.cv3tbAnalyzeSample.sampScaleFactor)

    mdacCodeSarRange = {}
    mdacCodeAvgBit = {}
    #loop over measurements
    for measNum in runData:
      measInfo = runData[measNum]
      if "data" not in measInfo or "attrs" not in measInfo:
        continue
      measData = measInfo["data"]
      measAttrs = measInfo["attrs"]
      if chId not in measData:
        continue
      measNumVal = measNum.split("_")
      measNumVal = int(measNumVal[1])
      dacVal = int(measNumVal)
      if "dac_aval" in measAttrs :
        dacVal = int(measAttrs["dac_aval"])
      #dacVolt = (65535-dacVal)*0.0366 #mV
      dacVolt = (dacVal)*0.0366 #mV
      chWf = measData[chId]
      
      mdacCodes = {}
      for samp in chWf :
        sampBin = self.cv3tbAnalyzeSample.convertIntTo32BitWord(samp)
        header = sampBin[0:2]
        clockCheck = sampBin[2:4]
        mdacBits = sampBin[4:12]
        sarBits = sampBin[12:32]
        val = self.cv3tbAnalyzeSample.getColutaSampleValue(sarBits,mdacBits)
        sarVal = self.cv3tbAnalyzeSample.getColutaSampleValue(sarBits,[0,0,0,0,0,0,0,0])
        sarVal = int(sarVal*self.cv3tbAnalyzeSample.sampScaleFactor/4.)
        print(measNum,"\t",mdacBits,"\t",sarBits,"\t",sarVal)
        if mdacBits not in mdacCodes :
          mdacCodes[mdacBits] = []
        mdacCodes[mdacBits].append( sarVal )
      #end sample loop
      if len(mdacCodes) != 1 : continue
      mdacCode = list(mdacCodes)[0]
      meanSarVal = np.mean(list(mdacCodes.values())[0])
      rmsSarVal = np.std(list(mdacCodes.values())[0])
      if mdacCode not in mdacCodeSarRange :
        mdacCodeSarRange[mdacCode] = {}
        mdacCodeSarRange[mdacCode]["dacVals"] = []
        mdacCodeSarRange[mdacCode]["sarVals"] = []
        mdacCodeSarRange[mdacCode]["sarValsRms"] = []
      mdacCodeSarRange[mdacCode]["dacVals"].append(dacVal)
      mdacCodeSarRange[mdacCode]["sarVals"].append(meanSarVal)
      mdacCodeSarRange[mdacCode]["sarValsRms"].append(rmsSarVal)
    #end measurement loop
    
    #look over sar values for each mdac code
    print( mdacCodeSarRange.keys() )
    for mdacCode in mdacCodeSarRange :
      xs = mdacCodeSarRange[mdacCode]["dacVals"][1:-1]
      ys = mdacCodeSarRange[mdacCode]["sarVals"][1:-1]
      ysErr = mdacCodeSarRange[mdacCode]["sarValsRms"][1:-1]
      
      lineResult = self.measureLinearity(xs=xs,ys=ys,ysErr=ysErr,lowLim=0,upLim=100000)
      if True and lineResult != None :
        slope = lineResult[0]
        intercept = lineResult[1]
        mdacCodeSarRange[mdacCode]["slope"] = slope
        mdacCodeSarRange[mdacCode]["intercept"] = intercept
      
      if False :
        fig, axes = plt.subplots(1,2,figsize=(10, 6))
        axes[0].plot(xs,ys,".")
        if True and lineResult != None :
          slope = lineResult[0]
          intercept = lineResult[1]
          slopeErr = lineResult[2]
          interceptErr = lineResult[3]
          chiSq = lineResult[4]
          resid_x = lineResult[5]
          resid_y = lineResult[6]
          resid_yRms = lineResult[7]
        
          print( lineResult)
        
          X_plotFit = np.linspace(np.min(xs),np.max(xs),1000)
          Y_plotFit = X_plotFit*slope + intercept
          axes[0].plot(X_plotFit,Y_plotFit,"-")
          #.set_title("COLUTA SAmple Value vs DAC Voltage")
          #axes.set_xlim(self.lowRun,self.highRun)
          #axes.set_ylim(0,32768)
          #axes[0][0].set_xlim(0,2400)
          #axes[0][0].text(0.05, 0.95, textstr, transform=axes[0][0].transAxes, fontsize=14, verticalalignment='top')
          axes[1].plot(resid_x,resid_y)
          axes[1].set_xlabel('DAC Voltage [mV]', horizontalalignment='right', x=1.0)
          axes[1].set_ylabel('Fit Residual Value [ADC]', horizontalalignment='center', x=1.0)
        axes[0].set_xlabel('DAC Voltage [mV]', horizontalalignment='right', x=1.0)
        axes[0].set_ylabel('Avg. SAR Sample Value [ADC]', horizontalalignment='center', x=1.0)
        plt.show()
    
    #review results
    for mdacCode in mdacCodeSarRange :
      slope = mdacCodeSarRange[mdacCode]["slope"]
      intercept = mdacCodeSarRange[mdacCode]["intercept"]
      minDac = np.min(mdacCodeSarRange[mdacCode]["dacVals"])
      maxDac = np.max(mdacCodeSarRange[mdacCode]["dacVals"])
      print(mdacCode,"\t",minDac,"\t",maxDac,"\t",slope,"\t",intercept)

    return

def main():
  if len(sys.argv) != 2 :
    print("ERROR, program requires filename argument")
    return
  fileName = sys.argv[1]
  cv3tbAnalyzeFile = CV3TB_ANALYZE_DACSCAN(fileName)
  cv3tbAnalyzeFile.openFile()
  
#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
