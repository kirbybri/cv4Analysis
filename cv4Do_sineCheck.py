import os
import json
from pathlib import Path
import sys
import numpy as np
from math import *
import pickle
import matplotlib.pyplot as plt
from iminuit import Minuit

from cv4ProcessFile import CV4_PROCESS_FILE
from cv4AnalyzeWaveform import CV4_ANALYZE_WAVEFORM
from cv4AnalyzeDacScan import CV4_ANALYZE_DACSCAN

timeVals = []
sampVals = []
sampErrs = []

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
    
def plotQuick2(vals_x,vals_y,label=""):
    fig, axes = plt.subplots(1,1,figsize=(13, 6))
    axes.plot(vals_x,vals_y,".")
    axes.set_xlabel('Sine Amp [ADC]', fontsize=20)
    axes.set_ylabel('3rd Harmonic Size [dB]', fontsize=20)
    axes.set_title("4.985MHz Sine Wave 3rd Harmonic VS Amplitude", fontsize=20)
    axes.tick_params(axis="x", labelsize=12)
    axes.tick_params(axis="y", labelsize=12)
    #axes[0].set_xlim(0,4000)
    #axes[0].set_ylim(17087-50,17087+50)    
    axes.grid()
    fig.suptitle(label, fontsize=16)
    fig.tight_layout()
    
    plt.show()    

def sineFunc(time,amp,freq,phase,offset):
    return amp*np.sin( time*freq + phase ) + offset

def getSimSineWave():
    sampFreq = 40.E+6 #Hz
    sampPeriod = 1/float(sampFreq) #s
    sineFreq = 5.024E+6 #Hz
    #sineFreq = 198.4E+3 #Hz
    sineRadFreq = sineFreq*2.*np.pi
    numSamp = 6250
    sineAmp = 10000.0 #?
    sinePhase = np.pi/3. #rad
    sineOffset = 0.
    noise = 1.5
    
    for num in range(0,numSamp,1):
        sampTime = num*sampPeriod
        #sampVal = sineAmp*np.sin( sampTime*sineRadFreq + sinePhase ) + sineOffset
        sampVal = sineFunc(sampTime,sineAmp,sineRadFreq,sinePhase,sineOffset)
        sampNoise = noise*np.random.randn()
        timeVals.append(sampTime)
        sampVals.append(sampVal+sampNoise)
        sampErrs.append(1.5)
    return timeVals,sampVals

def getResid(amp,freq,phase,offset):
    resid = []
    for sampNum,samp in enumerate(sampVals):
        sampTime = timeVals[sampNum]
        resid.append( samp - sineFunc(sampTime,amp,freq,phase,offset) )
        #print(sampNum,samp, sineFunc(sampTime,amp,freq,phase,offset), samp - sineFunc(sampTime,amp,freq,phase,offset) )
    return resid

def LSQ(amp,freq,phase,offset):
    ls = 0
    resids = getResid(amp,freq,phase,offset)
    for resid in resids:
        ls += resid*resid
    return ls / 5. / 5.

def SINAD(fourier):
    sum2 = 0
    for normBin in fourier:
        if normBin==1: continue
        sum2 += normBin**2
    return -10*np.log10(sum2)

def ENOB(fourier):
    return (SINAD(fourier)-1.76)/6.02 

def getFftWaveform(vals,sampFreq=40.):
    print("num Samples",len(vals))
    fft_wf = np.fft.fft(vals)
    fftWf_x = []
    fftWf_y = []
    psd = []
    psd_x = []
    phase = []
    for sampNum,samp in enumerate(fft_wf) :
        if sampNum > float( len(fft_wf) ) / 2. : continue
        freqVal = sampFreq * sampNum/float(len(fft_wf)) #MHz
        sampVal = np.abs(samp)
        if sampNum == 0 :
            sampVal = 0
        fftWf_x.append(freqVal)
        fftWf_y.append(sampVal)
        phase.append(np.angle(samp))
    if np.max(fftWf_y) <= 0 :
        return psd_x,psd,0,0
    fourier_fftWf_y = fftWf_y/np.max(fftWf_y)
    for sampNum,samp in enumerate(fourier_fftWf_y) :
        if sampNum == 0 :
            continue
        else:
            psd_x.append( fftWf_x[sampNum] )
            psd.append( 20*np.log10(samp) )
    sinad = SINAD(fourier_fftWf_y)
    enob = ENOB(fourier_fftWf_y)
    print ("ENOB\t",enob)
    return psd_x,psd,phase,sinad,enob

def plotVals(timeVals,sampVals,psd_x,psd):
    numSamp = 200
    #timeVals = timeVals[0:numSamp]
    #sampVals = sampVals[0:numSamp]
    print("# Samples ", len(sampVals) )
    fig, axes = plt.subplots(1,2,figsize=(10, 6))
    axes[0].plot(timeVals,sampVals,".")
    axes[0].set_xlabel('Sample Time [ns]', fontsize=20)
    axes[0].set_ylabel('Sample Value [ADC]',  fontsize=20)
    axes[0].tick_params(axis="x", labelsize=12)
    axes[0].tick_params(axis="y", labelsize=12)
    #axes[0].set_xlim(0,1000*25)
    axes[0].set_title("DATA")

    axes[1].plot(psd_x,psd,"-")
    axes[1].set_xlabel('Frequency [MHz]', fontsize=20)
    axes[1].set_ylabel('PSD [dB]', fontsize=20)
    axes[1].tick_params(axis="x", labelsize=12)
    axes[1].tick_params(axis="y", labelsize=12)

    axes[1].set_title("PSD")
    #axes[1].set_xlim(0,1)

    fig.tight_layout()
    plt.show()    

def plotResids(times,resids,psd_x,psd):
    numSamp = 200
    #timeVals = timeVals[0:numSamp]
    #sampVals = sampVals[0:numSamp]
    print("# Samples ", len(resids) )
    fig, axes = plt.subplots(1,2,figsize=(10, 6))
    axes[0].plot(times,resids,".")
    axes[0].set_xlabel('Sample Time [ns]', fontsize=20)
    axes[0].set_ylabel('Residual [ADC]',  fontsize=20)
    axes[0].tick_params(axis="x", labelsize=12)
    axes[0].tick_params(axis="y", labelsize=12)
    #axes[0].set_xlim(0,1000*25)
    axes[0].set_title("Fit Residuals")

    axes[1].plot(psd_x,psd,"-")
    axes[1].set_xlabel('Frequency [MHz]', fontsize=20)
    axes[1].set_ylabel('PSD [dB]', fontsize=20)
    axes[1].tick_params(axis="x", labelsize=12)
    axes[1].tick_params(axis="y", labelsize=12)

    axes[1].set_title("Fit Residuals PSD")
    #axes[1].set_xlim(0,1)

    fig.tight_layout()
    plt.show()    

def doSineFit(vals):
      for num,val in enumerate(vals):
        timeVals.append(num*25.E-9)
        sampVals.append(val)
      print(len(vals),len(sampVals))    
      psd_x,psd,phase,sinad,enob = getFftWaveform(vals)
      #plotVals([x*25E+9 for x in timeVals][0:400],sampVals[0:400],psd_x,psd)
      
      #infer sine parameters
      maxPsd = psd[1]
      maxFreq = psd_x[1]
      maxPhase = phase[1]
      for freqNum,freq in enumerate(psd_x):
        #print(freqNum,freq,phase[freqNum])
        if psd[freqNum] > maxPsd :
            maxPsd = psd[freqNum]
            maxFreq = freq
            maxPhase = phase[freqNum]
      maxFreq = maxFreq*1.E+6 #convert to Hz
      ampEst = (np.max(sampVals) - np.min(sampVals))/2.
      offsetEst = np.mean(sampVals)
      m = Minuit(LSQ, amp=ampEst, freq=maxFreq*2.*np.pi, phase=maxPhase-np.pi/2.,offset=offsetEst)
      #print( m.fixed )
      #m.fixed[0] = True
      #m.fixed[1] = True  
      #m.fixed[2] = True  
      m.fixed[3] = True
      m.errors[0] = 10.
      m.errors[1] = 10.
      m.errors[2] = 0.1
      #print( m.fixed )
      print(m.values)
      m.simplex()  # run optimiser 
      print(m.values)
      m.migrad()  # run optimiser
      print(m.values)  # x: 2, y: 3, z: 4
      m.hesse()   # run covariance estimator
      print(m.errors)  # x: 1, y: 1, z: 1
      #plotQuick(plot_x,plot_y,label="COLUTA Ch1,5.005MHz Sine Wave,8000 Samples")
      #print("AVG",np.mean(plot_y),"RMS",np.std(plot_y))
      fitAmp = m.values[0]
      fitFreq = m.values[1]
      fitPhase = m.values[2]
      fitOffset = m.values[3]
      print(fitAmp)
      resids = getResid(fitAmp,fitFreq,fitPhase,fitOffset)
      psd_x_resid,psd_resid,phase_resid,sinad_resid,enob_resid = getFftWaveform(resids)
      print( "RESIDS",np.mean(resids),np.std(resids))
      plotResids([x*25E+9 for x in timeVals][0:400],resids[0:400],psd_x_resid,psd_resid)    

def findMaxPsd(chFFt_x,chFft_y,lowFreq,highFreq):
  maxPsd = None
  maxPsd_x = None
  for psdNum,freq in enumerate(chFFt_x) :
    if freq < lowFreq : continue
    if freq > highFreq : continue
    psdVal = chFft_y[psdNum]
    #print( psdNum , freq , psdVal )
    if maxPsd == None :
      maxPsd = psdVal
      maxPsd_x = freq
    if psdVal > maxPsd :
      maxPsd = psdVal
      maxPsd_x = freq
  #print("MAX",maxPsd,maxPsd_x)
  return maxPsd,maxPsd_x        
        
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
    amp_x = []
    psdDiff = []
    for measNum in cv4AnalyzeWaveform.runResultsDict["results"] :
      awgAmp = float(cv4AnalyzeWaveform.runResultsDict["results"][measNum]["attrs"]['awgAmp'])
      #meanvals, stdvals, maxval,minval,enob,sinad = cv4AnalyzeWaveform.viewWaveform(chId=chanName,measNum=measNum,doPrint=False,doPlot=False)
      #print( measNum , enob )
      #plot_x.append(awgAmp)
      #plot_y.append(enob)
      #continue
      #fit below
      #timeVals,sampVals = getSineWave()
      timeVals.clear()
      sampVals.clear()
      vals = cv4AnalyzeWaveform.getMeasChData(chId=chanName,measNum=measNum)
      doSineFit(vals)
      return
      continue
    
    
      for num,val in enumerate(vals):
        timeVals.append(num*25.E-9)
        sampVals.append(val)
      
      #plotVals([x*25E+9 for x in timeVals][0:400],sampVals[0:400],psd_x,psd)
      #continue
      #print(len(vals),len(sampVals))    
      #psd_x,psd,phase,sinad,enob = getFftWaveform(vals)
      #maxPsd,maxPsd_x = findMaxPsd(psd_x,psd,14.5,15.5)
      #print(maxPsd,maxPsd_x)
      #amp_x.append( np.std(vals)*np.sqrt(2)*2 )
      #psdDiff.append(maxPsd)
      #plotVals([x*25E+9 for x in timeVals][0:400],sampVals[0:400],psd_x,psd)
    plotQuick2(amp_x,psdDiff)

  return

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
