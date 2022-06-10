import json
import sys
import subprocess
import numpy as np
from math import *
#import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import statsmodels.api as sm
import scipy.stats

from iminuit import Minuit
from iminuit.cost import LeastSquares

#BEGIN SLICE_ANALYZE_WAVEFORM CLASS
class TEMPLATE_FIT(object):

  #__INIT__#
  def __init__(self):
    self.data_x = None
    self.data_y = None
    self.temp_x = None
    self.temp_y = None
    self.pred_x = None
    self.pred_y = None
    self.predData_x = None
    self.predData_y = None
    self.pedVal = None
    return
    
  #required
  #method to generate new hypthesis from template given fit parameters
  #method to get hypoithesis preedicted value at specific time (sample number)
  
  def makePred(self,base = 0, timeOffset = 0 , amp = 0):
    self.pred_x = []
    self.pred_y = []
    for sampNum, samp in enumerate(self.temp_x):
      val = self.temp_y[sampNum]*amp + base
      self.pred_x.append( self.temp_x[sampNum] + timeOffset )
      self.pred_y.append( val )
    
    self.predData_x = []
    self.predData_y = []
    for sampNum,sampTime in enumerate(self.data_x) :
      val = self.getPredValue(sampTime)
      self.predData_x.append( sampTime )
      self.predData_y.append( val )
    
    return None
    
  def getPredValue(self,time=0):
    foundPred = False
    predSampNum = None
    #print(time, self.pred_x[0])
    if time <= self.pred_x[0] :
      return self.pred_y[0]
    for sampNum,sampTime in enumerate(self.pred_x):
      if sampNum == 0 : 
        continue
      if time > sampTime :
        continue
      foundPred = True
      predSampNum = sampNum
      break
    if foundPred == False :
      return self.pred_y[0]
    #if here, found predicted elements corresponding to requested time
    predTimeNext = self.pred_x[predSampNum]
    predTimePrev = self.pred_x[predSampNum-1] #will always work because sampNum == 0 skipped
    predValNext = self.pred_y[predSampNum]
    predValPrev = self.pred_y[predSampNum-1] #will always work because sampNum == 0 skipped
    timeDiff = time - predTimePrev
    if predValNext == predValPrev :
      return predValPrev
    if predTimeNext == predTimePrev :
      print("ERROR IN TEMPLATE")
      return predValPrev
    predVal = predValPrev + timeDiff*(predValNext - predValPrev)/float(predTimeNext - predTimePrev )
    return predVal
    
  def makePredGetVal(self,timeVals=[0],base = 0, timeOffset = 0 , amp = 1):
    #print("REQUESTED",time,base,timeOffset,amp)
    self.makePred(base,timeOffset,amp)
    results = []
    for time in timeVals :
      results.append( self.getPredValue(time) )
    return results
    
  def doTemplateFit(self,plot_x=[],plot_y=[],avg_pulse_x=[],avg_pulse_y=[],pedVal=0,pedRms=1.):
    self.data_x = plot_x
    self.data_y = plot_y
    self.data_y_err = [pedRms for x in self.data_y]
    self.temp_x = avg_pulse_x
    self.temp_y = avg_pulse_y
    self.pedVal = pedVal
    self.makePred(0,0,1)
    #self.makePred(base = 0, timeOffset = -5 , amp = 2)
    
    #print("START")
    least_squares = LeastSquares(self.data_x, self.data_y, self.data_y_err, self.makePredGetVal)
    m = Minuit(least_squares, base=0, timeOffset=0,amp=1)  # starting values for α and β
    m.fixed["base"] = True
    m.migrad()  # finds minimum of least_squares function
    
    #test_x = []
    #test_y = []
    #for sampNum, sampTime in enumerate(self.data_x):
    #  test_x.append(sampTime)
    #  test_y.append( self.makePredGetVal([sampTime],0,0,1) )
    
    #for p, v, e in zip(m.parameters, m.values, m.errors):
    #  #fit_info.append(f"{p} = ${v:.3f} \\pm {e:.3f}$")
    #  print(p,v,e)
    #print("AMP",m.values["amp"])
    
    self.amp = round(m.values["amp"],6)
    
    if False :
      fig, axes = plt.subplots(1,1,figsize=(13, 8))
      axes.plot(plot_x,plot_y,"o",label="Samples")
      axes.plot(avg_pulse_x,avg_pulse_y,"-",label="Average")
      axes.plot(self.pred_x,self.pred_y,"-",label="Template Prediction")
      #axes.plot(self.predData_x,self.predData_y,"x",label="Data Prediction")
      #axes.plot(test_x,test_y,"-",label="TEST")
      axes.set_xlabel('Sample [Number]', fontsize=20)
      axes.set_ylabel('Sample Value [ADC]', fontsize=20)
      axes.set_title("Sample Value vs Sample Number", fontsize=20)
      axes.legend(fontsize=12)
      fig.tight_layout()
      plt.show()
    return None
    

#END CLASS

def main():
  print("HELLO")
  return

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
