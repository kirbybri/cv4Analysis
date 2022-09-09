import numpy as np
from math import *
import matplotlib.pyplot as plt
import json
import statsmodels.api as sm

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
      #print(xs[num],"\t",lowLim,"\t",upLim)
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

def main():
    #10MHz case
    x_data = [0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19]
    y_data = [13.714069790054815, 14.616268369865573, 15.43505399534083, 16.184347409060457, 16.87488971241477, 17.5151338433963, 18.11184622783493, 18.67052419658239, 19.195692595660866, 19.691119026858793]
    y_data_err = [0.5 for x in y_data]

    lowVal = 0.05
    highVal = 0.195
    lineResult = measureLinearity(xs=x_data,ys=y_data,ysErr=y_data_err,lowLim=lowVal,upLim=highVal)
    if lineResult == None :
      print("FIT FAILED")
      return
    X_plotFit = []
    Y_plotFit = []
    resid_x = []
    resid_y = []
    resid_yRms = []
    textstr = ""
    if True :
      #slope, intercept, slopeErr, interceptErr, chiSq,resid_x,resid_y,resid_yRms = linResult
      slope = lineResult[0]
      intercept = lineResult[1]
      slopeErr = lineResult[2]
      interceptErr = lineResult[3]
      chiSq = lineResult[4]
      resid_x = lineResult[5]
      resid_y = lineResult[6]
      resid_yRms = lineResult[7]
      X_plotFit = np.linspace(np.min(x_data),np.max(x_data),1000)
      Y_plotFit = X_plotFit*slope + intercept
      textstr = '\n'.join((
        r'$m=%.3f\pm%.3f$' % (slope, slopeErr, ),
        r'$b=%.2f\pm%.2f$' % (intercept,interceptErr, ),
        r'$\chi^2=%.2f$' % (chiSq, )
      ))
      #print("GAIN ", slope/-0.0366 , "ADCs/mV" )
      print("SLOPE", slope , "V/DAC" )
      cutResid = []
      for residNum,resid_xVal in enumerate(resid_x) :
        cutResid.append( resid_y[residNum] )
      print("RESID RMS ",np.std(cutResid))
    
    fig, axes = plt.subplots(1,1,figsize=(10, 9))
    axes.plot(x_data,y_data,".")
    axes.plot(X_plotFit,Y_plotFit,"-")
    axes.set_xlabel('Disturbance Amplitude [ADC]', horizontalalignment='right', x=1.0, fontsize=16)
    axes.set_ylabel('PSD Difference', horizontalalignment='center', x=1.0, fontsize=16)
    axes.tick_params(axis="x", labelsize=16)
    axes.tick_params(axis="y", labelsize=16)
    #axes[0].set_xlim(0,2400)
    axes.grid()
    axes.text(0.05, 0.95, textstr, transform=axes.transAxes, fontsize=14, verticalalignment='top')
    fig.tight_layout()
    plt.show()
    
    

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
