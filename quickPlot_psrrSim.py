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
    #0.1MHz, ch1
    x_data = [0.05, 0.1, 0.15000000000000002, 0.2, 0.25, 0.30000000000000004, 0.35000000000000003, 0.4, 0.45, 0.5, 0.55, 0.6000000000000001, 0.65, 0.7000000000000001, 0.75, 0.8, 0.8500000000000001, 0.9, 0.9500000000000001]
    y_data = [12.1, 14.48, 18.04, 20.56, 22.51, 24.1, 25.45, 26.61, 27.64, 28.56, 29.39, 30.14, 30.84, 31.48, 32.08, 32.65, 33.17, 33.67, 34.14]
    lowVal = 0.55
    highVal = 0.75
    targetDiff = 30.5

    #0.1MHz, ch4
    x_data =  [0.05, 0.1, 0.15000000000000002, 0.2, 0.25, 0.30000000000000004, 0.35000000000000003, 0.4, 0.45, 0.5, 0.55, 0.6000000000000001, 0.65, 0.7000000000000001, 0.75, 0.8, 0.8500000000000001, 0.9, 0.9500000000000001]
    y_data =  [11.27, 14.72, 18.22, 20.71, 22.64, 24.22, 25.56, 26.72, 27.74, 28.66, 29.49, 30.24, 30.94, 31.58, 32.18, 32.74, 33.27, 33.76, 34.23]
    lowVal = 0.15
    highVal = 0.35
    targetDiff = 22.0

    #1MHz, ch1
    x_data =  [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19]
    y_data = [1.54, 2.43, 4.24, 6.58, 8.38, 9.93, 11.28, 12.46, 13.5, 14.43, 15.27, 16.04, 16.74, 17.4, 18.0, 18.57, 19.1, 19.61, 20.08]
    lowVal = 0.02
    highVal = 0.065
    targetDiff = 7.6

    #1MHz, ch4
    x_data =  [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19]
    y_data = [1.4, 2.51, 4.74, 6.77, 8.57, 10.11, 11.45, 12.61, 13.65, 14.57, 15.41, 16.17, 16.87, 17.52, 18.12, 18.69, 19.22, 19.72, 20.19]
    lowVal = 0.005
    highVal = 0.03
    targetDiff = 1.8

    #5MHz, ch1
    x_data = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19]
    y_data = [5.58, 7.6, 9.27, 10.7, 11.94, 13.03, 14.01, 14.89, 15.68, 16.42, 17.09, 17.72, 18.3, 18.85, 19.37, 19.85, 20.31, 20.75, 21.17]
    lowVal = 0.04
    highVal = 0.08
    targetDiff = 13.7

    #5MHz, ch4
    x_data = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19]
    y_data = [1.83, 3.56, 5.42, 7.68, 9.48, 10.96, 12.21, 13.3, 14.27, 15.14, 15.93, 16.65, 17.32, 17.94, 18.51, 19.06, 19.57, 20.05, 20.5]
    lowVal = 0.085
    highVal = 0.125
    targetDiff = 15.5

    #0.05MHz, ch1
    x_data = [0.05, 0.1, 0.15000000000000002, 0.2, 0.25, 0.30000000000000004, 0.35000000000000003, 0.4, 0.45, 0.5, 0.55, 0.6000000000000001, 0.65, 0.7000000000000001, 0.75, 0.8, 0.8500000000000001, 0.9, 0.9500000000000001]
    y_data = [12.1, 14.53, 18.05, 20.56, 22.51, 24.1, 25.44, 26.61, 27.63, 28.55, 29.38, 30.14, 30.84, 31.48, 32.08, 32.64, 33.17, 33.67, 34.14]
    lowVal = 0.55
    highVal = 0.75
    targetDiff = 30.6

    #0.05MHz, ch4
    x_data = [0.05, 0.1, 0.15000000000000002, 0.2, 0.25, 0.30000000000000004, 0.35000000000000003, 0.4, 0.45, 0.5, 0.55, 0.6000000000000001, 0.65, 0.7000000000000001, 0.75, 0.8, 0.8500000000000001, 0.9, 0.9500000000000001]
    y_data = [11.27, 14.9, 18.36, 20.82, 22.73, 24.3, 25.63, 26.78, 27.8, 28.71, 29.53, 30.28, 30.97, 31.61, 32.21, 32.77, 33.29, 33.79, 34.26]
    lowVal = 0.45
    highVal = 0.65
    targetDiff = 29.4

    #0.25MHz, ch1
    x_data = [0.05, 0.1, 0.15000000000000002, 0.2, 0.25, 0.30000000000000004, 0.35000000000000003, 0.4, 0.45, 0.5, 0.55, 0.6000000000000001, 0.65, 0.7000000000000001, 0.75, 0.8, 0.8500000000000001, 0.9, 0.9500000000000001]
    y_data = [8.84, 14.8, 18.27, 20.73, 22.65, 24.22, 25.55, 26.7, 27.72, 28.63, 29.45, 30.2, 30.89, 31.54, 32.13, 32.69, 33.22, 33.71, 34.18]
    lowVal = 0.55
    highVal = 0.75
    targetDiff = 32.1

    #0.25MHz, ch4
    x_data = [0.05, 0.1, 0.15000000000000002, 0.2, 0.25, 0.30000000000000004, 0.35000000000000003, 0.4, 0.45, 0.5, 0.55, 0.6000000000000001, 0.65, 0.7000000000000001, 0.75, 0.8, 0.8500000000000001, 0.9, 0.9500000000000001]
    y_data = [9.27, 14.93, 18.36, 20.82, 22.73, 24.3, 25.63, 26.78, 27.79, 28.7, 29.53, 30.28, 30.97, 31.61, 32.21, 32.77, 33.29, 33.79, 34.26]
    lowVal = 0.03
    highVal = 0.16
    targetDiff = 12.4

    #0.5MHz, ch1
    x_data = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19]
    y_data = [2.25, 2.25, 4.68, 6.97, 8.8, 10.32, 11.62, 12.75, 13.76, 14.66, 15.48, 16.23, 16.92, 17.56, 18.16, 18.71, 19.24, 19.73, 20.2]
    lowVal = 0.1
    highVal = 0.19
    targetDiff = 16.6

    #0.5MHz, ch4
    x_data = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19]
    y_data = [2.22, 2.85, 4.56, 6.38, 8.36, 10.0, 11.37, 12.55, 13.6, 14.53, 15.37, 16.13, 16.84, 17.49, 18.1, 18.66, 19.2, 19.7, 20.17]
    lowVal = 0.02
    highVal = 0.065
    targetDiff = 7.1

    #0.75MHz, ch1
    x_data = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19]
    y_data = [1.78, 2.44, 5.06, 7.21, 8.97, 10.5, 11.8, 12.92, 13.92, 14.81, 15.62, 16.36, 17.04, 17.68, 18.27, 18.82, 19.34, 19.83, 20.29]
    lowVal = 0.02
    highVal = 0.065
    targetDiff = 8.8

    #0.75MHz, ch4
    x_data = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19]
    y_data = [1.66, 3.06, 4.87, 6.94, 8.81, 10.37, 11.69, 12.84, 13.86, 14.77, 15.59, 16.34, 17.03, 17.67, 18.26, 18.82, 19.34, 19.84, 20.3]
    lowVal = 0.005
    highVal = 0.035
    targetDiff = 1.5

    y_data_err = [0.5 for x in y_data]

    #lowVal = np.min(x_data)
    #highVal = np.max(x_data)
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
      X_plotFit = np.linspace(lowVal,highVal,1000)
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
      print("PRED AMP", (targetDiff - intercept)/slope )
    
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
