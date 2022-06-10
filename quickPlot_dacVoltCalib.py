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

    calibData = {}

    data = {}
    with open('outputData_doDacVoltTest_CV4_154_DACp_ch1Disconn.json', 'r') as fp:
      temp = json.load(fp)  
      data = temp
    if "CV4_DACp" not in data :
      print("NO CV4_DACp")
      return
    calibData["CV4_DACp"] = data["CV4_DACp"]

    data = {}
    with open('outputData_doDacVoltTest_CV4_154_DACn_ch1Disconn.json', 'r') as fp:
      temp = json.load(fp)  
      data = temp
    if "CV4_DACn" not in data :
      print("NO CV4_DACn")
      return
    calibData["CV4_DACn"] = data["CV4_DACn"]
    
    x_data = []
    y_data = []
    x_data_err = []
    y_data_err = []
    for dacVal in calibData["CV4_DACp"] :
      if dacVal not in calibData["CV4_DACn"] :
        print("MISSING DAC VAL ",dacVal)
        return
      volt =  float(calibData["CV4_DACp"][dacVal] ) - float(calibData["CV4_DACn"][dacVal])
      print(dacVal,"\t",calibData["CV4_DACp"][dacVal],"\t",calibData["CV4_DACn"][dacVal] ,"\t",volt )
      x_data.append( int(dacVal) )
      y_data.append( float(volt)*1000 ) #mV
      x_data_err.append( 0 )
      y_data_err.append( 0.00004*1000 )
    #lowVal = 0
    #highVal = 65536
    lowVal = 200
    highVal = 65400
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
    
    fig, axes = plt.subplots(2,1,figsize=(10, 9))
    axes[0].plot(x_data,y_data,".")
    axes[0].plot(X_plotFit,Y_plotFit,"-")
    axes[0].set_xlabel('DAC Code', horizontalalignment='right', x=1.0, fontsize=16)
    axes[0].set_ylabel('DAC Voltage [mV,diff]', horizontalalignment='center', x=1.0, fontsize=16)
    axes[0].set_title("DAC Voltage vs DAC Code", fontsize=16)
    axes[0].tick_params(axis="x", labelsize=16)
    axes[0].tick_params(axis="y", labelsize=16)
    #axes[0].set_xlim(0,2400)
    axes[0].grid()
    axes[0].text(0.05, 0.95, textstr, transform=axes[0].transAxes, fontsize=14, verticalalignment='top')


    axes[1].plot(resid_x,resid_y,".")
    axes[1].set_xlabel('DAC Code', horizontalalignment='right', x=1.0, fontsize=16)
    axes[1].set_ylabel('Fit Residual [mV,diff]', horizontalalignment='center', x=1.0, fontsize=16)
    #axes[1][0].set_xlim(0,2400)
    #axes[1][0].set_ylim(-25,25)
    axes[1].tick_params(axis="x", labelsize=16)
    axes[1].tick_params(axis="y", labelsize=16)
    axes[1].grid()
    axes[1].set_title("DAC Voltage Fit Residuals vs DAC Code", fontsize=16)

    fig.suptitle("Board 153, Ch1+2 Connected", fontsize=16)
    fig.tight_layout()
    plt.show()
    
    

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
