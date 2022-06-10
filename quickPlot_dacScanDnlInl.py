import numpy as np
from math import *
import matplotlib.pyplot as plt
import json

def main():


    data = {}
    with open('outputData_CV4_dacScanHist.json', 'r') as fp:
      temp = json.load(fp)  
      data = temp
    
    minCode = 0+1000
    maxCode = 32767-1000
    
    avgVal = 0
    avgValSimple = 0
    numBins = 0
    for code in data :
      codeVal = int (code)
      if codeVal < minCode : continue
      if codeVal > maxCode : continue
      count = data[code] 
      avgVal = avgVal + (count - avgVal)/float(numBins+1)
      numBins = numBins + 1
      #print(code,count)
    
    minCode = 0
    maxCode = 32767
    dnl_x = []
    dnl_y = []
    inl_y = []
    inlSum = 0
    for code in data :
      codeVal = int (code)
      if codeVal < minCode : continue
      if codeVal > maxCode : continue
      count = data[code] 
      
      codeDnl = (count - avgVal)/avgVal
      
      dnl_x.append(codeVal)
      dnl_y.append(codeDnl)
      inlSum = inlSum + codeDnl
      inl_y.append(inlSum)
      
      if codeDnl <= -0.99 : 
        print( code, count , codeDnl )  
    


    """
    fig_dnl, axes_dnl = plt.subplots(1,1,figsize=(10, 6))
    axes_dnl.plot(dnl_x,dnl_y,marker='o',linestyle='')
    axes_dnl.set_xlabel('COLUTA Code [ADC]', fontsize=20)
    axes_dnl.set_ylabel('DNL [15-Bit LSB]', fontsize=20)
    axes_dnl.tick_params(axis="x", labelsize=12)
    axes_dnl.tick_params(axis="y", labelsize=12)
    axes_dnl.set_title("Run 1244 Normal Mode DAC Scan DNL vs COLUTA Code Value" )
    fig_dnl.tight_layout()
    """
    fig_dnlHist, axes_dnlHist = plt.subplots(1,1,figsize=(10, 6))
    axes_dnlHist.hist(dnl_y,bins = np.arange(min(dnl_y), max(dnl_y)+1, 0.01))
    axes_dnlHist.set_xlabel('DNL [15-Bit LSB]', fontsize=20)
    axes_dnlHist.set_ylabel('Number of Codes', fontsize=20)
    axes_dnlHist.set_title("Run 1244 Normal DAC Scan DNL Distribution" )
    #axes_dnlHist.text(0.85, 0.75, textstr_dnl, transform=axes_dnl.transAxes, fontsize=14, verticalalignment='top')
    axes_dnlHist.tick_params(axis="x", labelsize=12)
    axes_dnlHist.tick_params(axis="y", labelsize=12)
    axes_dnlHist.set_yscale('log')
    axes_dnlHist.set_xlim(-1.1,1.1)
    axes_dnlHist.grid()
    fig_dnlHist.tight_layout()
    """
    fig_inl, axes_inl = plt.subplots(1,1,figsize=(10, 6))
    axes_inl.plot(dnl_x,inl_y,marker='o',linestyle='')
    axes_inl.set_xlabel('COLUTA Code [ADC]', fontsize=20)
    axes_inl.set_ylabel('INL [15-bit LSB]', fontsize=20)
    axes_inl.tick_params(axis="x", labelsize=12)
    axes_inl.tick_params(axis="y", labelsize=12)
    axes_inl.set_title("Run 1244 Normal Mode DAC Scan INL vs COLUTA Code Value" )
    fig_inl.tight_layout()    
    """
    
    plt.show()     
    return
    
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
