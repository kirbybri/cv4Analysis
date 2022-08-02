import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
from fpdf import FPDF
from pathlib import Path
import os.path

class CV4_PROCESS_FILE(object):
  #__INIT__#
  def __init__(self,fileName=None):
    self.fileName = fileName
    self.hdf5File = None
    self.runResultsDict = None
    self.runNo = 0
    self.runNameBase = "RadRun"

  def processMeasurement(self, meas):
    #loop over measurement group members, process waveform data
    resultsDict = {}
    for group_key in meas.keys() :
      mysubgroup = meas[group_key] #group object
      if group_key == "COLUTA" :
        for subgroup_key in mysubgroup.keys() :        
          channelNum = subgroup_key
          samples = mysubgroup[subgroup_key][()]
          resultsDict[channelNum] = {'wf' : samples }
      if group_key == "Histogram" :
        for subgroup_key in mysubgroup.keys() :
          channelNum = subgroup_key
          binCount  = mysubgroup[channelNum]['bin_count'][()]
          binNumber = mysubgroup[channelNum]['bin_number'][()]
          resultsDict[channelNum] = {'binCount' : binCount , "binNumber" : binNumber }
    return resultsDict
  
  def getRunNo(self,fileName): #extract run # from file name
    fileNamePath = Path(fileName)
    if len( fileNamePath.parts  ) == 0 :
      print("INVALID filename")
      return None
    fileName = fileNamePath.parts[-1]
    nameSplit = fileName.split('_')
    if len(nameSplit) < 2 :
      return None
    if nameSplit[0] != self.runNameBase :
      return None
    return nameSplit[1]

  def getReqAttrs(self,meas=None):
    if meas == None :
      return None
    measAttrs = {}
    for attr in meas.attrs :
      measAttrs[attr] = meas.attrs[attr]
    return measAttrs
    
  def getFileAttrs(self):
    if self.hdf5File == None :
      return None
    fileAttrs = {}
    for attr in self.hdf5File.attrs.keys() :
      fileAttrs[attr] = self.hdf5File.attrs[attr]
    return fileAttrs

  def processFile(self):
    self.openFile()
    if self.hdf5File == None :
      print("ERROR: file not open")
      return None

    self.runResultsDict = {}
    self.runResultsDict['file'] = self.hdf5File.filename
    self.runResultsDict['run'] = self.runNo

    print("CV4PROCESSFILE : Process file")
    measResultsDict = {}
    fileAttrs = self.getFileAttrs()
    for measNumVal, measNum in enumerate( self.hdf5File.keys() ):
      meas = self.hdf5File[measNum]
      measAttrs = self.getReqAttrs(meas)
      measVal = measNum.split("_")
      if len(measVal ) != 2 : continue
      measVal = int(measVal[1])
      if measNumVal % 100 == 0 :
        print( "Measurement","\t",measNum,"\t",measNumVal)        
      measData = self.processMeasurement(meas)
      if measData == None :
        print("Missing waveform data, will not process measurement")
        continue
      measResultsDict[measNum] = {'data':measData.copy(),'attrs':measAttrs.copy()}
    print("CV4PROCESSFILE : Process file done")
    self.runResultsDict['results'] = measResultsDict.copy()
    self.runResultsDict['fileAttrs'] = fileAttrs
    self.hdf5File.close()
    return

  def openFile(self):
    if self.fileName == None:
      print("ERROR: no file name supplied")
      return None
    if os.path.isfile(self.fileName) == False :
      print("ERROR: file does not exist ",self.fileName)
      return None
    self.hdf5File = None
    try:
      self.hdf5File = h5py.File(self.fileName, "r") #file object
    except:
      print("ERROR: couldn't open file",self.fileName)
      return None
    #check for required attributes
    self.runNo = self.getRunNo( self.hdf5File.filename )
    if self.runNo == None :
      print("ERROR: couldn't get run number, assigning run number 0")
      self.runNo = 0
    return None
#END CV4_PROCESS_FILE CLASS
    
class CV4_PROCESS_RADTESTSINE(object):
  #__INIT__#
  def __init__(self, fileName = None):
    self.fileName = fileName
    self.runResultsDict = None
    self.gotResults_sine = False
    self.gotResults_table = False
    self.fileNameBase = "default"
    self.plotFileName = 'sinePlot' + self.fileNameBase + '.png'
    self.tableFileName = 'sineTable' + self.fileNameBase + '.png'
    self.measResults = {}
    
  def SINAD(self,fourier):
    sum2 = 0
    for normBin in fourier:
      if normBin==1: continue
      sum2 += normBin**2
    return -10*np.log10(sum2)

  def ENOB(self,fourier):
    return (self.SINAD(fourier)-1.76)/6.02

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
      return psd_x,psd,0,0 #error, should not happen

    fourier_fftWf_y = fftWf_y/np.max(fftWf_y)
    for sampNum,samp in enumerate(fourier_fftWf_y) :
      if sampNum == 0 :
        continue
      else:
        psdVal = 0
        if samp > 0 :
          psdVal = 20*np.log10(samp)
        psd_x.append( fftWf_x[sampNum] )
        psd.append( psdVal )
    sinad = self.SINAD(fourier_fftWf_y)
    enob = self.ENOB(fourier_fftWf_y)
    return psd_x,psd,sinad,enob

  def getMeasData(self,measNum=None,measInfo=None):
    if measNum == None or measInfo == None:
      print("ERROR getMeasData, invalid input")
      return None
    if "data" not in measInfo:
      print("ERROR getMeasData, measurement does not contain data,",measNum)
      return None
    measData = measInfo["data"]
    if "attrs" not in measInfo:
      print("ERROR getMeasData, measurement does not contain metadata,",measNum)
      return None
    measAttrs = measInfo["attrs"]
    return measData,measAttrs

  def makeTable(self):
    reqChannels = ["channel1","channel3","channel5","channel7"]
    reqAmps = [0.5,1.5,2.5,3.15]
    #define table
    data = []
    row0 = ["0.5V","1.5V","2.5V","3.15V"]
    data.append(row0)
    for chan in reqChannels :
      rowData = []
      rowData.append(chan)
      for amp in reqAmps :
        meanVal = None
        stdVal = None
        if chan in self.measResults :
          if amp in self.measResults[chan] :
            meanVal = round(self.measResults[chan][amp]["meanEnob"],3)
            stdVal = round(self.measResults[chan][amp]["stdEnob"],3)
        if meanVal != None : 
          rowData.append( str(meanVal)+" +- "+str(stdVal) )
        else :
          rowData.append( "" )
      data.append(rowData)
    # Pop the headers from the data array
    column_headers = data.pop(0)
    row_headers = [x.pop(0) for x in data]
    
    cell_text = []
    for row in data:
      cell_text.append(row)
      
    title_text = 'FFT Results'
    plt.figure(linewidth=2, tight_layout={'pad':1}, figsize=(5,2) )
    the_table = plt.table(cellText=cell_text, rowLabels=row_headers, rowLoc='right', colLabels=column_headers, loc='center')                 
    # Hide axes
    ax = plt.gca()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    # Hide axes border
    plt.box(on=None)
    self.tableFileName = 'sineTable' + self.fileNameBase + '.png'
    plt.savefig(self.tableFileName, bbox_inches='tight', dpi=150 )
    self.gotResults_table = True

  def makePlot(self):
    #dimensions of summary plots 
    numRows = 2
    numCols = 4

    #dict of required results, organized by plot panel
    reqPlotDict = {}
    reqPlotDict[(0,0)] = {'ch':"channel1",'amp':3.15,'data':"wf","title":"CH1,1Vpp"}
    reqPlotDict[(0,1)] = {'ch':"channel3",'amp':3.15,'data':"wf","title":"CH3,1Vpp"}
    reqPlotDict[(0,2)] = {'ch':"channel5",'amp':3.15,'data':"wf","title":"CH5,1Vpp"}
    reqPlotDict[(0,3)] = {'ch':"channel7",'amp':3.15,'data':"wf","title":"CH7,1Vpp"}
    reqPlotDict[(1,0)] = {'ch':"channel1",'amp':3.15,'data':"psd_y","title":"CH1,1Vpp PSD"}
    reqPlotDict[(1,1)] = {'ch':"channel3",'amp':3.15,'data':"psd_y","title":"CH3,1Vpp PSD"}
    reqPlotDict[(1,2)] = {'ch':"channel5",'amp':3.15,'data':"psd_y","title":"CH5,1Vpp PSD"}
    reqPlotDict[(1,3)] = {'ch':"channel7",'amp':3.15,'data':"psd_y","title":"CH7,1Vpp PSD"}
              
    #plot required data, leave pad empty if data missing
    fig, axes = plt.subplots(numRows,numCols,figsize=(14, 6))
    for row in range(0,numRows,1):
      for col in range(0,numCols,1):
        reqData = reqPlotDict[(row,col)]
        vals_x = []
        vals_y = []
        if reqData['ch'] in self.measResults :
          if reqData['amp'] in self.measResults[ reqData['ch'] ] :
            if len( self.measResults[ reqData['ch'] ][ reqData['amp'] ]["measResults"]  ) > 0 :
              if reqData['data'] in self.measResults[ reqData['ch'] ][ reqData['amp'] ]["measResults"][ 0 ] :
                if reqData['data'] == "wf" :
                  vals_y = self.measResults[ reqData['ch'] ][ reqData['amp'] ]["measResults"][ 0 ][ reqData['data'] ]
                  vals_x = [x*25 for x in range(0,len(vals_y),1)]
                if reqData['data'] == "psd_y" :
                  vals_x = self.measResults[ reqData['ch'] ][ reqData['amp'] ]["measResults"][ 0 ][ "psd_x" ]
                  vals_y = self.measResults[ reqData['ch'] ][ reqData['amp'] ]["measResults"][ 0 ][ "psd_y" ]
        axes[row][col].plot(vals_x,vals_y,".")
        if reqData['data'] == "wf" :
          axes[row][col].set_xlabel('Sample Time [ns]', horizontalalignment='right', x=1.0)
          axes[row][col].set_ylabel('ADC CODE [ADC]', horizontalalignment='left', x=1.0)
        if reqData['data'] == "psd_y" :
          axes[row][col].set_xlabel('Frequency [MHz]', horizontalalignment='right', x=1.0)
          axes[row][col].set_ylabel('PSD [dB]', horizontalalignment='left', x=1.0)
        axes[row][col].set_title( reqData['title'] )

    plotTitle = "Summary: " + str(self.fileNameBase)
    fig.suptitle(plotTitle, fontsize=16)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    self.gotResults_sine = True
    self.plotFileName = 'sinePlot' + self.fileNameBase + '.png'
    plt.savefig(self.plotFileName, bbox_inches='tight', dpi=150 )
    return

  def processMeasurement(self,measNum=None,measInfo=None):
    #get measurement data objects
    measDataObject = self.getMeasData(measNum,measInfo)
    if measDataObject == None:
      print("Could not parse measurement data objects,",measNum)
      return None
    measData,measAttrs = measDataObject
    
    #first check if sine wave measurement
    if "measType" not in measAttrs :
      print("ERROR processMeasurement, required field measType missing from metadata,",measNum)
      return None
    measType = measAttrs["measType"]
    if measType != "FFT" :
      return None
    #have sine wave measurement here, check for other attributes and get performance measurements
    reqAttrs = ["measTime",'measChan','awgAmp','awgFreq']
    for attr in reqAttrs :
      if attr not in measAttrs :
        print("ERROR processMeasurement, required field missing from metadata,",measNum,attr)
        return None
    chName = measAttrs["measChan"]
    amp = measAttrs['awgAmp']
    if chName not in measData :
      print("ERROR processMeasurement, required channel data missing from measurement")
      return None
    chData = measData[chName]
    if "wf" not in chData :
      print("ERROR processMeasurement, required waveform data missing from measurement")
      return None
    vals = chData["wf"]
    result = self.getFftWaveform( vals )
    if result == None :
      print("ERROR processMeasurement, could not get FFT")
      return None
    psd_x,psd,sinad,enob = result
    #save performance measurements
    if chName not in self.measResults:
      self.measResults[chName] = {}
    if amp not in self.measResults[chName] :
      self.measResults[chName][amp] = {"measResults":[]}
    self.measResults[chName][amp]["measResults"].append( {"wf":vals,"psd_x":psd_x,"psd_y":psd,"enob":enob} )
    return None
  
  def getAvgResults(self):
    if self.measResults == None :
      return None
    for chName in self.measResults :
      for amp in self.measResults[chName] :
        enobs = []
        for measNum,meas in enumerate(self.measResults[chName][amp]["measResults"]) :
          enobs.append( self.measResults[chName][amp]["measResults"][measNum]["enob"] )
        meanEnob = np.mean( enobs )
        stdEnob = np.std( enobs )
        self.measResults[chName][amp]["meanEnob"] = meanEnob
        self.measResults[chName][amp]["stdEnob"] = stdEnob
    return None

  def getFileNameBase(self):
    self.fileNameBase = "default"
    if self.fileName == None :
      print("NO filename")
      return None
    fileNamePath = Path(self.fileName)
    if len( fileNamePath.parts  ) == 0 :
      print("INVALID filename")
      return None
    self.fileNameBase = fileNamePath.parts[-1]    
    return None

  def makeSummary(self):
    self.getAvgResults()
    self.makePlot()
    self.makeTable()

    pdf = FPDF(format='letter')
    pdf.add_page()
    pdf.set_font("helvetica", 'B', size=20)
    pdf.set_fill_color(255,255,255)
    pdf.rect(0,0,pdf.w,pdf.h,"F") #fill page background white

    pdf.cell(200, 10, txt="RadBoard Test Summary", align='C',new_x="LMARGIN",new_y="NEXT")
    pdf.set_font("helvetica", 'B', size=16)
    pdf.ln(4)  
    pdf.cell(60, 5, txt="File name: " + self.fileNameBase , align='L',new_x="LMARGIN",new_y="NEXT")
    pdf.ln(4)
    pdf.cell(60, 5, txt="Board: " + str(self.measBoard) , align='L',new_x="LMARGIN",new_y="NEXT")
    pdf.ln(4)
    pdf.cell(60, 5, txt="ASIC: " + str(self.measASIC) , align='L',new_x="LMARGIN",new_y="NEXT")
    pdf.ln(4)
  
    text = "Sine Wave Test Waveforms and PSDs"
    pdf.cell(200,5,txt=text,align='L',new_x="LMARGIN",new_y="NEXT")
    if self.gotResults_sine == True :
        pdf.image(self.plotFileName, w=180)
    pdf.ln(4)
    text = "Sine Wave Test Average ENOB"
    pdf.cell(200,5,txt=text,align='L',new_x="LMARGIN",new_y="NEXT")
    if self.gotResults_table == True :
        pdf.image(self.tableFileName, w=180)
    pdf.ln(4)
    summaryFileName = 'summary_' + self.fileNameBase + '.pdf'
    pdf.output(summaryFileName)
    print("PRODUCED SUMMARY FILE",summaryFileName)
    return None
  
  def getFileAttrs(self):
    if self.runResultsDict == None :
      print("ERROR, no data recovered from file ,exiting")
      return None
    if "fileAttrs" not in self.runResultsDict:
      print("ERROR, missing file attributes,exiting")
      return None
    if "measBoard" not in self.runResultsDict["fileAttrs"] or "measASIC" not in self.runResultsDict["fileAttrs"] :
      print("ERROR, missing required file attributes,exiting")
      return None
    self.measBoard = self.runResultsDict["fileAttrs"]["measBoard"]
    self.measASIC = self.runResultsDict["fileAttrs"]["measASIC"]
    return None
    
  def processFileData(self):
    if self.runResultsDict == None :
      print("ERROR, no data recovered from file ,exiting")
      return None
    if "results" not in self.runResultsDict:
      print("ERROR, no data recovered from file ,exiting")
      return None
    self.getFileNameBase()
    self.getFileAttrs()
    for cnt, (measNum, measInfo) in enumerate( self.runResultsDict["results"].items() ):
      self.processMeasurement(measNum,measInfo)
    #self.getAvgResults()
    #self.makePlot()
    #self.makeTable()
    #self.makeSummary()
    return None
#END CV4_PROCESS_RADTESTSINE CLASS

def main():
  print("HELLO, running cv4tbProcessRadTestAnalogPerfData")

  if len(sys.argv) < 2 :
    print("ERROR, cv4tbProcessRadTestAnalogPerfData requires filename(s)")
    return
  fileNames = sys.argv[1:]

  print("PROCESSING FILES",fileNames)
  cv4ProcessFile = CV4_PROCESS_FILE(fileNames[0])
  cv4ProcessRadTestSine = CV4_PROCESS_RADTESTSINE(fileNames[0])
 
  for fileName in fileNames :
    print("PROCESS FILE",fileName)
    cv4ProcessFile.fileName = fileName
    cv4ProcessFile.processFile()
    if cv4ProcessFile.runResultsDict == None:
      print("ERROR, could not process ",fileName)
      return None
    print("DONE PROCESS FILE")
    print("PROCESS SINE WAVE TEST")
    cv4ProcessRadTestSine.runResultsDict = cv4ProcessFile.runResultsDict
    cv4ProcessRadTestSine.processFileData()
    print("DONE PROCESS SINE WAVE TEST")
  cv4ProcessRadTestSine.makeSummary()
  return

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
