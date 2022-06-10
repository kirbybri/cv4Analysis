import h5py
import numpy as np
from math import *
import json
import sys
import pickle
import os.path
import struct

class CV4_PROCESS_FILE(object):

  #__INIT__#
  def __init__(self,fileName=None):
    self.fileName = fileName
    self.hdf5File = None
    self.runResultsDict = None
    #self.weights = [32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1] #generic binary weights for 16-bit words
    self.weights = [ np.power(2,num) for num in range(15,-1,-1) ]
    self.maxNumSamples = 10000
    self.limitNumSamples = False
    self.runNo = 0

  #analysis process applied to each "measurement"
  def processMeasurement(self, meas):
    #loop over measurement group members, process waveform data
    resultsDict = {}
    for group_key in meas.keys() :
      mysubgroup = meas[group_key] #group object
      for subgroup_key in mysubgroup.keys() :        
        channelNum = subgroup_key
        samples = mysubgroup[subgroup_key][()] #dataset objectraw_data_array = raw_data[()]
        resultsDict[channelNum] = {'wf' : samples }
        #resultsList.append( resultsDict )
    return resultsDict

  #extract run # from file name
  def getRunNo(self,fileName):
    pathSplit = fileName.split('/')
    nameSplit = pathSplit[-1].split('_')
    if len(nameSplit) < 2 :
      return None
    if nameSplit[0] != 'Run' :
      return None
    return nameSplit[1]

  def getReqAttrs(self,meas=None):
    if meas == None :
      return None
    measAttrs = {}
    for attr in meas.attrs :
      #print(attr,"\t",meas.attrs[attr])
      measAttrs[attr] = meas.attrs[attr]
    return measAttrs

  #run analysis on input file
  def processFile(self):
    self.openFile()
    if self.hdf5File == None :
      print("ERROR: file not open")
      return None

    #define run results dict
    self.runResultsDict = {}
    self.runResultsDict['file'] = self.hdf5File.filename
    self.runResultsDict['run'] = self.runNo

    print("CV4PROCESSFILE : Process File")
    #loop through measurements, store results in dict
    measResultsDict = {}
    for measNumVal, measNum in enumerate(self.hdf5File.keys() ):
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

    #measurement results stored in dict
    self.runResultsDict['results'] = measResultsDict.copy()
    self.hdf5File.close()
    return

  #output results dict to json file
  def outputFile(self):
    if self.runResultsDict == None:
      return None
    pathSplit = self.fileName.split('/')[-1]
    #jsonFileName = 'output_cv3tbProcessFile_' + pathSplit + '.json'
    #with open( jsonFileName , 'w') as outfile:
    #  json.dump( self.runResultsDict, outfile, indent=4)
    pickleFileName = 'output_cv4ProcessFile_' + pathSplit + '.pickle'
    print("Output file ",pickleFileName)
    pickle.dump( self.runResultsDict, open( pickleFileName, "wb" ) )
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
      print("ERROR: couldn't get run number")
      self.runNo = 0
    return None

  #dump file structure to terminal
  def dumpFile(self):
    self.openFile()
    if self.hdf5File == None :
      print("ERROR: file not open")
      return None
    #loop through measurements, store results in dict
    measResultsDict = {}
    print("All HDF5 keys")
    print("\t",self.hdf5File.keys())
    for measNumVal, measNum in enumerate(self.hdf5File.keys() ):
      meas = self.hdf5File[measNum]
      #measAttrs = self.getReqAttrs(meas)
      print("MEAS",measNum)
      print("\tMEAS KEYS = ", meas.keys() )
      print("\tMEAS ATTRS = ", meas.attrs.keys() )
      for attr_key in meas.attrs.keys() :
        myattr = meas.attrs[attr_key]
        print("\tATTR",attr_key,"\t",myattr)
      for group_key in meas.keys() :
        mysubgroup = meas[group_key] #group object
        print("\tGROUP",group_key)
        print("\t\tGROUP KEYS = ", mysubgroup.keys() )
        print("\t\tGROUP ATTRS = ",mysubgroup.attrs.keys())
        for subgroup_key in mysubgroup.keys() :
          mysubsubgroup = mysubgroup[subgroup_key] #group object
          print("\t\t\tSUBGROUP",subgroup_key)
          if group_key == "Calibration" :
            print("\t\t\t\tSUBGROUP KEYS",mysubsubgroup.keys())
            print("\t\t\t\tSUBGROUP ATTRS",mysubsubgroup.attrs.keys())
          print("\t\t\t\tSUBSUBGROUP", mysubsubgroup )
          #for num,element in enumerate(mysubsubgroup) :
          #  binElement = '{0:032b}'.format(element)
          #  print( "\t\t\t\t\t",num,element,binElement )
          #  if num > 100 : break
          #continue
          if len(mysubsubgroup) == 1 :
            print("\t\t\t\t\tMEMBER",mysubsubgroup )
          else :
            for member in mysubsubgroup :
              print("\t\t\t\t\tMEMBER",member )
              if group_key == "Calibration" :
                print( len(member) )
                continue
              else:
                #for num,element in enumerate(member) :
                #  print("\t\t\t\t\t\t",num,"\t",element)
                #  if num > 10 : break
                break
            #print("\t\t\t\t\tMEMBER",mysubsubgroup[0] )
          #print("\t\t\t\t\tMEMBER", len(mysubsubgroup))
          #for member in mysubsubgroup :
          #  print("\t\t\t\t\tMEMBER", member )
          colutaNum = group_key
          channelNum = subgroup_key

  def outputTextFiles(self):
    self.openFile()
    if self.hdf5File == None :
      print("ERROR: file not open")
      return None
    #loop through measurements, store results in dict
    measResultsDict = {}
    for measNumVal, measNum in enumerate(self.hdf5File.keys() ):
      meas = self.hdf5File[measNum]
      for group_key in meas.keys() :
        mysubgroup = meas[group_key] #group object
        if group_key != "COLUTA" : continue
        numFiles = 0 
        for subgroup_key in mysubgroup.keys() :
          mysubsubgroup = mysubgroup[subgroup_key] #group object
          chList = ["channel1","channel2","channel3","channel4","channel5","channel6","channel7","channel8"]
          if subgroup_key not in chList : continue
          samples = mysubsubgroup
          pathSplit = self.fileName.split('/')[-1]
          outputFileName = "data/output_cv4ProcessFile_wfData_" + str(pathSplit) + "_" + str(measNum) + ".txt"
          print("Writing",outputFileName)
          f = open(outputFileName, "w")
          for num,samp in enumerate(samples) :
            #if num == 0 : continue
            #if num > 100 : break
            binsamp = '{0:032b}'.format(samp)
            f.write( str(binsamp) + "\n")
          f.close()
          #numFiles = numFiles + 1
          #if numFiles > 100:
          #  break
          #end element/sample loop
        #end subgroup 


def main():
  print("HELLO")
  if len(sys.argv) != 2 :
    print("ERROR, program requires filename as argument")
    return
  fileName = sys.argv[1]
  processFile = CV4_PROCESS_FILE(fileName)
  processFile.openFile()
  processFile.dumpFile()
  #processFile.processFile()
  #processFile.outputFile()

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
