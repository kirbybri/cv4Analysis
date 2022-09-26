import numpy as np
from math import *
import matplotlib.pyplot as plt

def calcPsrr(dataOn,dataOff,dIn) :
    result = {}
    for chanName in dataOn :
      if chanName not in dataOff : continue
      fftOn = dataOn[chanName]["maxPsd"]
      fftOff = dataOff[chanName]["maxPsd"]
      diff = fftOn - fftOff
      rmsOn  = fftOn/8000.*sqrt(2)
      rmsOff = fftOff/8000.*sqrt(2)
      rmsDiff = diff/8000.*sqrt(2)
      #print(chanName, rmsOn , rmsOff , rmsDiff)
      if rmsOn > rmsOff :
        disturb = sqrt(rmsOn*rmsOn - rmsOff*rmsOff)
      else :
        disturb = 0
      if disturb == 0 :
        print(chanName, rmsOn , rmsOff , 0)
        continue
      convFactor = 2000. / 32767.
      disturbAmp = 2*sqrt(2)*disturb
      disturbAmpDiff = 2*sqrt(2)*rmsDiff
      psrr = 20*log10( dIn / ( disturbAmp * convFactor) )
      psrrDiff = 20*log10( dIn / ( disturbAmpDiff * convFactor) )
      #print(chanName, rmsOn , rmsOff , rmsDiff, disturb , psrr)
      print(chanName, psrr , psrrDiff)
      if psrr > 95 : psrr = 95
      if psrrDiff > 95 : psrrDiff = 95

      #result[chanName] = psrr
      result[chanName] = psrrDiff
    return result

def main():

    #115kHz no disturbance 
    dIn_115 = 16.8
    result_115_off = {'channel1': {'maxPsd': 125.26925841206713, 'freq': 0.115}, 'channel2': {'maxPsd': 120.50189343841151, 'freq': 0.115}, 'channel3': {'maxPsd': 120.46701410366452, 'freq': 0.115}, 'channel4': {'maxPsd': 137.57756619061348, 'freq': 0.115}, 'channel5': {'maxPsd': 140.93093197214904, 'freq': 0.115}, 'channel6': {'maxPsd': 116.48244357343854, 'freq': 0.115}, 'channel7': {'maxPsd': 129.4334606795059, 'freq': 0.115}, 'channel8': {'maxPsd': 123.581046940801, 'freq': 0.115}}
    result_115_on = {'channel1': {'maxPsd': 1063.4261476095169, 'freq': 0.115}, 'channel2': {'maxPsd': 794.3708030725462, 'freq': 0.115}, 'channel3': {'maxPsd': 1451.4198701167418, 'freq': 0.115}, 'channel4': {'maxPsd': 2188.5233363501247, 'freq': 0.115}, 'channel5': {'maxPsd': 1273.7770041746728, 'freq': 0.115}, 'channel6': {'maxPsd': 1123.5164100745599, 'freq': 0.115}, 'channel7': {'maxPsd': 1332.7127869568544, 'freq': 0.115}, 'channel8': {'maxPsd': 119.96647388686588, 'freq': 0.115}}

    #305kHz
    dIn_305 = 33.
    result_305_off = {'channel1': {'maxPsd': 101.31604116128933, 'freq': 0.305}, 'channel2': {'maxPsd': 107.25622994487266, 'freq': 0.305}, 'channel3': {'maxPsd': 101.1647514689482, 'freq': 0.305}, 'channel4': {'maxPsd': 108.71931846334114, 'freq': 0.305}, 'channel5': {'maxPsd': 106.27730113686034, 'freq': 0.305}, 'channel6': {'maxPsd': 102.49084154925612, 'freq': 0.305}, 'channel7': {'maxPsd': 111.96914698182604, 'freq': 0.305}, 'channel8': {'maxPsd': 110.55520643875128, 'freq': 0.305}}
    result_305_on = {'channel1': {'maxPsd': 1953.0504284436847, 'freq': 0.305}, 'channel2': {'maxPsd': 718.768322367164, 'freq': 0.305}, 'channel3': {'maxPsd': 2375.824504493591, 'freq': 0.305}, 'channel4': {'maxPsd': 4569.856576430726, 'freq': 0.305}, 'channel5': {'maxPsd': 2571.8693648854246, 'freq': 0.305}, 'channel6': {'maxPsd': 1761.2966165819673, 'freq': 0.305}, 'channel7': {'maxPsd': 2232.844023700281, 'freq': 0.305}, 'channel8': {'maxPsd': 595.4527570875981, 'freq': 0.305}}

    #505kHz
    dIn_505 = 14.
    result_505_off = {'channel1': {'maxPsd': 97.53627075190172, 'freq': 0.505}, 'channel2': {'maxPsd': 98.09336824722298, 'freq': 0.505}, 'channel3': {'maxPsd': 100.4343534435908, 'freq': 0.505}, 'channel4': {'maxPsd': 106.78835367287249, 'freq': 0.505}, 'channel5': {'maxPsd': 107.39627529710042, 'freq': 0.505}, 'channel6': {'maxPsd': 99.50961354826032, 'freq': 0.505}, 'channel7': {'maxPsd': 109.1727101549279, 'freq': 0.505}, 'channel8': {'maxPsd': 97.6847144710903, 'freq': 0.505}}
    result_505_on = {'channel1': {'maxPsd': 1589.37382338563, 'freq': 0.505}, 'channel2': {'maxPsd': 680.554290266049, 'freq': 0.505}, 'channel3': {'maxPsd': 1647.6144746775478, 'freq': 0.505}, 'channel4': {'maxPsd': 2678.396890879443, 'freq': 0.505}, 'channel5': {'maxPsd': 1163.930528488098, 'freq': 0.505}, 'channel6': {'maxPsd': 1012.1094369664384, 'freq': 0.505}, 'channel7': {'maxPsd': 1545.8496671931644, 'freq': 0.505}, 'channel8': {'maxPsd': 206.11256755301693, 'freq': 0.505}}

    #755kHz
    dIn_755 = 9.
    result_755_off = {'channel1': {'maxPsd': 93.06433564823917, 'freq': 0.755}, 'channel2': {'maxPsd': 101.39250255169496, 'freq': 0.755}, 'channel3': {'maxPsd': 106.41296898203368, 'freq': 0.755}, 'channel4': {'maxPsd': 106.0953969949402, 'freq': 0.755}, 'channel5': {'maxPsd': 97.56670953522634, 'freq': 0.755}, 'channel6': {'maxPsd': 89.12047979329382, 'freq': 0.755}, 'channel7': {'maxPsd': 100.08366899277989, 'freq': 0.755}, 'channel8': {'maxPsd': 99.47368719183501, 'freq': 0.755}}
    result_755_on = {'channel1': {'maxPsd': 1364.0999188449125, 'freq': 0.755}, 'channel2': {'maxPsd': 819.2179341919353, 'freq': 0.755}, 'channel3': {'maxPsd': 1055.788973012035, 'freq': 0.755}, 'channel4': {'maxPsd': 1814.4685437840878, 'freq': 0.755}, 'channel5': {'maxPsd': 862.0671721532941, 'freq': 0.755}, 'channel6': {'maxPsd': 604.5599254638433, 'freq': 0.755}, 'channel7': {'maxPsd': 1018.2106232643129, 'freq': 0.755}, 'channel8': {'maxPsd': 269.28478193818, 'freq': 0.755}}
    
    #1.005MHz
    dIn_1005 = 10.
    result_1005_off = {'channel1': {'maxPsd': 86.20070664870187, 'freq': 1.005}, 'channel2': {'maxPsd': 94.75596974290139, 'freq': 1.005}, 'channel3': {'maxPsd': 102.56256847418773, 'freq': 1.005}, 'channel4': {'maxPsd': 105.55970019388063, 'freq': 1.005}, 'channel5': {'maxPsd': 98.08940019467387, 'freq': 1.005}, 'channel6': {'maxPsd': 98.2857991245614, 'freq': 1.005}, 'channel7': {'maxPsd': 107.1417818707071, 'freq': 1.005}, 'channel8': {'maxPsd': 103.85472265586301, 'freq': 1.005}}
    result_1005_on = {'channel1': {'maxPsd': 630.1930514514356, 'freq': 1.005}, 'channel2': {'maxPsd': 591.0343765498644, 'freq': 1.005}, 'channel3': {'maxPsd': 670.2497319626289, 'freq': 1.005}, 'channel4': {'maxPsd': 1233.893421143232, 'freq': 1.005}, 'channel5': {'maxPsd': 1200.3657617783185, 'freq': 1.005}, 'channel6': {'maxPsd': 598.7109736913815, 'freq': 1.005}, 'channel7': {'maxPsd': 1075.0134663590566, 'freq': 1.005}, 'channel8': {'maxPsd': 406.9986250363121, 'freq': 1.005}}

    #2MHz
    dIn_2005 = 20.8
    result_2005_off = {'channel1': {'maxPsd': 98.14539555258963, 'freq': 2.005}, 'channel2': {'maxPsd': 105.34466418155492, 'freq': 2.005}, 'channel3': {'maxPsd': 101.85093619718965, 'freq': 2.005}, 'channel4': {'maxPsd': 98.85220163912878, 'freq': 2.005}, 'channel5': {'maxPsd': 105.7046581089303, 'freq': 2.005}, 'channel6': {'maxPsd': 99.77371283469314, 'freq': 2.005}, 'channel7': {'maxPsd': 96.77940913918384, 'freq': 2.005}, 'channel8': {'maxPsd': 95.9517505551115, 'freq': 2.005}}
    result_2005_on = {'channel1': {'maxPsd': 328.9130419872623, 'freq': 2.005}, 'channel2': {'maxPsd': 425.7153207840925, 'freq': 2.005}, 'channel3': {'maxPsd': 202.5700890261862, 'freq': 2.005}, 'channel4': {'maxPsd': 535.5393196364163, 'freq': 2.005}, 'channel5': {'maxPsd': 614.4838856580786, 'freq': 2.005}, 'channel6': {'maxPsd': 183.13638929650438, 'freq': 2.005}, 'channel7': {'maxPsd': 495.85981913850355, 'freq': 2.005}, 'channel8': {'maxPsd': 337.3424215269003, 'freq': 2.005}}

    #3MHz
    dIn_3005 = 19
    result_3005_off = {'channel1': {'maxPsd': 97.7219816881218, 'freq': 3.005}, 'channel2': {'maxPsd': 101.19948306766844, 'freq': 3.005}, 'channel3': {'maxPsd': 94.0812632774159, 'freq': 3.005}, 'channel4': {'maxPsd': 104.39085099334514, 'freq': 3.005}, 'channel5': {'maxPsd': 94.59251625330288, 'freq': 3.005}, 'channel6': {'maxPsd': 91.47769620215402, 'freq': 3.005}, 'channel7': {'maxPsd': 92.86954449773448, 'freq': 3.005}, 'channel8': {'maxPsd': 96.2532244519965, 'freq': 3.005}}
    result_3005_on = {'channel1': {'maxPsd': 191.51604751919993, 'freq': 3.005}, 'channel2': {'maxPsd': 226.13125267028033, 'freq': 3.005}, 'channel3': {'maxPsd': 101.6219399331544, 'freq': 3.005}, 'channel4': {'maxPsd': 206.0024813289132, 'freq': 3.005}, 'channel5': {'maxPsd': 251.02577394674478, 'freq': 3.005}, 'channel6': {'maxPsd': 120.50753340130008, 'freq': 3.005}, 'channel7': {'maxPsd': 191.01556959148675, 'freq': 3.005}, 'channel8': {'maxPsd': 204.167403322671, 'freq': 3.005}}

    #4MHz
    dIn_3985 = 20.6
    result_3985_off = {'channel1': {'maxPsd': 101.07708663812485, 'freq': 3.985}, 'channel2': {'maxPsd': 96.25410282654164, 'freq': 3.985}, 'channel3': {'maxPsd': 100.90829426858329, 'freq': 3.985}, 'channel4': {'maxPsd': 99.91122316856972, 'freq': 3.985}, 'channel5': {'maxPsd': 99.61740416243813, 'freq': 3.985}, 'channel6': {'maxPsd': 87.79751114879761, 'freq': 3.985}, 'channel7': {'maxPsd': 95.96680704530192, 'freq': 3.985}, 'channel8': {'maxPsd': 102.58659595980386, 'freq': 3.985}}
    result_3985_on = {'channel1': {'maxPsd': 133.7091487774229, 'freq': 3.985}, 'channel2': {'maxPsd': 134.6277463439083, 'freq': 3.985}, 'channel3': {'maxPsd': 98.28967751844448, 'freq': 3.985}, 'channel4': {'maxPsd': 123.71481964489175, 'freq': 3.985}, 'channel5': {'maxPsd': 137.17538463141872, 'freq': 3.985}, 'channel6': {'maxPsd': 98.17070030274421, 'freq': 3.985}, 'channel7': {'maxPsd': 120.90226353746817, 'freq': 3.985}, 'channel8': {'maxPsd': 140.044190593653, 'freq': 3.985}}

    #5.005MHz
    dIn_5005 = 18.4
    result_5005_off = {'channel1': {'maxPsd': 84.44616464228199, 'freq': 5.005}, 'channel2': {'maxPsd': 89.95805298823043, 'freq': 5.005}, 'channel3': {'maxPsd': 96.01354774713104, 'freq': 5.005}, 'channel4': {'maxPsd': 93.1100468308548, 'freq': 5.005}, 'channel5': {'maxPsd': 105.1590913231004, 'freq': 5.005}, 'channel6': {'maxPsd': 93.13861063796539, 'freq': 5.005}, 'channel7': {'maxPsd': 95.69690437411332, 'freq': 5.005}, 'channel8': {'maxPsd': 99.69970330344928, 'freq': 5.005}}
    result_5005_on = {'channel1': {'maxPsd': 124.54919186327078, 'freq': 5.005}, 'channel2': {'maxPsd': 116.92606342384482, 'freq': 5.005}, 'channel3': {'maxPsd': 104.37262721653869, 'freq': 5.005}, 'channel4': {'maxPsd': 98.76740979986751, 'freq': 5.005}, 'channel5': {'maxPsd': 100.44385974585508, 'freq': 5.005}, 'channel6': {'maxPsd': 101.25978348698462, 'freq': 5.005}, 'channel7': {'maxPsd': 103.29113091418826, 'freq': 5.005}, 'channel8': {'maxPsd': 122.99105556641383, 'freq': 5.005}}

    #8.005MHz
    dIn_8005 = 23.4
    result_8005_off =  {'channel1': {'maxPsd': 88.21325353372366, 'freq': 8.005}, 'channel2': {'maxPsd': 95.6558527680718, 'freq': 8.005}, 'channel3': {'maxPsd': 103.41426154286428, 'freq': 8.005}, 'channel4': {'maxPsd': 101.82637272804546, 'freq': 8.005}, 'channel5': {'maxPsd': 100.66681070437113, 'freq': 8.005}, 'channel6': {'maxPsd': 92.64572731629336, 'freq': 8.005}, 'channel7': {'maxPsd': 96.90241658828813, 'freq': 8.005}, 'channel8': {'maxPsd': 91.05278631153185, 'freq': 8.005}}
    result_8005_on = {'channel1': {'maxPsd': 173.98869624175683, 'freq': 8.005}, 'channel2': {'maxPsd': 127.78244905000003, 'freq': 8.005}, 'channel3': {'maxPsd': 94.57846222496514, 'freq': 8.005}, 'channel4': {'maxPsd': 98.98117852480583, 'freq': 8.005}, 'channel5': {'maxPsd': 105.43153591916162, 'freq': 8.005}, 'channel6': {'maxPsd': 107.44039585959607, 'freq': 8.005}, 'channel7': {'maxPsd': 104.56008833447306, 'freq': 8.005}, 'channel8': {'maxPsd': 103.73170118758493, 'freq': 8.005}}

    #15.005MHz
    dIn_15005 = 47
    result_15005_off = {'channel1': {'maxPsd': 88.6973595443514, 'freq': 15.005}, 'channel2': {'maxPsd': 100.89255374231371, 'freq': 15.005}, 'channel3': {'maxPsd': 95.19673157222037, 'freq': 15.005}, 'channel4': {'maxPsd': 101.16032092929508, 'freq': 15.005}, 'channel5': {'maxPsd': 105.04254997359219, 'freq': 15.005}, 'channel6': {'maxPsd': 90.1925140427447, 'freq': 15.005}, 'channel7': {'maxPsd': 90.14646565985325, 'freq': 15.005}, 'channel8': {'maxPsd': 92.63836491294202, 'freq': 15.005}}
    result_15005_on ={'channel1': {'maxPsd': 384.8151188725788, 'freq': 15.005}, 'channel2': {'maxPsd': 159.95633486704264, 'freq': 15.005}, 'channel3': {'maxPsd': 174.1764097367352, 'freq': 15.005}, 'channel4': {'maxPsd': 176.81778963675413, 'freq': 15.005}, 'channel5': {'maxPsd': 141.973220739069, 'freq': 15.005}, 'channel6': {'maxPsd': 155.2381747322524, 'freq': 15.005}, 'channel7': {'maxPsd': 162.32123706179232, 'freq': 15.005}, 'channel8': {'maxPsd': 130.06933166201418, 'freq': 15.005}}

    #19.005MHz
    dIn_19005 = 54.8
    result_19005_off = {'channel1': {'maxPsd': 83.3299563511681, 'freq': 19.005}, 'channel2': {'maxPsd': 97.97375918782299, 'freq': 19.005}, 'channel3': {'maxPsd': 95.24280827870516, 'freq': 19.005}, 'channel4': {'maxPsd': 95.55505400694582, 'freq': 19.005}, 'channel5': {'maxPsd': 98.10646839576278, 'freq': 19.005}, 'channel6': {'maxPsd': 97.47999280268765, 'freq': 19.005}, 'channel7': {'maxPsd': 93.11779768704939, 'freq': 19.005}, 'channel8': {'maxPsd': 95.05501548332478, 'freq': 19.005}}
    result_19005_on = {'channel1': {'maxPsd': 89.1670280686576, 'freq': 19.005}, 'channel2': {'maxPsd': 97.12315654137652, 'freq': 19.005}, 'channel3': {'maxPsd': 94.16648352232174, 'freq': 19.005}, 'channel4': {'maxPsd': 92.69889936191515, 'freq': 19.005}, 'channel5': {'maxPsd': 93.3386720635052, 'freq': 19.005}, 'channel6': {'maxPsd': 101.31165488219273, 'freq': 19.005}, 'channel7': {'maxPsd': 95.1578059971207, 'freq': 19.005}, 'channel8': {'maxPsd': 86.96965579129618, 'freq': 19.005}}

    
    allResults = {}
    print("115kHz")
    allResults[115] = calcPsrr(result_115_on,result_115_off,dIn_115)
    print("305kHz")
    allResults[305] = calcPsrr(result_305_on,result_305_off,dIn_305)
    print("505kHz")
    allResults[505] = calcPsrr(result_505_on,result_505_off,dIn_505)
    print("755kHz")
    allResults[755] = calcPsrr(result_755_on,result_755_off,dIn_755)
    print("1005kHz")
    allResults[1005] = calcPsrr(result_1005_on,result_1005_off,dIn_1005)
    print("2005kHz")
    allResults[2005] = calcPsrr(result_2005_on,result_2005_off,dIn_2005)
    print("3005kHz")
    allResults[3005] = calcPsrr(result_3005_on,result_3005_off,dIn_3005)
    print("3985kHz")
    allResults[3985] = calcPsrr(result_3985_on,result_3985_off,dIn_3985)
    print("5005kHz")
    allResults[5005] = calcPsrr(result_5005_on,result_5005_off,dIn_5005)
    print("8005kHz")
    allResults[8005] = calcPsrr(result_8005_on,result_8005_off,dIn_8005)
    print("15005kHz")
    allResults[15005] = calcPsrr(result_15005_on,result_15005_off,dIn_15005)
    print("19005kHz")
    allResults[19005] = calcPsrr(result_19005_on,result_19005_off,dIn_19005)

    chResults = {}
    for freq in allResults:
      for chanName in allResults[freq] :
        if chanName not in chResults :
          chResults[chanName] = {"freq":[],"psrr":[]}
        chResults[chanName]["freq"].append(freq)
        chResults[chanName]["psrr"].append(allResults[freq][chanName] )

    #carrier method results
    freqCarrier = [115,305,505,755,2005,3005,8005,15005]
    psrrCarrier = [58.1,56.2,52.6,53.3,78,78.2,94.8,83.6]

    chLabels = {"channel1":"Channel 1","channel2":"Channel 2","channel3":"Channel 3","channel4":"Channel 4","channel5":"Channel 5","channel6":"Channel 6","channel7":"Channel 7","channel8":"Channel 8"}

    fig, axes = plt.subplots(1,1,figsize=(10, 6))
    for chanName in chResults :
      chLabel = chLabels[chanName]
      axes.plot(chResults[chanName]["freq"],chResults[chanName]["psrr"],".-",label=chLabel)
    axes.plot(freqCarrier,psrrCarrier,".-",label="Channel 2 Carrier")
    axes.set_xlabel('Frequency [kHz]', fontsize=20)
    axes.set_ylabel('PSRR [dB]',  fontsize=20)
    axes.tick_params(axis="x", labelsize=12)
    axes.tick_params(axis="y", labelsize=12)
    axes.set_title("PSRR vs Frequency", fontsize=20)
    #plt.xlim(1, 100)
    plt.ylim(00, 100)
    axes.set_xscale("log")
    axes.grid()
    axes.legend()

    fig.tight_layout()
    plt.show()
    
    

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
