import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.signal import lfilter

class SIM_PULSE(object):
    #__INIT__#
    def __init__(self):
        self.fileName = "default"
        self.pulseTimes = None
        self.triangleShapeVals = None
        self.pulseShapeVals = None
        self.sampTimes = None
        self.noiseVals = None
    
    def makePulseShape(self):
        ext_freq = 1200e6  # arbitrary waveform generator frequency, Hz
        sim = float(1 / ext_freq)  # simulation time resolution
        rise = 2.4e-9  # rise time
        width = 400e-9  # width of pulse
        period = float(1600e-9 - 1 * sim)  # total length of signal
        amp = 1  # amplitude of signal
        delay = 10*sim # 0e-9
        tau = 21.5e-9  # RC time constant
        time1 = np.linspace(0, period, math.floor(period / sim))  # time array
        pulse = np.zeros(len(time1))

        # Form signal wave
        i = 0
        for t in time1:
            if t < rise + delay and t > delay:
                pulse[i] = amp * (t - delay) / rise  # rising ramp
            elif t < width + delay and t > delay:
                pulse[i] = amp * (width - t + delay) / (width - rise)  # falling ramp
            else:
                pulse[i] = 0  # zeros elsewhere
            i += 1
            
        # bilinear transform from laplace to z
        d = 2 * tau / sim
        c = d / ((1 + d) ** 3)
        f = (1 - d) / (1 + d)
        b_old = [1, 1, -1, -1]
        b = [c * j for j in b_old]  # transfer function numerator
        a = [1, 3 * f, 3 * f ** 2, f ** 3]  # transfer function denominator
        filtered_pulse = lfilter(b, a, pulse)

        self.pulseTimes = time1
        self.triangleShapeVals = pulse
        self.pulseShapeVals = filtered_pulse
        return None
    
    def simplePlot(self,vals_x,vals_y):
        fig, axes = plt.subplots(1,1,figsize=(10, 6))
        axes.plot(vals_x,vals_y,"-.",label="")
        axes.grid()
        fig.tight_layout()
        plt.show()
        return None

    def printVals(self):
        for sampNum,sampTime in enumerate(self.pulseTimes):
            time = round(sampTime,5)
            val1 = round(self.triangleShapeVals[sampNum],5)
            val2 = round(self.pulseShapeVals[sampNum],5)
            print(sampNum,"\t",time,"\t",val1,"\t",val2)
        for sampNum,sampTime in enumerate(self.sampTimes):
            time = round(sampTime,5)
            val1 = round(self.noiseVals[sampNum],5)
            print(sampNum,"\t",time,"\t",val1)
        return None
    
    def makeNoise(self):
        sampPeriod = 25E-9
        numSamp = 8000
        #time1 = np.linspace(0, period, sampPeriod)
        pedMean = 6500.
        noiseRms = 1.5
        self.sampTimes = np.zeros(numSamp)
        self.noiseVals = np.zeros(numSamp)
        for sampNum in range(0,numSamp,1) :
            sampTime = sampPeriod*sampNum
            sampNoise = noiseRms*np.random.randn()
            self.sampTimes[sampNum] = sampTime
            self.noiseVals[sampNum] = sampNoise + pedMean
        return None
    
    def test(self):
        self.makePulseShape()
        self.makeNoise()
        self.printVals()
        #self.simplePlot(self.pulseTimes,self.pulseShapeVals)
        #self.simplePlot(self.sampTimes,self.noiseVals)

        return None
            
def main():
    print("HELLO")
    simPulse = SIM_PULSE()
    simPulse.test()
            
#-------------------------------------------------------------------------
if __name__ == "__main__":
    main()