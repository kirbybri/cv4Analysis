import numpy as np
from math import *
import matplotlib.pyplot as plt


def main():

    x_data =     [0.05,0.1  ,0.25,0.5,0.75,1.  ,5   ,8   ,10  ,19]
    y_data_ch1 = [57  ,59   ,52  ,67 ,73  ,74.9,85.4,83  ,80.8,82.6]
    y_data_ch4 = [58.2,67.3 ,71  ,76 ,87  ,85.4,81.5,80.9,79.8,81]
    #y_err = [0.5, 0.5, 0.5, 0.5]

    fig, axes = plt.subplots(1,1,figsize=(10, 6))
    axes.plot(x_data,y_data_ch1,".-",label="ASIC 2 CH1")
    axes.plot(x_data,y_data_ch4,".-",label="ASIC 2 CH4")
    #axes.errorbar(x_data,y_data,y_err,linestyle = 'None',marker=".")
    axes.set_xlabel('Frequency [MHz]', fontsize=20)
    axes.set_ylabel('PSRR [dB]',  fontsize=20)
    axes.tick_params(axis="x", labelsize=12)
    axes.tick_params(axis="y", labelsize=12)
    axes.set_title("PSRR vs Frequency", fontsize=20)
    #plt.xlim(1, 100)
    plt.ylim(50, 90)
    axes.set_xscale("log")
    axes.grid()
    axes.legend()

    fig.tight_layout()
    plt.show()
    
    

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
