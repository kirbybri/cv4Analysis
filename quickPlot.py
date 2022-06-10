import numpy as np
from math import *
import matplotlib.pyplot as plt


def main():

    x_data = [1.00, 5.00, 8., 18.]
    y_data = [72.6  , 69.9     , 69.7  , 64.5    ]
    y_err = [0.5, 0.5, 0.5, 0.5]

    fig, axes = plt.subplots(1,1,figsize=(10, 6))
    #axes.plot(x_data,y_data,".")
    axes.errorbar(x_data,y_data,y_err,linestyle = 'None',marker=".")
    axes.set_xlabel('Frequency [MHz]', fontsize=20)
    axes.set_ylabel('SNR [dB]',  fontsize=20)
    axes.tick_params(axis="x", labelsize=12)
    axes.tick_params(axis="y", labelsize=12)
    axes.set_title("SNR vs Frequency")
    plt.xlim(1, 100)
    plt.ylim(0, 120)
    axes.set_xscale("log")

    fig.tight_layout()
    plt.show()
    
    

#-------------------------------------------------------------------------
if __name__ == "__main__":
  main()
