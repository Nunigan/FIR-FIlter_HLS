#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 14:57:38 2022

@author: nunigan
"""
from tkinter import X
import numpy as np
import matplotlib.pyplot as plt
import scipy
from datetime import datetime



class fir():
    def __init__(self, path):

        #------------------------------------------------
        # Create a signal for demonstration.
        #------------------------------------------------
        self.path = path
        self.sample_rate = 100.0
        self.nsamples = 1024
        self.t = np.arange(self.nsamples) / self.sample_rate
        self.x = (np.cos(2*np.pi*0.5*self.t) + 0.2*np.sin(2*np.pi*2.5*self.t+0.1) + \
                0.2*np.sin(2*np.pi*15.3*self.t) + 0.1*np.sin(2*np.pi*16.7*self.t + 0.1) + \
                    0.1*np.sin(2*np.pi*23.45*self.t+.8))/2

        #------------------------------------------------
        # Create a FIR filter and apply it to x.
        #------------------------------------------------
        
        # The Nyquist rate of the signal.
        self.nyq_rate = self.sample_rate / 2.0
        
        # The desired width of the transition from pass to stop,
        # relative to the Nyquist rate.  We'll design the filter
        # with a 5 Hz transition width.
        self.width = 5.0/self.nyq_rate
        
        # The desired attenuation in the stop band, in dB.
        self.ripple_db = 60.0
        
        # Compute the order and Kaiser parameter for the FIR filter.
        self.N, self.beta = scipy.signal.kaiserord(self.ripple_db, self.width)
        
        # The cutoff frequency of the filter.
        self.cutoff_hz = 10.0
        
        # Use firwin with a Kaiser window to create a lowpass FIR filter.
        self.taps = scipy.signal.firwin(self.N, self.cutoff_hz/self.nyq_rate, window=('kaiser', self.beta))
        
        # Use lfilter to filter x with the FIR filter.
        self.filtered_x =scipy.signal.lfilter(self.taps, 1.0, self.x)

        self.shift_reg = np.zeros(self.N)

        np.savetxt("taps.csv", self.taps, delimiter=",")

    def comp(self):
        acc_arr = []
        out = np.zeros(self.nsamples)
        for j in range(self.nsamples):
            acc = 0
            for i in range(self.N-1, 0, -1):
                self.shift_reg[i] = self.shift_reg[i - 1]
                acc += self.shift_reg[i] * self.taps[i]
                acc_arr.append(acc)
            acc += self.x[j] * self.taps[0]
            acc_arr.append(acc)

            self.shift_reg[0] = self.x[j]
            out[j] = acc
        return out, acc_arr

    def fir(x, h):
        nsamples = len(x)
        N = len(h)
        y = np.zeros(len(x))
        shift_reg = np.zeros(len(N))

        for j in range(nsamples):
            acc = 0
            for i in range(N-1, 0, -1):
                shift_reg[i] = shift_reg[i - 1]
                acc += shift_reg[i] * h[i]
            acc += x[j] * h[0]

            shift_reg[0] = x[j]
            y[j] = acc
        return y
        
    def plot(self, filtered_x):
        #------------------------------------------------
        # Plot the FIR filter coefficients.
        #------------------------------------------------
        
        if filtered_x == []:
            filtered_x = self.filtered_x
        
        
        plt.figure(1, figsize=(8,4))
        plt.plot(self.taps, 'bo-', linewidth=2)
        plt.title('Filter Coefficients (%d taps)' % self.N)
        plt.grid(True)
        
        #------------------------------------------------
        # Plot the magnitude response of the filter.
        #------------------------------------------------
        
        plt.figure(2, figsize=(8,4))
        plt.clf()
        w, h = scipy.signal.freqz(self.taps, worN=8000)
        plt.plot((w/np.pi)*self.nyq_rate, np.absolute(h), linewidth=2)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Gain')
        plt.title('Frequency Response')
        plt.ylim(-0.05, 1.05)
        plt.grid(True)
        
        # Upper inset plot.
        ax1 = plt.axes([0.42, 0.6, .45, .25])
        plt.plot((w/np.pi)*self.nyq_rate, np.absolute(h), linewidth=2)
        plt.xlim(0,8.0)
        plt.ylim(0.9985, 1.001)
        plt.text(4,1.0005,"Passband Ripple",bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
        plt.grid(True)
        
        # Lower inset plot
        ax2 = plt.axes([0.42, 0.25, .45, .25])
        plt.plot((w/np.pi)*self.nyq_rate, np.absolute(h), linewidth=2)
        plt.xlim(12.0, 20.0)
        plt.text(16, 0.002,"Stopband Ripple",bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
        plt.ylim(0.0, 0.0025)
        plt.grid(True)
        

        #------------------------------------------------
        # Plot the original and filtered signals.
        #------------------------------------------------
        
        # The phase delay of the filtered signal.
        delay = 0.5 * (self.N-1) / self.sample_rate
        
        plt.figure(3, figsize=(8,4))
        # Plot the original signal.
        plt.plot(self.t, self.x)
        # Plot the filtered signal, shifted to compensate for the phase delay.
        plt.plot(self.t-delay, filtered_x, 'r-')
        # Plot just the "good" part of the filtered signal.  The first N-1
        # samples are "corrupted" by the initial conditions.
        plt.plot(self.t[self.N-1:]-delay, filtered_x[self.N-1:], 'g', linewidth=4)
        
        plt.xlabel('t')
        plt.grid(True)
        
        plt.show()


    def wirte_c_arrays(self,I,W):

        with open(self.path+'/fir.h', 'w') as file:
            file.writelines('/* fir HLS, autogenerated File, Michael Schmid, {} */ \n \n'.format(datetime.now().strftime("%d/%m/%Y, %H:%M:%S")))
             
            file.writelines('#include <stdint.h> \n')
            file.writelines('#include "hls_math.h" \n')
            file.writelines('#include "ap_fixed.h"\n')

            file.writelines('\n')

            file.writelines('#define NUM_TAPS {}\n'.format(self.N))
            file.writelines('#define SIZE {}\n'.format(self.nsamples))
            file.writelines('#define W {}\n'.format(W))
            file.writelines('#define I {}\n'.format(I))
            file.writelines('\n')

            file.writelines('static const float taps[{}] = \n'.format(self.N))
            file.writelines('  {')
             
            for j in range(self.N):
                if j == self.N -1:
                    file.writelines('{}'.format(self.taps[j]))
                else:
                    if j % 4 == 0: # linebreak
                        file.writelines('\n        {},'.format(self.taps[j]))
                    else:
                        file.writelines('{},'.format(self.taps[j]))
            file.writelines('  };\n')
            file.writelines('\n')


            file.writelines('static const float x[{}] = \n'.format(self.nsamples))
            file.writelines('  {')
             
            for j in range(self.nsamples):
                if j == self.nsamples -1:
                    file.writelines('{}'.format(self.x[j]))
                else:
                    if j % 4 == 0: # linebreak
                        file.writelines('\n        {},'.format(self.x[j]))
                    else:
                        file.writelines('{},'.format(self.x[j]))
            file.writelines('  };\n')
            file.writelines('\n')

            file.writelines('static const float filtered_x[{}] = \n'.format(self.nsamples))
            file.writelines('  {')
             
            for j in range(self.nsamples):
                if j == self.nsamples -1:
                    file.writelines('{}'.format(self.filtered_x[j]))
                else:
                    if j % 4 == 0: # linebreak
                        file.writelines('\n        {},'.format(self.filtered_x[j]))
                    else:
                        file.writelines('{},'.format(self.filtered_x[j]))
            file.writelines('  };\n')
            file.writelines('\n')



            file.writelines('static const ap_fixed<W,I> taps_fixed[{}] = \n'.format(self.N))
            file.writelines('  {')
             
            for j in range(self.N):
                if j == self.N -1:
                    file.writelines('{}'.format(self.taps[j]))
                else:
                    if j % 4 == 0: # linebreak
                        file.writelines('\n        {},'.format(self.taps[j]))
                    else:
                        file.writelines('{},'.format(self.taps[j]))
            file.writelines('  };\n')
            file.writelines('\n')


# test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if __name__ == '__main__':

    fir_obj = fir(".")
    # fir_obj.plot([])
    fir_obj.wirte_c_arrays(I=1,W=12)
    out, arr = fir_obj.comp()
    fir_obj.plot(out)
