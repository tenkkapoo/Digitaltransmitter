"""
=========
Transmitter
=========

Transmitter model template The System Development Kit
Used as a template for all TheSyDeKick Entities.

Current docstring documentation style is Numpy
https://numpydoc.readthedocs.io/en/latest/format.html

This text here is to remind you that documentation is important.
However, youu may find it out the even the documentation of this 
entity may be outdated and incomplete. Regardless of that, every day 
and in every way we are getting better and better :).

Initially written by Marko Kosunen, marko.kosunen@aalto.fi, 2017.

"""

import os
import sys
if not (os.path.abspath('../../thesdk') in sys.path):
    sys.path.append(os.path.abspath('../../thesdk'))

from thesdk import *
import pdb
import numpy as np
import math
import matplotlib.pyplot as plot
from scipy import signal
from scipy import interpolate
from plot_PSD import plot_PSD

class Transmitter(thesdk):
    @property
    def _classfile(self):
        return os.path.dirname(os.path.realpath(__file__)) + "/"+__name__

    def __init__(self,*arg): 
        self.print_log(type='I', msg='Inititalizing %s' %(__name__)) 
        self.proplist = ['f_c', 'Fs', 'Inter_factor', 'amplitude_bits', 'phase_bits'];    # Properties that can be propagated from parent
        self.IOS=Bundle()            # Pointer for input data
        self.IOS.Members['A']=IO()   # Pointer for input data
        self.IOS.Members['B']=IO()   # Baseband sampling frequency
        self.IOS.Members['C']=IO()   # Baseband bandwidth
        self.IOS.Members['Z']= IO()  # Output
        self.model='py';             # Can be set externally, but is not propagated
        self.par= False              # By default, no parallel processing
        self.queue= []               # By default, no parallel processing
        self.f_c=700e6               # Carrier frequency
        self.Fs=2400e6               # Sampling rate
        self.Inter_factor=64         # Interpolation factor. Minimum 

        if len(arg)>=1:
            parent=arg[0]
            self.copy_propval(parent,self.proplist)
            self.parent =parent;

        self.init()

    def init(self):
        pass #Currently nohing to add


    def main(self):
        '''Guideline. Isolate python processing to main method.
        
        To isolate the interna processing from IO connection assigments, 
        The procedure to follow is
        1) Assign input data from input to local variable
        2) Do the processing
        3) Assign local variable to output

        '''
        inval=self.IOS.Members['A'].Data
        Fs_baseband=self.IOS.Members['B']
        BW=self.IOS.Members['C']
        #pdb.set_trace()

        def interpolator(time, I, Q, Fs, Fs_baseband, BW):
            Fs_inter=Fs_baseband
            #Loop for upsampling by two and half-band filtering np.log2(self.Inter_factor) times.
            for i in range(0, int(np.log2(self.Inter_factor))): 
                Fs_inter=2*Fs_inter
                print(Fs_inter)
                I=np.transpose(np.matrix(np.array(I).reshape(1,-1)))
                I=np.pad(I,((0,0),(0,int(1))),'constant',constant_values=0)
                I=np.concatenate(I)
                I=np.trim_zeros(I)
                Q=np.transpose(np.matrix(np.array(Q).reshape(1,-1)))
                Q=np.pad(Q,((0,0),(0,int(1))),'constant',constant_values=0)
                Q=np.concatenate(Q)
                Q=np.trim_zeros(Q)
                BW=int(BW) 
                b=signal.remez(int(61), [0, 0.1*Fs_inter, 0.28*Fs_inter, 0.5*Fs_inter], [1, 0], Hz=int(Fs_inter))
                I=np.fft.ifft(np.fft.fft(I.reshape((-1,1))[:,0],len(I))*np.fft.fft(b,len(I))).reshape((-1,1))
                Q=np.fft.ifft(np.fft.fft(Q.reshape((-1,1))[:,0],len(Q))*np.fft.fft(b,len(Q))).reshape((-1,1))
                w, h=signal.freqz(b, [1])
                fig = plot.figure()
                ax = fig.add_subplot(111)
                ax.plot(0.5*Fs_inter*w/np.pi, 20*np.log10(np.abs(h)))
                ax.set_xlim(0, 0.5*Fs)
                plot.show(block=False)
            """Fs_inter=self.Inter_factor*Fs_baseband
            b=signal.remez(int(181), [0, BW/2, 2.1*BW/2, 0.5*Fs_inter], [1, 0], Hz=int(Fs_inter))
            I_fil=np.fft.ifft(np.fft.fft(I.reshape((-1,1))[:,0],len(I))*np.fft.fft(b,len(I))).reshape((-1,1))
            Q_fil=np.fft.ifft(np.fft.fft(Q.reshape((-1,1))[:,0],len(Q))*np.fft.fft(b,len(Q))).reshape((-1,1))
            w, h=signal.freqz(b, [1])
            fig = plot.figure()
            ax = fig.add_subplot(111)
            ax.plot(0.5*Fs_inter*w/np.pi, 20*np.log10(np.abs(h)))
            ax.set_xlim(0, 0.5*Fs)
            plot.show(block=False)"""
            newtime=np.arange(0, len(I))/Fs_inter
            return newtime, np.real(I.reshape(1, -1)[0]), np.real(Q.reshape(1, -1)[0])


        def quantization(I, Q, bits):
            I_step_size=(max(I)-min(I))/(2**bits-1)
            Q_step_size=(max(Q)-min(Q))/(2**bits-1)
            I_quantized=I_step_size*(np.floor(I/I_step_size)+1/2)
            Q_quantized=Q_step_size*(np.floor(Q/Q_step_size)+1/2)
            """binsI=np.linspace(min(I), max(I), num=2**bits, endpoint=False)
            I_quantized=np.digitize(I, binsI)
            binsQ=np.linspace(min(Q), max(Q), num=2**bits, endpoint=False)
            Q_quantized=np.digitize(Q, binsQ)"""
            return I_quantized, Q_quantized


        def cartesian(time, I, Q, f_c, bits):
            w_c=2*math.pi*f_c
            output=I*np.cos(w_c*time)-Q*np.sin(w_c*time)
            output_stepsize=(max(output)-min(output))/(2**bits-1)
            output_quantized=output_stepsize*(np.floor(output/output_stepsize)+1/2)
            return output_quantized

        def polar(time, I, Q, f_c, amplitude_bits, phase_bits):
            w_c=2*math.pi*f_c
            r=np.sqrt(I**2+Q**2)
            """binsr=np.linspace(min(r), max(r), num=2**amplitude_bits, endpoint=False)
            r_quantized=np.digitize(r, binsr)"""
            r_step_size=(max(r)-min(r))/(2**amplitude_bits-1)
            r_quantized=r_step_size*(np.floor(r/r_step_size)+1/2)
            print(len(np.unique(r_quantized)))
            phi=np.arctan2(Q, I)
            phi_step_size=(max(phi)-min(phi))/(2**phase_bits-1)
            phi_quantized=phi_step_size*(np.floor(phi/phi_step_size)+1/2)
            print(len(np.unique(phi_quantized)))
            output=r_quantized*np.cos(w_c*time+phi_quantized)
            return output

        def outphasing(time, I, Q, f_c):
            w_c=2*math.pi*f_c
            r=np.sqrt(I**2+Q**2)
            phi1=np.arctan2(Q, I)+np.arccos(r/max(r))
            phi2=np.arctan2(Q, I)-np.arccos(r/max(r))
            S1=1/2*np.cos(w_c*time+phi1)
            S2=1/2*np.cos(w_c*time+phi2)
            output=S1+S2
            return output

        def multileveloutphasing(time, I, Q, f_c):
            w_c=2*math.pi*f_c
            r=np.sqrt(I**2+Q**2)
            A_MAX=np.max(2)
            A_MO=np.ceil(r*A_MAX)
            theta_MO=np.arccos(r*A_MAX/A_MO)
            phi1=np.arctan2(Q, I)+theta_MO
            phi2=np.arctan2(Q, I)-theta_MO
            S1=1/2*np.cos(w_c*time+phi1)
            S2=1/2*np.cos(w_c*time+phi2)
            output=(A_MO/(2*A_MAX))*(S1+S2)
            return output

        #self.IOS.Members['Z'].Data=inval[:, 1]+1j*inval[:, 2]
        self.Fs=Fs_baseband*self.Inter_factor
        I, Q=quantization(inval[:, 1], inval[:, 2], self.amplitude_bits)  #Input quantization
        time, I, Q=interpolator(inval[:, 0], I, Q, self.Fs, Fs_baseband, BW)
        #I_quant, Q_quant = quantization(I, Q, 6)
        #self.IOS.Members['Z'].Data=I+1j*Q
        out=cartesian(time, I, Q, self.f_c, self.amplitude_bits)
        #print(len(np.unique(out)))
        #out=polar(time, I, Q, self.f_c, self.amplitude_bits, self.phase_bits)
        #out=outphasing(time, I, Q, self.f_c)
        #out=multileveloutphasing(time, I, Q, self.f_c)
        if self.par:
            self.queue.put(out)
        self.IOS.Members['Z'].Data=out

    
 
    def run(self,*arg):
        '''Guideline: Define model depencies of executions in `run` method.

        '''
        if len(arg)>0:
            self.par=True      #flag for parallel processing
            self.queue=arg[0]  #multiprocessing.queue as the first argument
        if self.model=='py':
            self.main()

if __name__=="__main__":
    import matplotlib.pyplot as plt
    from  Transmitter import *
    from  Transmitter.controller import controller as Transmitter_controller
    import pdb
    import math
    time=np.arange(0, 10, 0.1)
    indata=np.array([time, np.cos(time), np.sin(time)])
    models=[ 'py']
    duts=[]
    for model in models:
        d=Transmitter()
        duts.append(d) 
        d.model=model
        d.IOS.Members['A'].Data=indata
        d.init()
        d.run()
    input()

