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
from plot_PSD import plot_PSD

class Transmitter(thesdk):
    @property
    def _classfile(self):
        return os.path.dirname(os.path.realpath(__file__)) + "/"+__name__

    def __init__(self,*arg): 
        self.print_log(type='I', msg='Inititalizing %s' %(__name__)) 
        self.proplist = [ 'Rs' ];    # Properties that can be propagated from parent
        self.Rs =  100e6;            # Sampling frequency
        self.IOS=Bundle()            # Pointer for input data
        self.IOS.Members['A']=IO()   # Pointer for input data
        self.IOS.Members['B']=IO()   # Baseband sampling frequency
        self.IOS.Members['C']=IO()   # Baseband bandwidth
        self.IOS.Members['Z']= IO()  # Output
        self.model='py';             # Can be set externally, but is not propagated
        self.par= False              # By default, no parallel processing
        self.queue= []               # By default, no parallel processing
        self.f_c=1000e6

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

        def interpolator(time, I, Q, f_c, BW):
            newtime=np.arange(0, time[-1], 1/(2*(f_c+BW)))
            new_I1=np.interp(newtime, time, I)
            new_Q1=np.interp(newtime, time, Q) 
            b, a=signal.iirfilter(15, BW/2, btype='low', fs=2*(f_c+BW/2))
            new_I2=signal.lfilter(b, a, new_I1)
            new_Q2=signal.lfilter(b, a, new_Q1)
            return newtime, new_I2, new_Q2

        def quantization(I, Q, bits):
            I_step_size=(max(I)-min(I))/(2**bits)
            Q_step_size=(max(Q)-min(Q))/(2**bits)
            I_quantized=I_step_size*np.floor(2**(bits-1)*I/max(I))+2**-1
            return I_quantized


        def cartesian(time, I, Q, f_c):
            w_c=2*math.pi*f_c
            output=I*np.cos(w_c*time)-Q*np.sin(w_c*time)
            return output

        def polar(time, I, Q, f_c):
            w_c=2*math.pi*f_c
            r=np.sqrt(I**2+Q**2)
            phi=np.arctan(Q/I)
            output=r*np.cos(w_c*time+phi)
            return output

        def outphasing(time, I, Q, f_c):
            w_c=2*math.pi*f_c
            r=np.sqrt(I**2+Q**2)
            phi1=np.arctan(I/Q)+np.arctan(r/np.max(r))
            phi2=np.arctan(I/Q)-np.arctan(r/np.max(r))
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
            phi1=np.arctan(I/Q)+theta_MO
            phi2=np.arctan(I/Q)-theta_MO
            S1=1/2*np.cos(w_c*time+phi1)
            S2=1/2*np.cos(w_c*time+phi2)
            output=(A_MO/(2*A_MAX))*(S1+S2)
            return output

        #self.IOS.Members['Z'].Data=inval[:, 1]+1j*inval[:, 2]
        f=plot.figure(1)
        plot.plot(inval[:, 0], inval[:, 1])
        f.show()
        I = quantization(inval[:, 1], inval[:, 2], 4)
        G=plot.figure(2)
        plot.plot(inval[:, 0], I)
        G.show()
        time, I, Q=interpolator(inval[:, 0], inval[:, 1], inval[:, 2], self.f_c, BW)
        print(max(I))
        print(min(I))
        print(max(Q))
        print(min(Q))
        #out=multileveloutphasing(time, I, Q, self.f_c)
        if self.par:
            self.queue.put(out)
        #self.IOS.Members['Z'].Data=out

    
 
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

