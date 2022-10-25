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
        self.IOS.Members['Z']= IO()
        self.model='py';             # Can be set externally, but is not propagated
        self.par= False              # By default, no parallel processing
        self.queue= []               # By default, no parallel processing

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
        time=inval[:, 0]
        I=inval[:, 1]
        Q=inval[:, 2]
        #self.IOS.Members['Z'].Data=I+1j*Q
        #pdb.set_trace()

        def cartesian(time, I, Q):
            w_c=2*math.pi*30e9
            output=I*np.cos(w_c*time)-Q*np.sin(w_c*time)
            return output

        def polar(time, I, Q):
            w_c=1000e6/2*math.pi
            r=np.sqrt(I**2+Q**2)
            phi=np.arctan(I/Q)
            output=r*np.cos(w_c*time+np.pi/2)
            return output

        def outphasing(time, I, Q):
            w_c=2*math.pi1000e6
            r=np.sqrt(I**2+Q**2)
            phi1=np.arctan(I/Q)+np.arctan(r/np.max(r))
            phi2=np.arctan(I/Q)-np.arctan(r/np.max(r))
            S1=1/2*np.cos(w_c*time+phi1)
            S2=1/2*np.cos(w_c*time+phi2)
            output=S1+S2
            plot.plot(time, output)
            plot.show()
            return output

        def multileveloutphasing(time, I, Q):
            w_c=1
            r=np.sqrt(I**2+Q**2)
            A_MAX=np.max(2)
            A_MO=math.ceil(r*A_MAX)
            theta_MO=np.arccos(r*A_MAX/A_MO)
            phi1=np.arctan(I/Q)+theta_MO
            phi2=np.arctan(I/Q)-theta_MO
            S1=1/2*np.cos(w_c*time+phi1)
            S2=1/2*np.cos(w_c*time+phi2)
            output=(A_MO/(2*A_MAX))*(S1+S2)
            plot.plot(time, output)
            plot.show()
            return output

        out=cartesian(time, I, Q)
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

