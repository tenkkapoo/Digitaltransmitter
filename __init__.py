"""
==================
Inverter testbench
==================

This class demonstrates how you can construct testbeches for your model
Entities to gain more versatilyty that can be provided by self tests in the 
mainguard.

Initially written by Marko Kosunen, marko.kosunen@aalto.fi, 2022.

"""

import os
import sys
if not (os.path.abspath('../../thesdk') in sys.path):
    sys.path.append(os.path.abspath('../../thesdk'))

from thesdk import *

from Transmitter import *
from NR_signal_generator import NR_signal_generator
from plot_PSD import plot_PSD
from math import pi


import numpy as np

class transmitter_testbench(thesdk):
    @property
    def _classfile(self):
        return os.path.dirname(os.path.realpath(__file__)) + "/"+__name__

    def __init__(self,*arg): 
        self.print_log(type='I', msg='Inititalizing %s' %(__name__)) 
    
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
        s_source=NR_signal_generator()
        d=Transmitter(self)
        d.IOS.Members['A']=s_source.IOS.Members['out']
        s_source.run_gen()
        d.IOS.Members['B']=s_source.s_struct['Fs']
        d.IOS.Members['C']=s_source.BW
        d.model="py"
        d.init()
        d.run()
        _,ACLR=plot_PSD(signal=d.IOS.Members['Z'].Data, BW=s_source.BW, BW_conf=s_source.BW_conf,
        ACLR_BW=s_source.ACLR_BW, Fc=0, Fs=s_source.s_struct['Fs'], double_sided=False, zoom_plt=0,
        no_plot=0)
        """_,ACLR=plot_PSD(signal=d.IOS.Members['Z'].Data, BW=s_source.BW, BW_conf=s_source.BW_conf,
        ACLR_BW=s_source.ACLR_BW, Fc=1000e6, Fs=1000e6+100e6, double_sided=True, zoom_plt=0,
        no_plot=0)"""

    def run(self,*arg):
        '''Guideline: Define model depencies of executions in `run` method.

        '''
        if self.model=='py':
            self.main()


if __name__=="__main__":
    from  inverter_testbench import *
    import pdb
    tb=transmitter_testbench()
    tb.models=['py']
    #tb.configuration='parallel'
    tb.run()
    #This is here to keep the images visible
    #For batch execution, you should comment the following line 
    input()

