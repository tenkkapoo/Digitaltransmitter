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
from receiver import *
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

        def ACPR_bit_sweep():
            s_source=NR_signal_generator()
            d=Transmitter(self)
            d.IOS.Members['A']=s_source.IOS.Members['out']
            s_source.run_gen()
            d.IOS.Members['B']=s_source.s_struct['Fs']
            d.IOS.Members['C']=s_source.BW
            d.model="py"
            d.init()
            d.f_c=700e6
            d.Fs=2400e6
            d.baseband_bits=2
            d.amplitude_bits=2
            d.phase_bits=2
            d.Inter_factor=128
            d.A_MAX=2
            i=0
            ACLRdict={"cartesian": 0, "polar": 0, "outphasing": 0, "multileveloutphasing1": 0, "multileveloutphasing2": 0, "multileveloutphasing3": 0, "multileveloutphasing4": 0}
            ACLRdict2={"cartesian": 0, "polar": 0, "outphasing": 0, "multileveloutphasing1": 0, "multileveloutphasing2": 0, "multileveloutphasing3": 0, "multileveloutphasing4": 0}
            keyslist=list(ACLRdict.keys())
            types=["cartesian", "polar", "outphasing", "multileveloutphasing"]
            for d.transmittertype in types:  #Loop for amplitude-bits
                ACLRarray_1ADJC=[]
                ACLRarray_2ADJC=[]
                """if d.transmittertype=="multileveloutphasing":
                    d.A_MAX+=1"""
                print(d.transmittertype)
                for d.bits in range(3, 13, 1):
                    d.run()
                    x, y, ACLR=plot_PSD(signal=d.IOS.Members['Z'].Data, BW=s_source.BW, BW_conf=s_source.BW_conf,
                    ACLR_BW=s_source.ACLR_BW, Fc=700e6, Fs=d.Fs, double_sided=False, zoom_plt=0,
                    no_plot=1)
                    ACLRarray_1ADJC.append(ACLR[1])
                    ACLRarray_2ADJC.append(ACLR[0])
                ACLRdict[keyslist[i]]=ACLRarray_1ADJC
                ACLRdict2[keyslist[i]]=ACLRarray_2ADJC
                i+=1
                d.init()
            """for d.transmittertype in types: #Loop for phase-bits
                ACLRarray_1ADJC=[]
                ACLRarray_2ADJC=[]
                print(d.transmittertype)
                for d.phase_bits in range(2, 14, 2):
                    d.run()
                    x, y, ACLR=plot_PSD(signal=d.IOS.Members['Z'].Data, BW=s_source.BW, BW_conf=s_source.BW_conf,
                    ACLR_BW=s_source.ACLR_BW, Fc=700e6, Fs=d.Fs, double_sided=False, zoom_plt=0,
                    no_plot=1)
                    ACLRarray_1ADJC.append(ACLR[1])
                    ACLRarray_2ADJC.append(ACLR[0])
                ACLRdict[d.transmittertype]=ACLRarray_1ADJC
                d.init()"""
            print(ACLRdict["outphasing"])
            bits=np.arange(3, 13, 1)
            f=plot.figure(1)
            plot.plot(bits, ACLRdict["cartesian"], linestyle="-.", label="Karteesinen")
            plot.plot(bits, ACLRdict["polar"], linestyle="--", label="Polaarinen")
            plot.plot(bits, ACLRdict["outphasing"], linestyle="-.", label="Poisvaiheistus")
            plot.plot(bits, ACLRdict["multileveloutphasing1"], linestyle=":", label="Monitasoinen poisvaiheistus A_MAX=2")
            """plot.plot(bits, ACLRdict["multileveloutphasing2"], linestyle=":", label="Monitasoinen poisvaiheistus A_MAX=2")
            plot.plot(bits, ACLRdict["multileveloutphasing3"], linestyle=":", label="Monitasoinen poisvaiheistus A_MAX=3")
            plot.plot(bits, ACLRdict["multileveloutphasing4"], linestyle=":", label="Monitasoinen poisvaiheistus A_MAX=4")"""
            plot.xlabel("Bitit")
            plot.ylabel("ACLR ensimmäinen kanava (dB)")
            plot.xticks(bits, bits)
            plot.legend()
            f.show()
            g=plot.figure(2)
            plot.plot(bits, ACLRdict2["cartesian"], linestyle="-.", label="Karteesinen")
            plot.plot(bits, ACLRdict2["polar"], linestyle="--", label="Polaarinen")
            plot.plot(bits, ACLRdict2["outphasing"], linestyle="-.", label="Poisvaiheistus")
            plot.plot(bits, ACLRdict2["multileveloutphasing1"], linestyle=":", label="Monitasoinen poisvaiheistus A_MAX=2")
            """plot.plot(bits, ACLRdict2["multileveloutphasing2"], linestyle=":", label="Monitasoinen poisvaiheistus A_MAX=2")
            plot.plot(bits, ACLRdict2["multileveloutphasing3"], linestyle=":", label="Monitasoinen poisvaiheistus A_MAX=3")
            plot.plot(bits, ACLRdict2["multileveloutphasing4"], linestyle=":", label="Monitasoinen poisvaiheistus A_MAX=4")"""
            plot.xlabel("Bitit")
            plot.ylabel("ACLR toinen kanava (dB)")
            plot.xticks(bits, bits)
            plot.legend()
            g.show()

        def EVM_bit_sweep():
            s_source=NR_signal_generator()
            d=Transmitter(self)
            r=receiver(self)
            d.IOS.Members['A']=s_source.IOS.Members['out']
            d.IOS.Members['Z']=r.IOS.Members['sig_in']
            s_source.run_gen()
            d.IOS.Members['B']=s_source.s_struct['Fs']
            d.IOS.Members['C']=s_source.BW
            s_source.IOS.Members['in_dem']=r.IOS.Members['out']
            d.model="py"
            d.init()
            d.f_c=700e6
            d.Fs=2400e6
            d.baseband_bits=2
            d.amplitude_bits=2
            d.phase_bits=2
            d.Inter_factor=128
            d.A_MAX=1
            """r.init()
            r.Fc=-700e6
            r.Rs_sim=d.Fs
            r.Rs_bb=s_source.s_struct['Fs']
            r.Factor=128"""
            EVMdict={"cartesian": 0, "polar": 0, "outphasing": 0, "multileveloutphasing": 0}
            types=["cartesian", "polar", "outphasing", "multileveloutphasing"]
            for d.transmittertype in types:  #Loop for amplitude-bits
                EVMarray=[]
                print(d.transmittertype)
                for d.bits in range(3, 13, 1):
                    d.run()
                    r.init()
                    r.Fc=-700e6
                    r.Rs_sim=d.Fs
                    r.Rs_bb=s_source.s_struct['Fs']
                    r.Factor=128
                    r.run()
                    s_source.run_dem()
                    s_source.run_EVM()
                    EVMarray.append(s_source.EVM[0][0]*100)
                EVMdict[d.transmittertype]=EVMarray
            bits=np.arange(3,13, 1)
            plot.plot(bits, EVMdict["cartesian"],  linestyle="-", label="Karteesinen")
            plot.plot(bits, EVMdict["polar"], linestyle="--", label="Polaarinen")
            plot.plot(bits, EVMdict["outphasing"], linestyle="-.", label="Poisvaiheistus")
            plot.plot(bits, EVMdict["multileveloutphasing"], linestyle=":", label="Monitasoinen poisvaiheistus A_MAX=1")
            plot.xlabel("Bitit")
            plot.ylabel("EVM (%)")
            plot.xticks(bits, bits)
            plot.legend()
            plot.show()
            """for d.A_MAX in np.arange(1, 5, 1):
                        EVMarray=[]
                        print(d.transmittertype)
                        for d.bits in range(3, 13, 1):
                            d.run()
                            r.init()
                            r.Fc=-700e6
                            r.Rs_sim=d.Fs
                            r.Rs_bb=s_source.s_struct['Fs']
                            r.Factor=128
                            r.run()
                            s_source.run_dem()
                            s_source.run_EVM()
                            EVMarray.append(s_source.EVM[0][0]*100)
                        EVMdict[d.transmittertype]=EVMarray"""

        ACPR_bit_sweep()
        #EVM_bit_sweep()
        """s_source=NR_signal_generator()
        d=Transmitter(self)
        r=receiver(self)
        d.IOS.Members['A']=s_source.IOS.Members['out']
        #d.IOS.Members['Z']=r.IOS.Members['sig_in']
        r.IOS.Members['sig_in']=d.IOS.Members['Z']
        s_source.run_gen()
        d.IOS.Members['B']=s_source.s_struct['Fs']
        d.IOS.Members['C']=s_source.BW
        s_source.IOS.Members['in_dem']=r.IOS.Members['out']
        d.model="py"
        d.init()
        d.f_c=700e6
        d.Fs=2400e6
        d.baseband_bits=4
        d.amplitude_bits=12
        d.phase_bits=12
        d.Inter_factor=128
        d.transmittertype="cartesian"
        d.run()
        r.init()
        r.Fc=-700e6
        r.Rs_sim=d.Fs
        r.Rs_bb=s_source.s_struct['Fs']
        r.Factor=128
        r.run()
        s_source.run_dem()
        s_source.run_EVM()
        text=f'EVM:{np.round(s_source.EVM[0][0]*100,2)}%'
        for i in range(0,s_source.BW.size):
            for j in range(0,len(s_source.BWP[i])):
                if  s_source.EVM[i].any():
                    if s_source.EVM[i][j]!=0:
                        plt.figure()
                        plt.plot(s_source.rxDataSymbols[i][j].real,s_source.rxDataSymbols[i][j].imag,'o')
                        plt.plot(s_source.cnstl[i][j].real,s_source.cnstl[i][j].imag,'o')
                        plt.text(0.025,0.975,text,usetex=plt.rcParams['text.usetex'],
                                horizontalalignment='left',verticalalignment='top',
                                multialignment='left',fontsize=plt.rcParams['legend.fontsize'],
                                fontweight='normal',transform=plt.gca().transAxes,
                                bbox=dict(boxstyle='square,pad=0',fc='#ffffffa0',ec='none'))
                        plt.title("Carrier "+str(i+1)+", frame "+str(j+1))
                        plt.show(block=False)
        print(s_source.EVM[0][0])
        _,ACLR=plot_PSD(signal=s_source.IOS.Members['out'].Data[:, 1]+1j*s_source.IOS.Members['out'].Data[:, 2], BW=s_source.BW, BW_conf=s_source.BW_conf,
        ACLR_BW=s_source.ACLR_BW, Fc=0, Fs=s_source.s_struct['Fs'], double_sided=True, zoom_plt=0,
        no_plot=0)
        _,ACLR=plot_PSD(signal=d.IOS.Members['Z'].Data, BW=s_source.BW, BW_conf=s_source.BW_conf,
        ACLR_BW=s_source.ACLR_BW, Fc=700e6, Fs=d.Fs, double_sided=False, zoom_plt=0,
        no_plot=0)
        _,ACLR2=plot_PSD(signal=r.IOS.Members['out'].Data[:, 1]+1j*r.IOS.Members['out'].Data[:, 2], BW=s_source.BW, BW_conf=s_source.BW_conf,
        ACLR_BW=s_source.ACLR_BW, Fc=0, Fs=s_source.s_struct['Fs'], double_sided=True, zoom_plt=0,
        no_plot=0)
        #_,EVM=s_source.measEVMdownlink(BW=s_source.BW, BWP=s_source.BWP, cnstl=s_source.s_struct['cnstl'], dem=r.IOS.Members['out'].Data[:, 1]+1j*r.IOS.Members['out'].Data[:, 2])
        _,ACLR=plot_PSD(signal=d.IOS.Members['Z'].Data, BW=s_source.BW, BW_conf=s_source.BW_conf,
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

