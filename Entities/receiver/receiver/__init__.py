"""
========
Receiver
========

Inverter model template The System Development Kit
Used as a template for all TheSyDeKick Entities.

Current docstring documentation style is Numpy
https://numpydoc.readthedocs.io/en/latest/format.html

This text here is to remind you that documentation is iportant.
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
from rtl import *
from spice import *
import plot_PSD as plot_func
import pdb
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
class receiver(rtl,spice,thesdk):
    @property
    def _classfile(self):
        return os.path.dirname(os.path.realpath(__file__)) + "/"+__name__

    def __init__(self,*arg):
        """ The attributes and paremeters of the class

        Attributes
        ----------
        Fc : integer
            Carrier frequency
        factor : integer
            Factor for resampling. Preferably 2^x for stabel functionality
        am_resolution_tran : integer
            Number of bits for amplitude resolution
        phase_resolution_tran : integer
            Number of bits for phase resolution
        modulation : string
            Modulation type of transmitter. (cartesian , polar , outphasing , multi_outphasing)
        mode : string
            Operation mode of chip, transmitter or reciever. (tran, rec)
        """

        self.print_log(type='I', msg='Inititalizing %s' %(__name__)) 
        self.proplist = [ 'Fc' ,'Rs_sim','Rs_bb','new_model' ];    # Properties that can be propagated from parent
        self.IOS=Bundle()
        #self.IOS.Members['t_in']=IO() # Pointer for input data
        self.IOS.Members['sig_in']=IO() # Pointer for input data
        #self.IOS.Members['in']=IO() # Pointer for input data
        self.IOS.Members['out']=IO() # Pointer for input data
        #self.IOS.Members['additional_interp_factor']=IO()
        self.model='py';             # Can be set externally, but is not propagated
        self.par= False              # By default, no parallel processing
        self.queue= []               # By default, no parallel processing
        self.IOS.Members['control_write']= IO()

        #self.factor=1
        #self.pm_interp_facto=10
        self.Fc=0
        self.Rs_sim=4e9
        self.Rs_bb=200e6
        self.new_model=True

        # File for control is created in controller
        
        if len(arg)>=1:
            parent=arg[0]
            self.copy_propval(parent,self.proplist)
            self.parent =parent;

            self.Fc=-1*self.Fc
       
        self.init()

    def init(self):
        pass #Currently nohing to add

    def init_params(self,**kwargs):
        """ Method for setting variable values from outside of object
            Parameters
            ----------
            Same as in __init__

            Example
            -------
            x.init_param(modulation="cartesian",Fc=1e9)

        """

        if 'modulation' in kwargs:
            self.modulation=kwargs.get('modulation')
        if 'factor' in kwargs:
            self.factor=kwargs.get('factor')
        if 'Fc' in kwargs:
            self.Fc=kwargs.get('Fc')
        if 'am_resolution' in kwargs:
            self.am_resolution=kwargs.get('am_resolution')
        if 'phase_resolution' in kwargs:
            self.phase_resolution=kwargs.get('phase_resolution')
        if 'mode' in kwargs:
            self.mode=kwargs.get('mode')






    def main(self):
        '''Guideline. Isolate python processing to main method.
        
        To isolate the interna processing from IO connection assigments, 
        The procedure to follow is
        1) Assign input data from input to local variable
        2) Do the processing
        3) Assign local variable to output

        '''
        #pdb.set_trace()

       # signal=self.IOS.Members['in'].Data
       # if signal is None:
       #     self.t=self.IOS.Members['t_in'].Data
       #     im_zeros=np.zeros(len(self.t))
       #     s=self.IOS.Members['sig_in'].Data
       #     self.s=np.transpose(np.vstack((self.t, s, im_zeros)))
       # #self.Fs=self.IOS.Members['B'].Data
       # #self.Fc=self.IOS.Members['C'].Data
       # else:
       #     self.t=signal[:,0]
       #     self.s=signal

        #pdb.set_trace()
        self.s=self.IOS.Members['sig_in'].Data
        t=np.arange(0,len(self.s))/self.Rs_sim
        im_zeros=np.zeros(len(t))
        self.s_input=np.transpose(np.vstack((t, self.s, im_zeros)))
        # Debug values that should be deleted when whole chain is changed.
        #self.Fs_RF=1/(self.t[1])
        #add_factor=self.IOS.Members['additional_interp_factor'].Data
        #if add_factor is None:
        #    self.Fs_BB=1/(self.t[1])*self.factor
        #else:
        #    self.Fs_BB=1/(self.t[1])*self.factor*(1/add_factor)






        self.interp, self.decim=self.def_factor()

        self.print_log(type='I', msg='Downconverting...')
        self.s=self.demodulate()
        self.print_log(type='I', msg='..Done!')
        self.print_log(type='I', msg='Decimating...')
        #self.s=self.interpolate()
        self.s=self.resample()
        self.print_log(type='I', msg='..Done!')

        self.IOS.Members['out'].Data=self.s


       

    def run(self,*arg):
        '''Guideline: Define model depencies of executions in `run` method.

        '''
        if len(arg)>0:
            self.par=True      #flag for parallel processing
            self.queue=arg[0]  #multiprocessing.queue as the first argument
        if self.model=='py':
            self.main()
   

    def def_factor(self):
        """ Method for calculating interpolation and decimation factors for resampling. 
            Factors are smallest possible integers based on overall resampling factor.

            Example
            -------
            self.def_factor()

        """


        #Fs_org=1/self.t[1]
        #a=max(Fs_org*self.factor,Fs_org)
        #b=min(Fs_org*self.factor,Fs_org)
        a=max(self.Rs_sim,self.Rs_bb)
        b=min(self.Rs_sim,self.Rs_bb)
        while b!=0:
            temp=b
            b=a%b
            a=temp

        interp=int(self.Rs_bb/a)
        decim=int(self.Rs_sim/a)

        self.print_log(type="I", msg=f"Interp factor is: {interp}")
        self.print_log(type="I", msg=f"Decim factor is: {decim}")
        return interp, decim


    



    def demodulate(self):
        """ Method for selecting and calling wanted modulation type
            Example
            -------
            self.demodulate()

        """

 
        #t_vec=self.s[:,0]
        #Fs=1/t_vec[1]
        Fs=self.Rs_sim
        self.print_log(type='I', msg=f"Sampling Frequency Rs_sim used in Remez: {self.Rs_sim}")
        t_vec=np.arange(0,len(self.s))/Fs
        #pdb.set_trace()
        order=100 
        b=sig.remez(int(order+1),[0,Fs*0.24,Fs*0.28,0.5*Fs],[1,0],Hz=int(Fs))
        #w, h = sig.freqz(b, [1], worN=2000)
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        #ax.plot(0.5*Fs*w/np.pi, 20*np.log10(np.abs(h)))
        #ax.set_xlim(0, 0.5*Fs)
        #plt.show(block=False)
        #b=np.pad(b,(0,len(s)-len(b)),constant_values=0).reshape((1,-1))
        #s_fil=np.convolve(s.reshape((-1,1))[:,0],b,mode='same').reshape((-1,1))
        #s_fil=np.fft.ifft(np.fft.fft(s.reshape((-1,1))[:,0],len(s))*np.fft.fft(b,len(s))).reshape((-1,1))
        #s_raw=s
        #pdb.set_trace()
        #I=self.s[:,1]*np.cos(2*np.pi*t_vec*self.Fc)
        #Q=self.s[:,1]*np.sin(2*np.pi*t_vec*self.Fc)
        self.print_log(type='I', msg=f"Carrier frequency within receiver (demod): {self.Fc}")
        I=self.s*np.cos(2*np.pi*t_vec*self.Fc)
        Q=self.s*np.sin(2*np.pi*t_vec*self.Fc)

        #I=np.fft.ifft(np.fft.fft(I.reshape((-1,1))[:,0],len(I))*np.fft.fft(b,len(I))).reshape((-1,1))
        #Q=np.fft.ifft(np.fft.fft(Q.reshape((-1,1))[:,0],len(Q))*np.fft.fft(b,len(Q))).reshape((-1,1))
        I=sig.convolve(I,b,'same' ).reshape((-1,1)) 
        Q=sig.convolve(Q,b,'same').reshape((-1,1)) 

        s_out=np.transpose(np.vstack((t_vec, np.real(I[:,0]), np.real(Q[:,0]))))



        return s_out


    def resample(self):
        """ Method for interpolation and decimation of signal based on factors.
            Example
            -------
            self.interpolate()

        """
        #pdb.set_trace()
        Fs_org=self.Rs_sim

        s_fil=(self.s[:,1]+1j*self.s[:,2]).reshape((-1,1))
        #s_fil=self.s.reshape((-1,1))
        #pdb.set_trace()
        #test_decim=sig.decimate(s_fil,self.decim,ftype='fir')
        interp=self.interp
        decim=self.decim
        decim2=decim
        if interp >1:
            if not np.log2(interp).is_integer():
                t=np.arange(0,len(s_fil))/Fs
                s_out=np.transpose(np.vstack((np.real(s_fil[:,0]),np.imag(s_fil[:,0]))))
                s_out=np.transpose(np.vstack((t,s_out[:,0],s_out[:,1])))
                return s_out
            fact=int(np.log2(interp)/1)
            Fs=Fs_org*2
            order=70*fact/2
            #pdb.set_trace() 
            for i in range(0,int(fact)):
                s=np.transpose(np.matrix(np.array(s_fil).reshape(1,-1)))
                s=np.pad(s,((0,0),(0,int(1))),'constant',constant_values=0)
                s=np.concatenate(s)
                s=np.trim_zeros(s)
                
                Fcut=Fs/2
                #b=sig.firwin(int(order+1),Fs/4,fs=Fs)
                b=sig.remez(int(order+1),[0,Fs*0.24,Fs*0.28,0.5*Fs],[1,0],Hz=int(Fs))
                #w, h = sig.freqz(b, [1], worN=2000)
                #fig = plt.figure()
                #ax = fig.add_subplot(111)
                #ax.plot(0.5*Fs*w/np.pi, 20*np.log10(np.abs(h)))
                #ax.set_xlim(0, 0.5*Fs)
                #plt.show(block=False)
                #b=np.pad(b,(0,len(s)-len(b)),constant_values=0).reshape((1,-1))
                #s_fil=np.convolve(s.reshape((-1,1))[:,0],b,mode='same').reshape((-1,1))
                s_fil=np.fft.ifft(np.fft.fft(s.reshape((-1,1))[:,0],len(s))*np.fft.fft(b,len(s))).reshape((-1,1))

                #s_fil=s_fil[0::2,0].reshape((-1,1))
                # compensate group delay (order/2)
                s_fil=np.concatenate((s_fil[int(order/2):],s_fil[0:int(order/2)]))
                #s_out=s_fil
                self.s=s_fil
                Fs=Fs*2

     
        if decim>1:
            fact=int(np.log2(decim)/1)
            Fs=self.Rs_sim
            order=100
            #pdb.set_trace()
            while decim>10:
                Fs=Fs/2
                decim=decim/2
                b=sig.remez(int(order+1),[0,Fs*0.24,Fs*0.28,0.5*Fs],[1,0],Hz=int(Fs))
                #w, h = sig.freqz(b, [1], worN=2000)
                #fig = plt.figure()
                #ax = fig.add_subplot(111)
                #ax.plot(0.5*Fs*w/np.pi, 20*np.log10(np.abs(h)))
                #ax.set_xlim(0, 0.5*Fs)
                #plt.show(block=False)
                #b=np.pad(b,(0,len(s_fil)-len(b)),constant_values=0)
                #self.print_log(type='I', msg='Doing filtering with convolution...')
                s_fil=np.convolve(s_fil.reshape((-1,1))[:,0],b,mode='same').reshape((-1,1))
               
                self.print_log(type='I', msg='Doing filtering with FFT...')
                #s_fil=np.fft.ifft(np.fft.fft(s_fil.reshape((-1,1))[:,0],len(s_fil))*np.fft.fft(b,len(s_fil))).reshape((-1,1))

                #s_fil=s_fil[0::2,0].reshape((-1,1))
                # compensate group delay (order/2)
                #s_fil=np.concatenate((s_fil[int(order/2):],s_fil[0:int(order/2)]))
                #s_out=s_fil
                s_fil=s_fil[0::2]
                #Fs=Fs/2
            # Decimat by remaining factor
            #pdb.set_trace()
            decim=int(decim)
            Fs=Fs/decim
            b=sig.remez(int(order+1),[0,Fs*0.24,Fs*0.28,0.5*Fs],[1,0],Hz=int(Fs))
            #w, h = sig.freqz(b, [1], worN=2000)
            #fig = plt.figure()
            #ax = fig.add_subplot(111)
            #ax.plot(0.5*Fs*w/np.pi, 20*np.log10(np.abs(h)))
            #ax.set_xlim(0, 0.5*Fs)
            #plt.show(block=False)
            #b=np.pad(b,(0,len(s_fil)-len(b)),constant_values=0)
            #self.print_log(type='I', msg='Doing filtering with convolution...')
            s_fil=np.convolve(s_fil.reshape((-1,1))[:,0],b,mode='same').reshape((-1,1))
            
            #self.print_log(type='I', msg='Doing filtering with FFT...')
            #s_fil=np.fft.ifft(np.fft.fft(s_fil.reshape((-1,1))[:,0],len(s_fil))*np.fft.fft(b,len(s_fil))).reshape((-1,1))

            #s_fil=s_fil[0::2,0].reshape((-1,1))
            # compensate group delay (order/2)
            #s_fil=np.concatenate((s_fil[int(order/2):],s_fil[0:int(order/2)]))
            #s_out=s_fil
            s_fil=s_fil[0::decim]




        Fs=self.Rs_sim*interp/decim2
        t=np.arange(0,len(s_fil))/Fs
        self.fs_test=Fs
        s_out=np.transpose(np.vstack((np.real(s_fil[:,0]),np.imag(s_fil[:,0]))))
        # Normalize to 1
        test_matrix=np.concatenate((np.absolute(np.real(s_out)),np.absolute(np.imag(s_out))))
        max_value=test_matrix.max()
        s_out=s_out/max_value

        s_out=np.transpose(np.vstack((t,s_out[:,0],s_out[:,1])))

    
        return s_out











    def interpolate(self):
        """ Method for interpolation and decimation of signal based on factors.
            Example
            -------
            self.interpolate()

        """
        #pdb.set_trace()
        additional_factor=self.IOS.Members['additional_interp_factor'].Data
        #additional_factor=None
        if additional_factor is not None:
            #pdb.set_trace()
            s_fil=(self.s[:,1]+1j*self.s[:,2]).reshape((-1,1))
            #order=np.floor(70*additional_factor/40)
            order=70
            Fs=1/self.t[1]
            F_cut=Fs/(additional_factor)
            b=sig.remez(int(order+1),[0,F_cut*0.24,F_cut*0.28,0.5*F_cut],[1,0],Hz=int(F_cut))
            w, h = sig.freqz(b, [1], worN=2000)
            fig =plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(0.5*Fs*w/np.pi, 20*np.log10(np.abs(h)))
            ax.set_xlim(0, 0.5*Fs)
            plt.show(block=False)
            #b=np.pad(b,(0,len(s_fil)-len(b)),constant_values=0)
            s_fil=np.convolve(s_fil.reshape((-1,1))[:,0],b,mode='same').reshape((-1,1))
      
            #s_fil=np.fft.ifft(np.fft.fft(s_fil.reshape((-1,1))[:,0],len(s_fil))*np.fft.fft(b,len(s_fil))).reshape((-1,1))

            #s_fil=s_fil[0::2,0].reshape((-1,1))
            # compensate group delay (order/2)
            s_fil=np.concatenate((s_fil[int(order/2):],s_fil[0:int(order/2)]))
            #s_out=s_fil
            s_fil=s_fil[0::additional_factor]
            

            Fs_org=(1/self.t[1])/additional_factor
        else:
            Fs_org=1/self.t[1]

            s_fil=(self.s[:,1]+1j*self.s[:,2]).reshape((-1,1))
        #s_fil=self.s.reshape((-1,1))
        #pdb.set_trace()
        interp=self.interp
        decim=self.decim
        if interp >1:
            if not np.log2(interp).is_integer():
                t=np.arange(0,len(s_fil))/Fs
                s_out=np.transpose(np.vstack((np.real(s_fil[:,0]),np.imag(s_fil[:,0]))))
                s_out=np.transpose(np.vstack((t,s_out[:,0],s_out[:,1])))
                return s_out
            fact=int(np.log2(interp)/1)
            Fs=Fs_org*2
            order=70*fact/2
            #pdb.set_trace() 
            for i in range(0,int(fact)):
                s=np.transpose(np.matrix(np.array(s_fil).reshape(1,-1)))
                s=np.pad(s,((0,0),(0,int(1))),'constant',constant_values=0)
                s=np.concatenate(s)
                s=np.trim_zeros(s)
                
                Fcut=Fs/2
                #b=sig.firwin(int(order+1),Fs/4,fs=Fs)
                b=sig.remez(int(order+1),[0,Fs*0.24,Fs*0.28,0.5*Fs],[1,0],Hz=int(Fs))
                #w, h = sig.freqz(b, [1], worN=2000)
                #fig = plt.figure()
                #ax = fig.add_subplot(111)
                #ax.plot(0.5*Fs*w/np.pi, 20*np.log10(np.abs(h)))
                #ax.set_xlim(0, 0.5*Fs)
                #plt.show(block=False)
                #b=np.pad(b,(0,len(s)-len(b)),constant_values=0).reshape((1,-1))
                #s_fil=np.convolve(s.reshape((-1,1))[:,0],b,mode='same').reshape((-1,1))
                s_fil=np.fft.ifft(np.fft.fft(s.reshape((-1,1))[:,0],len(s))*np.fft.fft(b,len(s))).reshape((-1,1))

                #s_fil=s_fil[0::2,0].reshape((-1,1))
                # compensate group delay (order/2)
                s_fil=np.concatenate((s_fil[int(order/2):],s_fil[0:int(order/2)]))
                #s_out=s_fil
                self.s=s_fil
                Fs=Fs*2

     
        if decim>1:
            if not np.log2(decim).is_integer():
                t=np.arange(0,len(s_fil))/Fs
                s_out=np.transpose(np.vstack((np.real(s_fil[:,0]),np.imag(s_fil[:,0]))))
                s_out=np.transpose(np.vstack((t,s_out[:,0],s_out[:,1])))
                return s_out
            fact=int(np.log2(decim)/1)
            Fs=Fs_org/2
            order=70*fact/2
            for i in range(0,fact):

                b=sig.remez(int(order+1),[0,Fs*0.24,Fs*0.28,0.5*Fs],[1,0],Hz=int(Fs))
                #w, h = sig.freqz(b, [1], worN=2000)
                #fig = plt.figure()
                #ax = fig.add_subplot(111)
                #ax.plot(0.5*Fs*w/np.pi, 20*np.log10(np.abs(h)))
                #ax.set_xlim(0, 0.5*Fs)
                #plt.show(block=False)
                #b=np.pad(b,(0,len(s_fil)-len(b)),constant_values=0)
                #self.print_log(type='I', msg='Doing filtering with convolution...')
                #s_fil=np.convolve(s_fil.reshape((-1,1))[:,0],b,mode='same').reshape((-1,1))
                
                self.print_log(type='I', msg='Doing filtering with FFT...')
                s_fil=np.fft.ifft(np.fft.fft(s_fil.reshape((-1,1))[:,0],len(s_fil))*np.fft.fft(b,len(s_fil))).reshape((-1,1))

                #s_fil=s_fil[0::2,0].reshape((-1,1))
                # compensate group delay (order/2)
                s_fil=np.concatenate((s_fil[int(order/2):],s_fil[0:int(order/2)]))
                #s_out=s_fil
                s_fil=s_fil[0::2]
                Fs=Fs/2


        Fs=Fs_org*interp/decim
        t=np.arange(0,len(s_fil))/Fs
        self.fs_test=Fs
        s_out=np.transpose(np.vstack((np.real(s_fil[:,0]),np.imag(s_fil[:,0]))))
        # Normalize to 1
        test_matrix=np.concatenate((np.absolute(np.real(s_out)),np.absolute(np.imag(s_out))))
        max_value=test_matrix.max()
        s_out=s_out/max_value

        s_out=np.transpose(np.vstack((t,s_out[:,0],s_out[:,1])))

    
        return s_out

    def quant(self,**kwargs):
        """ Method for signal quantization.
            Parameters
            ----------
            sig : array
                Signal
            bits : integer
                Number of bits used
            low : integer
                Lower limit for signal value
            high : integer
                Upper limit for signal value
            p : binary (0,1)
                1 if signal is phase and values are 2pi periodic


            Example
            -------
            self.quant(sig=phase,bits=10,low=0,high=2*pi,p=1)

        """
        sig=self.s[:,1]
        bits=10
        low=0
        high=1
        p=0
        if 'sig' in kwargs:
            sig=kwargs.get('sig')
        if 'bits' in kwargs:
            bits=kwargs.get('bits')
        if 'low' in kwargs:
            low=kwargs.get('low')
        if 'high' in kwargs:
            high=kwargs.get('high')
        if 'p' in kwargs:
            p=kwargs.get('p')

        #pdb.set_trace()
        L=2**bits
        R=2
        d=R/L
        if p==0:
            levels=np.linspace(low,high,L)
        else:
            pos=sig>0
            sig=np.mod(sig,2*np.pi)
            sig[(sig==0) & (pos)]=2*np.pi
            #sig[sig<0]+=2*np.pi
            levels=np.linspace(low,high,L)

        # With roundings
        centers=(levels[1:]+levels[:-1])/2
        s_index=np.digitize(sig,centers)

        new_s=levels[s_index]

        # Without rounding => Next higher
        #s_index=np.digitize(sig,levels)
        #new_s=levels[s_index-1]

        s_out=np.transpose(np.vstack(new_s))

        return s_out
    

    def save_figs(self,**kwargs):
        figs=kwargs.get('figs',[plt.figure(n) for n in plt.get_fignums()])
        over_write=kwargs.get('over_write',True)
        pdb.set_trace()
        folder=kwargs.get('folder',"../figs")# self.picpath)#"../figs")
        if os.path.exists(folder) and os.path.isdir(folder):
            if over_write:
                self.print_log(type='I', msg='Over writing will delete ALL figures greated before. Press Enter to continu or stop simulation.')
                input()
                shutil.rmtree(folder)
                os.makedirs(folder)
            else:
                pass
        else:
            os.makedirs(folder)
        figs=[plt.figure(n) for n in plt.get_fignums()]
        for fig in figs:
            time=str(datetime.now())
            if len(fig.axes)>=1:
                title=fig.axes[0].get_title()
                fig.savefig(folder+"/"+title+time+".eps", format='eps')









    def define_io_conditions(self):
        '''This overloads the method is called by run_rtl method. It defines the read/write conditions for the files

        '''
        # Input A is read to verilog simulation after 'initdone' is set to 1 by controller
        self.iofile_bundle.Members['A'].verilog_io_condition='initdone'
        # Output is read to verilog simulation when all of the outputs are valid, 
        # and after 'initdone' is set to 1 by controller
        self.iofile_bundle.Members['Z'].verilog_io_condition_append(cond='&& initdone')






if __name__=="__main__":

    # SelfTest DO NOT work for now
    import matplotlib.pyplot as plt
    from  ht_chip import *
    from  ht_chip.controller import controller as ht_chip_controller
    import pdb
    length=1024
    rs=100e6
    Fs=1500
    Fc=5e6
    t=np.arange(0,length)/1000
    indata=np.random.randint(2,size=(length,3)).astype(float)#.reshape(-1,1);
    indata[:,0]=t
    #indata=np.random.randint(2,size=length)
    controller=ht_chip_controller()
    controller.Rs=rs
    #controller.reset()
    #controller.step_time()
    controller.start_datafeed()
    #pdb.set_trace()
    models=[ 'py']
    duts=[]
    BWP = np.array([[[3,7,0,1]]])
    QAM = "64QAM"
    osr = 1
    Fc =1*10**9
    time_resolution_multiplier =10
    factor = int(2**4/1)
    BW = np.array([200e6])
    #in_bits = np.array([[np.random.randint(2,size=4752*6)]])
    in_bits = np.array([["max"]])

    am_resolution_tran=10
    phase_resolution_tran=10
    am_resolution_rec=10
    phase_resolution_rec=10
    modulation='multi_outphasing'  # cartesian , polar , outphasing , multi_outphasing 
    #loop = True
    new_model = True

    for model in models:
        gen=NR_signal_generator()
        chip_transmitter=ht_chip()
        chip_transmitter.new_model = new_model
        chip_transmitter.init(Fc=Fc,modulation=modulation,factor=factor,am_resolution=am_resolution_tran,phase_resolution=phase_resolution_tran, mode='tran' )

        d=receiver()
        duts.append(d) 
        d.model=model
        #d.Rs=rs
        #d.preserve_iofiles=True
        #d.interactive_rtl=True
        #d.interactive_spice=True



        d.init_params(Fc=-Fc,modulation=modulation,factor=(1/factor),am_resolution=am_resolution_rec,phase_resolution=phase_resolution_rec,mode='rec' )
        chip_transmitter.IOS.Members['bb_data_in']=gen.IOS.Members['out']
        d.IOS.Members['sig_in']=chip_transmitter.IOS.Members['rf_out']
        gen.IOS.Members['in_dem']=d.IOS.Members['out'] 


        gen.run_gen()
        chip_transmitter.run()

        #pdb.set_trace()
        d.Rs_bb=gen.s_struct["Fs"]
        d.Rs_sim=chip_transmitter.Rs_sim

        d.run()
        print("Done receiving")
        gen.run_dem()
        gen.run_EVM()
        print(gen.EVM)
        s=gen.s_struct["s"]
        #:pdb.set_trace()
        plot_func.plot_PSD(signal=d.s_input,Fc=Fc,Fs=d.Rs_sim, BW=BW)
                   
        plot_func.plot_PSD(BW=BW,signal=d.s,Fc=0,Fs=d.Rs_bb)
        plot_func.plot_PSD(BW=BW,signal=s,Fc=1e6,Fs=d.Rs_bb)

        d.save_figs()
        input()
        
