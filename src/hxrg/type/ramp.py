#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 12:03:14 2020

@author: espressjo
"""
from numpy import zeros,int16,full,where,asarray,linspace,abs,cumsum,fliplr
from scipy import fft
from astropy.io import fits
from hxrg.fits2ramp.refpixcorr import refpxcorr
from hxrg.type.read import hxread
from os.path import join,isfile
from os import getcwd


class hxramp2():
    def __init__(self,**kwargs):
        
        self.i = 0
        self.H = None
        self.timestamp = [];
        self.inttime = 5.24288
        self.bias = zeros((4096,4096))
        if 'saturation' in kwargs:
            self.saturation = kwargs.get('saturation')
        else:
            self.saturation = 45000
        
        self.refpxeachread = False
        self.selfbias = True
        self.toponly = False
        self.notopnobottom = False
        self.oddeven = False
        self.nl = False
        if 'toponly' in kwargs:
            self.toponly = kwargs.get('toponly')
        if 'notopnobottom' in kwargs:
            self.notopnobottom = kwargs.get('notopnobottom')
        if 'oddeven' in kwargs:
            self.oddeven = kwargs.get('oddeven')
        if 'refpxeachread' in kwargs:
            self.refpxeachread = kwargs.get('refpxeachread')

    def nlcorr(self,array):
        array = asarray(array,dtype=float)
        
        return array+self.c3*(array**3)+self.c2*(array**2)
    def upload_nonlin(self,f:str):
        if isfile(f):
            hdu = fits.open(f)            
            self.c3 = hdu[0].data
            self.c2 = hdu[1].data
            self.nl = True
        else:
            print('%s not found, no nl-correction'%f)
    def upload_bias(self,fname:str):
        '''
        Upload a super bias
        '''
        if isfile(fname):
            self.bias = fits.getdata(fname)
            
            self.selfbias = False
        else:
            print("[Warning] bias file not found, will perform selfbias")            
    def unfold(self):
        '''
        Flip [l/r] every odd amp. This can be usefull to perform fft, 
        or IPC analysis. fit() method should be called before this.

        Returns
        -------
        None.

        '''
        for i in range(32):
            if (i+1)%2!=0:
                amp = self.a[:,i*128:(i+1)+128]
                self.a[:,i*128:(i+1)+128] = fliplr(amp)
                amp = self.b[:,i*128:(i+1)+128]
                self.b[:,i*128:(i+1)+128] = fliplr(amp)
                amp = self.n[:,i*128:(i+1)+128]
                self.n[:,i*128:(i+1)+128] = fliplr(amp)
    def __call__(self,fname:str):
        '''
        Open a ramp file.
        '''

        if '/' not in fname:
            fname = join(getcwd(),fname)
        hdu = fits.open(fname)
        self.H = hdu[0].header
        if 'HIERARCH ESO DET SEQ1 DIT' in self.H:
            self.total_integration = self.H['HIERARCH ESO DET SEQ1 DIT']
        self.a = asarray(hdu[1].data,dtype=float)
        self.b = asarray(hdu[2].data,dtype=float)
        self.n = hdu[3].data
    def init_from_read(self,read:hxread):
        '''
        Initialize a new ramp using a 1st read

        Parameters
        ----------
        read : hxread
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        self.shape = (read.x*1024,read.x*1024)
        self.sx = zeros(self.shape)
        self.sx2 = zeros(self.shape)
        self.n = zeros(self.shape,dtype=int16)#number of read used to perform fits2ramp
        self.sy = zeros(self.shape)
        self.sxy = zeros(self.shape)
        self.a = zeros(self.shape)#slope
        self.b = zeros(self.shape)#intercept
    def applyrefpxcorr(self,array):
        '''
        Apply ref pix. correction given the options

        Parameters
        ----------
        array : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        from hxrg.fits2ramp.refpixcorr import refpxcorr as rpc
        if self.toponly:
            return rpc(asarray(array,dtype=float)).refpxcorrtop(oddeven=self.oddeven)
        elif self.notopnobottom:
            return rpc(asarray(array,dtype=float)).refpxcorrside(oddeven=self.oddeven)
        else:
            return rpc(asarray(array,dtype=float)).refpxcorr(oddeven=self.oddeven)
    def get_time(self):
        return self.timestamp[:-1]
    def get_coeff(self):
        '''
        Return slope and intercept as hxread datatype

        Returns
        -------
        a : hxread
            slope.
        b : hxread
            intercept.

        '''
        a = hxread()
        a(self.a)
        a.H = self.H
        a.fname = 'coeff_a.fits'
        b = hxread()
        b(self.b)
        b.H = self.H        
        b.fname = 'coeff_b.fits'
        return a,b
    def write_to(self,fname):
        hdul = fits.PrimaryHDU([self.a,self.b,self.n])
        hdul.writeto(fname,overwrite=True)
    '''    
    def __lshift__(self, o:hxread):        
        if self.i==0:
            self.init_from_read(o)
            if self.selfbias:
                self.bias = o.im
            self.timestamp.append(self.inttime)
            self.i+=1
            self.total_integration = 0
            self.goodmask = full((o.w,o.h),True,dtype=bool)
            self.n+=self.goodmask
            self.H = o.H
            return#???
        self.total_integration +=self.inttime
        self.timestamp.append(self.inttime*(self.i+1))
        self.goodmask = (o.im <= self.saturation)*self.goodmask
        self.n+=self.goodmask    
        if self.nl:
            imc = self.nlcorr(o.im-self.bias)
        else:
            imc = (o.im-self.bias)
        if self.refpxeachread:
             imc = self.applyrefpxcorr(imc)
        
        self.sy[self.goodmask]+=imc[self.goodmask]
        self.sxy[self.goodmask]+=(imc[self.goodmask]*self.inttime*(self.i+1))
        self.i+=1
    '''
    def __lshift__(self, o:hxread):        
        if self.i==0:
            self.init_from_read(o)
            if self.selfbias:
                self.bias = o.im
            #self.timestamp.append(self.inttime)
            #self.i+=1
            self.total_integration = 0
            self.goodmask = full((o.w,o.h),True,dtype=bool)
            #self.n+=self.goodmask
            self.H = o.H
        self.i+=1
        self.total_integration +=self.inttime
        self.timestamp.append(self.inttime*(self.i))
        self.goodmask = (o.im <= self.saturation)*self.goodmask
        self.n+=self.goodmask    
        if self.nl:
            
            imc = self.nlcorr(o.im-self.bias)
        else:
            imc = (o.im-self.bias)
        if self.refpxeachread:
             imc = self.applyrefpxcorr(imc)
        
        self.sy[self.goodmask]+=imc[self.goodmask]
        self.sxy[self.goodmask]+=(imc[self.goodmask]*self.inttime*(self.i))
        

    def fit(self):
        self.sx=where(self.n>0,(cumsum(self.timestamp))[self.n-1],0)
        self.sx2=where(self.n>0,(cumsum(asarray(self.timestamp)**2))[self.n-1],0)
        valid = self.n>1
        self.b[valid] = (self.sx*self.sxy-self.sx2*self.sy)[valid]/(self.sx**2-self.n*self.sx2)[valid] # algebra of the linear fit
        self.a[valid] = (self.sy-self.n*self.b)[valid]/self.sx[valid]        
        if not self.refpxeachread:
            self.a = self.applyrefpxcorr(self.a)
            
        if self.refpxeachread:
            print("Ref. px. corrected after each read")
        else:
            print("Ref. px. corrected on the slope")
        if self.selfbias:
            print("Self biased")
        else:
            print("super bias used")
        
        if self.toponly:
            print("Only using top and side ref. px.")
        else:
            print("All ref. px. used")
        if self.notopnobottom:
            print("Only side ref. px. used")
        if self.oddeven:
            print("odd and even column processed separetly")
   
        if self.nl:
            print("nl correction applied.")
        self.b+=self.bias            
    def show(self):
        from matplotlib import pyplot as plt
        from astropy.visualization import (SqrtStretch,ImageNormalize,ZScaleInterval)
        f,ax = plt.subplots()
        ax.set_title('Coeff. a')
        norm = ImageNormalize(self.a,interval=ZScaleInterval(),stretch=SqrtStretch()) 
        ax.imshow(self.a,cmap='Greys_r',norm=norm)   
        
        f,ax = plt.subplots()
        ax.set_title('Coeff. b')
        norm = ImageNormalize(self.b,interval=ZScaleInterval(),stretch=SqrtStretch()) 
        ax.imshow(self.b,cmap='Greys_r',norm=norm) 
        plt.show()
class hxramp(refpxcorr):
    def __init__(self,rfpx = True,toponly = False,sideonly=False,refpxslope=True):
        #by default the reference pixels are used everytime an image is read. We can overide this with refpxslope flag to True
        #hxrg.__init__(self)
        self.i = 0
        self.refpxslope = refpxslope
        self.H = None
        self.rfpx = rfpx
        
        if not rfpx:
            print('[Warning] reference pixel will not be used')
        self.timestamp = [];
        self.saturation = 45000
        if toponly and sideonly:
            print("cannot initiate ramp with both toponly and sideonly. Pick one")
            exit(0)
        self.sideonly = sideonly
        self.toponly = toponly
    def unfold(self):
        #flip left-rright every odd amp
        for i in range(32):
            if (i+1)%2!=0:
                amp = self.a[:,i*128:(i+1)+128]
                self.a[:,i*128:(i+1)+128] = fliplr(amp)
                amp = self.b[:,i*128:(i+1)+128]
                self.b[:,i*128:(i+1)+128] = fliplr(amp)
                amp = self.n[:,i*128:(i+1)+128]
                self.n[:,i*128:(i+1)+128] = fliplr(amp)
    def __call__(self,e):
        if isinstance(e,str):
            if '/' not in e:
                f = join(getcwd(),e)
            f = e
            hdu = fits.open(f)
            self.H = hdu[0].header
            if 'HIERARCH ESO DET SEQ1 DIT' in self.H:
                self.total_integration = self.H['HIERARCH ESO DET SEQ1 DIT']
            self.a = asarray(hdu[1].data,dtype=float)
            self.b = asarray(hdu[2].data,dtype=float)
            self.n = hdu[3].data
    def init_from_read(self,read:hxread):
        self.shape = (read.x*1024,read.x*1024)
        self.sx = zeros(self.shape)
        self.sx2 = zeros(self.shape)
        self.n = zeros(self.shape,dtype=int16)#number of read used to perform fits2ramp
        self.sy = zeros(self.shape)
        self.sxy = zeros(self.shape)
        self.a = zeros(self.shape)#slope
        self.b = zeros(self.shape)#intercept
        
    def get_time(self):
        return self.timestamp[:-1]
    def get_coeff(self):
        a = hxread()
        a(self.a)
        a.H = self.H
        a.fname = 'coeff_a.fits'
        b = hxread()
        b(self.b)
        b.H = self.H
        
        b.fname = 'coeff_b.fits'
        return a,b
    
    def write_to(self,fname):
        hdul = fits.PrimaryHDU([self.a,self.b,self.n])
        hdul.writeto(fname,overwrite=True)
    def __lshift__(self, o:hxread):        
        if self.i==0:
            self.init_from_read(o)
            self.bias = hxread()
            self.bias(o.im)
            self.inttime = 5.24288#other.H['AEEXPT']
            self.timestamp.append(self.inttime)
            self.i+=1
            self.total_integration = 0
            self.goodmask = full((o.w,o.h),True,dtype=bool)
            self.n+=self.goodmask
            self.H = o.H
            return
        self.total_integration +=self.inttime
        self.timestamp.append(self.inttime*(self.i+1))
        self.goodmask = (o.im <= self.saturation)*self.goodmask
        self.n+=self.goodmask    
        imc = o-self.bias
        if self.rfpx:
            imc.refpx(toponly=self.toponly,sideonly=self.sideonly)
        self.sy[self.goodmask]+=imc.im[self.goodmask]
        self.sxy[self.goodmask]+=(imc.im[self.goodmask]*self.inttime*(self.i+1))
        self.i+=1        
    def fit(self):
        self.sx=where(self.n>0,(cumsum(self.timestamp))[self.n-1],0)
        self.sx2=where(self.n>0,(cumsum(asarray(self.timestamp)**2))[self.n-1],0)
        valid = self.n>1
        self.b[valid] = (self.sx*self.sxy-self.sx2*self.sy)[valid]/(self.sx**2-self.n*self.sx2)[valid] # algebra of the linear fit
        self.a[valid] = (self.sy-self.n*self.b)[valid]/self.sx[valid]
        
        self.b+=self.bias.im            
    def show(self):
        from matplotlib import pyplot as plt
        from astropy.visualization import (SqrtStretch,ImageNormalize,ZScaleInterval)
        f,ax = plt.subplots()
        ax.set_title('Coeff. a')
        norm = ImageNormalize(self.a,interval=ZScaleInterval(),stretch=SqrtStretch()) 
        ax.imshow(self.a,cmap='Greys_r',norm=norm)   
        
        f,ax = plt.subplots()
        ax.set_title('Coeff. b')
        norm = ImageNormalize(self.b,interval=ZScaleInterval(),stretch=SqrtStretch()) 
        ax.imshow(self.b,cmap='Greys_r',norm=norm) 
        plt.show()
if '__main__' in __name__:
    
    ramp = hxramp()
    ramp('/home/noboru/NIRPS_GEN_ORDERDEF253_0006.fits')
    