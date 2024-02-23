#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 11:00:56 2021

@author: noboru
"""
from os.path import isfile,join
from findread import findread
from findramp import findramp
from hxrg.fits2ramp.refpixcorr import refpxcorr as rpc
from astropy.io import fits
import numpy as np
from datetime import datetime
from astropy.stats import sigma_clipped_stats as sc
from scipy.ndimage import median_filter

class pipelinedBPM():
    def __init__(self):
        self.path = '/nirps_raw/characterization'
        self.pathtmp = '/nirps_raw/characterization/.tmp'
        self.nonlin = 'nonlin.fits'
        self.time_fmt = "%Y-%m-%dT%H:%M:%S"
        #data set of january 2 april 2021
        self.ls_dark_jan2april = [   '20210126142748',
                                     '20210126151555',
                                     '20210126151657',
                                     '20210126155614',
                                     '20210126154837',
                                     '20210126160627',
                                     '20210126160453',
                                     '20210126155749',
                                     '20210126145045',
                                     '20210126150244',
                                     '20210126151758',
                                     '20210126151859',
                                     '20210214144456',
                                     '20210214144738',
                                     '20210214154609',
                                     '20210214155616',
                                     '20210214160624',
                                     '20210214161631',
                                     '20210214162639',
                                     '20210314030404',
                                     '20210313185428',
                                     '20210314051727',
                                     '20210313203644',
                                     '20210314054817',
                                     '20210313204657',
                                     '20210314052740',
                                     '20210313205710',
                                     '20210314053753',
                                     '20210313210723',
                                     '20210313214821',
                                     '20210314063928']
        #data set for detetector characterization @ COPL
        
        self.ls_dark_april = [  '20210403135159',
                                '20210403140207',
                                '20210403141214',
                                '20210403142222',
                                '20210403143229',
                                '20210403144237',
                                '20210403145244',
                                '20210403150252',
                                '20210403151259',
                                '20210403152307',
                                '20210403153314',
                                '20210403154322',
                                '20210403155329',
                                '20210403160336',
                                '20210403161344',
                                '20210403162351',
                                '20210403163359',
                                '20210403164406',
                                '20210403165414',
                                '20210403170421']#unique ID
        
        #data set of january to april 2021
        self.ls_flat_jan2april = [  '20210217222559',
                                    '20210217223708',
                                    '20210217225947',
                                    '20210305112911',
                                    '20210305113344',
                                    '20210305114251',
                                    '20210305113818',
                                    '20210305114724']
        #data set for detetector characterization @ COPL
        
        self.ls_flat_april = [  '20210403171440',
                                '20210403175352',
                                '20210403183253',
                                '20210403191153',
                                '20210403195054',
                                '20210403203006',
                                '20210403210907',
                                '20210403214818',
                                '20210403222719',
                                '20210403230620']#unique ID for flat
        self.ls_flat = self.ls_flat_jan2april
        self.ls_dark = self.ls_dark_jan2april
        if not isfile(join(self.path,self.nonlin)):
            print("no nonlin.fits found.")
            self.nonlin = ''
    def _make_cds(self,fr1,fr2):
        '''
        For NIRPS we want to avoid using the bottom ref. pix. and 
        we do not want to use the odd-even option.
        Parameters
        ----------
        fr1 : STR
            File name of read 1
        fr2 : STR
            File name of read N
        Returns
        -------
        TYPE
            CDS np.array of float with corrected ref.px. using side and top.
        '''
        return rpc(np.asarray(fits.getdata(fr2),dtype=float)-np.asarray(fits.getdata(fr1),dtype=float)).refpxcorrtop(False)
    def hotpixels(self):
        fr = findread()
        cube = [];
        hp = np.zeros((4096,4096),dtype=int)
        for uid in self.ls_dark:
            ls = fr.findreads(uid)
            cube.append(self._make_cds(ls[0],ls[1]))
        dark = np.median(cube,axis=0)
        mean,_,std = sc(dark,sigma=5)
        hp[dark>(mean+5*std)]=1
        hdu = fits.PrimaryHDU(data=hp)
        hdu.header['DATE'] = (datetime.now().strftime(self.time_fmt),"Creation date")
        hdu.header['UNIQUEID'] = (datetime.now().strftime("%y%m%d%H%M%S"),"Unique identification number")
        hdu.header.add_comment("Built using:")
        for uid in self.ls_dark:
            hdu.header.add_comment(uid)
        hdu.writeto(join(self.pathtmp,"hotpixels.fits"),overwrite=True)
        print("hotpixels.fits write successfully.")
        return hp
    def badpixels(self):
        fr = findramp()
        freads = findread()
        #check slope on 1st ramp
        fname_ramp = fr.findramp(self.ls_flat[0])

        slope = fits.open(fname_ramp)['slope'].data
        mean,med,std = sc(slope[1800:2200,1800:2200],sigma=5)
        n = int(40000./(mean*5.5733))
        print("40,000ADU ~ n=%d"%n)
        cube_bpm = [];
        for uid in self.ls_flat:
            ls = freads.findreads(uid)
            if len(ls)-1 <n:
                cube_bpm.append(self._make_cds(ls[0], ls[-1]))     
            else:
                cube_bpm.append(self._make_cds(ls[0], ls[n]))
        im = np.median(cube_bpm,axis=0)
        filt =    median_filter(im,size=(20,20))
        filt/=np.median(filt)
        im/=filt
        mean,_,std = sc(im,sigma=5)
        bpm = np.zeros((4096,4096))
        bpm[im<(mean-5*std)]=1
        bpm[:4,:]=0
        bpm[-4:,:] = 0
        bpm[:,:4] = 0
        bpm[:,-4:] = 0
        
        _p = [128*i for i in range(32) if all([i!=0,i%2==0])]
        for p in _p:
            bpm[:,p-1]=0 
            bpm[:,p] = 0
            mn,_,st = sc(im[4:-4,p-1:p+1])
            _s1 = bpm[:,p-1:p+1]
            _s2 = im[:,p-1:p+1]
            _s1[_s2<(mn-5*st)]=1
            bpm[:,p-1:p+1] = _s1
        
        hdu = fits.PrimaryHDU(data=bpm)
        hdu.header['DATE'] = (datetime.now().strftime(self.time_fmt),"Creation date")
        hdu.header['UNIQUEID'] = (datetime.now().strftime("%y%m%d%H%M%S"),"Unique identification number")
        hdu.header.add_comment("Built using:")
        for uid in self.ls_flat:
            hdu.header.add_comment(uid)
        hdu.writeto(join(self.pathtmp,"badpixels.fits"),overwrite=True)
    def nonlinbpm(self):
        if self.nonlin=='':
            return
        nlbpm = fits.getdata(join(self.path,self.nonlin))[2]
        hdu = fits.PrimaryHDU(data=nlbpm)
        hdu.header['DATE'] = (datetime.now().strftime(self.time_fmt),"Creation date")
        hdu.header['UNIQUEID'] = (datetime.now().strftime("%y%m%d%H%M%S"),"Unique identification number")
        hdu.header.add_comment("Built using:")
        hdu.header.add_comment(fits.getheader(join(self.path,self.nonlin))['UNIQUEID'])
        hdu.writeto(join(self.pathtmp,'nl-bpm.fits'),overwrite=True)
    def gpm(self):
        gpm = np.ones((4096,4096),dtype=int)
        uid = [];
        if isfile(join(self.pathtmp,'hotpixels.fits')):
            print("Using hotpixels")
            hp = fits.getdata(join(self.pathtmp,'hotpixels.fits'))
            uid.append(fits.getheader(join(self.pathtmp,'hotpixels.fits'))['UNIQUEID'])
            gpm[hp==1]=0
        if isfile(join(self.pathtmp,'nl-bpm.fits')):
            print("Using nl-bpm")
            nl = fits.getdata(join(self.pathtmp,'nl-bpm.fits'))
            uid.append(fits.getheader(join(self.pathtmp,'nl-bpm.fits'))['UNIQUEID'])
            gpm[nl==1]=0
        if isfile(join(self.pathtmp,'badpixels.fits')):
            print("Using badpixels")
            bp = fits.getdata(join(self.pathtmp,'badpixels.fits'))
            uid.append(fits.getheader(join(self.pathtmp,'badpixels.fits'))['UNIQUEID'])
            gpm[bp==1]=0
        hdu = fits.PrimaryHDU(data=gpm)
        hdu.header['DATE'] = (datetime.now().strftime(self.time_fmt),"Creation date")
        hdu.header['UNIQUEID'] = (datetime.now().strftime("%y%m%d%H%M%S"),"Unique identification number")
        hdu.header.add_comment("Built using:")
        for ui in uid:
            hdu.header.add_comment(ui)
        hdu.writeto(join(self.path,'gpm.fits'),overwrite=True)
    def change_value(self,array,value=256):
        arr = np.zeros(array.shape)
        arr[array==1]=value
        return arr
        
    def make_eso_cmpl_file(self,withnl = False):
        hdu = fits.open(join(self.pathtmp,"badpixels.fits"))
        if withnl:
            nl = fits.getdata(join(self.path,'nonlin.fits'))
            hdu[0].data[nl==1]=1
        
        H = hdu[0].header
        H['HIERARCH ESO INS MODE'] = 'HA'
        H['HIERARCH ESO PRO CATG'] = 'BAD_PIXEL_MASK'
        H['HIERARCH ESO DET WIN1 BINX'] = 1 
        H['HIERARCH ESO DET WIN1 BINY'] = 1 
        H['HIERARCH ESO DET WIN1 NX'] = 4088 
        H['HIERARCH ESO DET WIN1 NY'] = 4088 
        H['MJD-OBS'] = 59335
        H['HIERARCH ESO QC EXT0 ROX0 ROY0 CONAD'] = 1/1.27
        #bad pixel map en HA
        Phdu = fits.PrimaryHDU(header=H)
        imext = fits.ImageHDU(self.change_value(np.flipud(np.rot90(hdu[0].data[4:-4,4:-4],-1)),value=8192))
        hdul = fits.HDUList([Phdu,imext])
        hdul.writeto(join(self.path,"badpixels_HA.fits"),overwrite=True)
        #bad pixel map en HE
        H['HIERARCH ESO INS MODE'] = 'HE'
        Phdu = fits.PrimaryHDU(header=H)
        imext = fits.ImageHDU(self.change_value(np.flipud(np.rot90(hdu[0].data[4:-4,4:-4],-1)),value=8192))
        hdul = fits.HDUList([Phdu,imext])
        hdul.writeto(join(self.path,"badpixels_HE.fits"),overwrite=True)
        
        
        
        #hot pixels
        hdu = fits.open(join(self.pathtmp,"hotpixels.fits"))
        H = hdu[0].header
        H['HIERARCH ESO INS MODE'] = 'HA'
        H['HIERARCH ESO PRO CATG'] = 'HOT_PIXEL_MASK'
        H['HIERARCH ESO DET WIN1 BINX'] = 1 
        H['HIERARCH ESO DET WIN1 BINY'] = 1 
        H['HIERARCH ESO DET WIN1 NX'] = 4088 
        H['HIERARCH ESO DET WIN1 NY'] = 4088 
        H['MJD-OBS'] = 59335
        
        H['HIERARCH ESO DET OUT1 GAIN'] = 1.27
        Phdu = fits.PrimaryHDU(header=H)
        imext = fits.ImageHDU(self.change_value(np.flipud(np.rot90(hdu[0].data[4:-4,4:-4],-1)),value=256))
        hdul = fits.HDUList([Phdu,imext])
        hdul.writeto(join(self.path,"hotpixels_HA.fits"),overwrite=True)
        #hot pixel map en HE
        H['HIERARCH ESO INS MODE'] = 'HE'
        Phdu = fits.PrimaryHDU(header=H)
        
        imext = fits.ImageHDU(self.change_value(np.flipud(np.rot90(hdu[0].data[4:-4,4:-4],-1)),value=256))
        hdul = fits.HDUList([Phdu,imext])
        hdul.writeto(join(self.path,"hotpixels_HE.fits"),overwrite=True)
        print("done")
        
if __name__=='__main__':
    
    bpm = pipelinedBPM()     
    print("Computing hot pixels")
    #bpm.hotpixels()
    print("Computing bad pixels")
    bpm.badpixels()
    print("Creating nl bpm")
    #bpm.nonlinbpm()
    bpm.make_eso_cmpl_file()
    bpm.gpm()
    
