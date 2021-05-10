#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 09:43:00 2021

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

class badpixelmap():
    def __init__(self):
        self.path = '/nirps_raw/characterization'
        self.pathtmp = '/nirps_raw/characterization/.tmp'
        self.nonlin = 'nonlin.fits'
        self.time_fmt = "%Y-%m-%dT%H:%M:%S"
        self.ls_dark = ['20210403135159',
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
        self.ls_flat = ['20210403171440',
                        '20210403175352',
                        '20210403183253',
                        '20210403191153',
                        '20210403195054',
                        '20210403203006',
                        '20210403210907',
                        '20210403214818',
                        '20210403222719',
                        '20210403230620']#unique ID for flat
        
        if not isfile(join(self.path,self.nonlin)):
            print("no nonlin.fits found.")
            self.nonlin = ''
    def _make_cds(self,fr1,fr2):
        
        return rpc(np.asarray(fits.getdata(fr2),dtype=float)-np.asarray(fits.getdata(fr1),dtype=float)).refpxcorr()
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
        hdu.writeto(join(self.pathtmp,'nl-bpm.fits'))
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
if __name__=='__main__':
    
       bpm = badpixelmap()     
       print("Computing hot pixels")
       bpm.hotpixels()
       print("Computing bad pixels")
       bpm.badpixels()
       print("Creating nl bpm")
       bpm.nonlinbpm()
       
       bpm.gpm()
       
       
       
    