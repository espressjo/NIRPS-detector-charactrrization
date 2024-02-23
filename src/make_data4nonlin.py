#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 08:24:58 2021

@author: noboru
"""

import numpy as np
from astropy.io import fits
from os.path import join,isfile
from tqdm import tqdm
from hxrg.fits2ramp.refpixcorr import refpxcorr as rpc
from datetime import datetime
from getconfig import getString

flats= ["NIRPS_2022-05-05T16_09_12_198.fits",
        "NIRPS_2022-05-05T16_30_17_336.fits",
        "NIRPS_2022-05-05T16_51_22_473.fits",
        "NIRPS_2022-05-05T17_12_27_612.fits",
        "NIRPS_2022-05-05T17_33_32_752.fits",
        "NIRPS_2022-05-05T18_15_43_030.fits",
        "NIRPS_2022-05-05T18_36_48_170.fits",
        "NIRPS_2022-05-05T19_18_58_450.fits",
        "NIRPS_2022-05-05T19_40_03_591.fits",
        "NIRPS_2022-05-05T20_01_08_734.fits",
        "NIRPS_2022-05-05T20_22_13_876.fits",
        "NIRPS_2022-05-05T20_43_19_017.fits",
        "NIRPS_2022-05-05T21_04_18_580.fits",
        "NIRPS_2022-05-05T21_25_29_302.fits",
        "NIRPS_2022-05-05T22_07_45_164.fits",
        "NIRPS_2022-05-05T22_49_55_447.fits"]

output = getString('OUTPUTPATH',default='') 
tmp = getString('TMPPATH',default='')
path = getString('DATAPATH',default='')

if not isfile(join(output,"superbias.fits")):
    print("You must run superbias.py first.")
    exit(0)

sbias_file = join(output,"superbias.fits")
#p_data = getString('DATAPATH',default="/nirps_raw/nirps/reads")

#load the bias file
hdu = fits.open(join(path,sbias_file))
superbias = np.asarray(hdu[0].data,dtype=float)
superbias_uid = hdu[0].header['UNIQUEID']
hdu.close()

def _with_bias():
    print("Making data")
    for fname in tqdm(flats):
        data = fits.getdata(join(path,fname))
        uid = fits.getheader(join(path,fname))["HIERARCH ESO DET RAMP ID"]
        cube = np.zeros(((data.shape[0]),4096,4096))
        for i in range((cube.shape[0])):
            cube[i,:,:] = rpc(np.asarray(data[i],dtype=float)-superbias).refpxcorr()
        hdu = fits.PrimaryHDU(data = cube)
        time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
        ID = datetime.now().strftime("%y%m%d%H%M%S")
        hdu.header['DATE'] = (time,"Creation date")
        hdu.header['UNIQUEID'] = (ID,"Unique id of this file")
        hdu.header['RUID'] = (uid,"Unique id of the ramp used to built this file")
        hdu.header['SBIASID'] = (superbias_uid,"Unique id of the super bias")
        hdu.header['SBIASF'] = (sbias_file,"Super bias filename")
        
        #save the files
        hdu.writeto(join(tmp,uid+".bias.fits"),overwrite=True)
        print("%s saved."%join(path,uid+".fits"))
'''
def _no_bias():
    print("Making data")
    for fname in tqdm(flats):
        data = fits.getdata(join(p_data,fname))
        uid = fits.getheader(join(p_data,fname))["HIERARCH ESO DET RAMP ID"]
        
        bias = np.asarray(data[0],dtype=float)
        
        cube = np.zeros((data.shape[0],4096,4096))
        
        for i in range(data.shape[0]-1):
            
            cube[i,:,:] = rpc(np.asarray(data[i+1],dtype=float)-bias).refpxcorr()
        hdu = fits.PrimaryHDU(data = cube)
        time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
        ID = datetime.now().strftime("%y%m%d%H%M%S")
        hdu.header['DATE'] = (time,"Creation date")
        hdu.header['UNIQUEID'] = (ID,"Unique id of this file")
        hdu.header['RUID'] = (uid,"Unique id of the ramp used to built this file")

        
        #save the files
        hdu.writeto(join(path,uid+".fits"))
        print("%s saved."%join(path,uid+".fits"))
'''
if '__main__' in __name__:
    
    _with_bias()

        
