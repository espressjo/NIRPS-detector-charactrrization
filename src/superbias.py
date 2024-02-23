#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 15:40:33 2021

@author: jonathan st-antoine

Since we don't want to use the bottom reference pixels, 
we need to recreate the ramp file to compute the superbias
"""
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import seaborn as sns
from os.path import basename,isfile,join
from tqdm import tqdm
from getconfig import getString
from hxrg.type.read import hxread
from hxrg.type.ramp import hxramp
from datetime import datetime

sns.set_theme()


class superbias():
    def __init__(self):
        self.path = getString('OUTPUTPATH',default='/nirps_raw/nirps/characterization') #place where to store the superbias file and the graphs
        self.tmppath = getString('TMPPATH',default='/nirps_raw/nirps/characterization/.tmp')#tmp path 
        self.p_data = getString('DATAPATH',default="/nirps_raw/nirps/reads")
        
    def makeramp(self,fname):
        if '/' not in fname:
            fname = join(self.p_data,fname)
        #given an uid, it will return a hxramp object
        cube = fits.getdata(fname)
        H = fits.getheader(fname)
        read = hxread()
        Ramp = hxramp(toponly = True)
        print("Working on %s"%basename(fname))
        for im in tqdm(cube):
            read(np.asarray(im,dtype=float))
            read.H = H
            Ramp<<read
        Ramp.fit()
        return Ramp
    def extractbias(self,darks):
        biases = [self.makeramp(fname).b for fname in darks]
        return np.asarray(biases)
            


def side(im):
    '''
    Return the median of side ref. pixels of an image

    '''
    l = np.append(im[:,0:4].ravel(),im[:,-4:].ravel())
    return np.median(l)

if '__main__' in __name__:
    
    #data to be used;
    darks = ["NIRPS_2022-05-05T15_31_01_590.fits",
            "NIRPS_2022-05-05T15_40_24_489.fits",
            "NIRPS_2022-05-05T15_49_47_388.fits",
            "NIRPS_2022-05-05T15_59_10_285.fits"]
    
    sb = superbias()
    
    fname = datetime.now().strftime("_%Y%m%d")
    time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    uid = datetime.now().strftime("%y%m%d%H%M%S")
    
    if isfile(join(sb.tmppath,"tmp.bias.fits")):
        if 'yes' in input("tmp.bias.fits found,\nDo you want to upload the file? [yes/no] "):
            data = fits.getdata(join(sb.tmppath,"tmp.bias.fits"))
        else:
            data = sb.extractbias(darks)
    else:
        data = sb.extractbias(darks)
    
    hdu = fits.PrimaryHDU(data = data)
    hdu.header.add_comment("Built using files:")
    for fname in darks:
        hdu.header.add_comment(fname)
    hdu.writeto(join(sb.tmppath,"tmp.bias.fits"),overwrite=True)
    
    bias_level = [side(im) for im in data]
    for i in range(len(data)):
        data[i]-=bias_level[i]
    
    #side ref. pixels bias level
    f,ax = plt.subplots(figsize=(14,14))
    name_x = [fname for fname in darks]
    ax.set_xticklabels(name_x,rotation=-45)
    ax.plot(name_x,bias_level,'d',markersize=3)
    ax.set(xlabel='File name',ylabel='Side ref. pixels bias level (ADU)')
    plt.tight_layout()
    
    f.savefig(join(sb.path,'bias%s.png'%fname))
    print('bias%s.png saved.'%fname)
    #show cross section of all data
    f,ax = plt.subplots(figsize=(14,14))
    y2 = [np.median(im,axis=0) for im in data]
    #ax.set_xticklabels(name_x,rotation=-45)
    for i in range(len(name_x)):
        ax.plot(y2[i],label=name_x[i])
        
    ax.set(xlabel='Pixels',ylabel='Bias cross section')
    ax.legend()
    plt.tight_layout()
    print('cross%s.png saved.'%fname)
    f.savefig(join(sb.path,'cross%s.png'%fname))

    superbias = np.median(data,axis=0)+np.median(bias_level)

    
    hdu = fits.PrimaryHDU(data = superbias)#np.uint16(superbias))
    names = [f for f in darks]
    for name in names:
        hdu.header.add_comment(name)
    hdu.header['DATE'] = (time,"Creation date") 
    hdu.header['UNIQUEID'] = (uid, "Unique id of this file")
    if isfile(join(sb.path,"superbias%s.fits"%fname)):
        if 'yes' in input('%s exist, do you want to overwrite it?[yes/no]'%("superbias%s.fits"%fname)):
            hdu.writeto(join(sb.path,"superbias.fits"),overwrite=True)
    else:
        hdu.writeto(join(sb.path,"superbias.fits"))
    
    
