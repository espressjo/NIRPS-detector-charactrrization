#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 16:08:29 2021

@author: noboru
"""

#capacitive coupling
from hxrg.type.ramp import hxramp
from hxrg.type.read import hxread

from numpy import zeros,argwhere
from astropy.io import fits
from os.path import join
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
cpath = '/nirps_raw/nirps/characterization'
p_data = "/nirps_raw/nirps/reads"
class ccoupling():
    def __init__(self):
        self.data = zeros((4096,4096))
        self.candidate = zeros((4096,4096))
    def get_data(self,fname):
        if '/' not in fname:
            fname = join(p_data,fname)
        R = hxramp(toponly = True)
        r = hxread()
        H = fits.getheader(fname)
        cube = fits.getdata(fname)
        print("Creating ramp file")
        for im in tqdm(cube):
            r(np.asarray(im,dtype=float))
            r.H = H
            R<<r
        R.fit()
        R.unfold()
        self.data = R.a
        return 
    def cap_coupling(self):
        
        #mask = zeros((4096,4096),dtype=bool)
        #mask[(self.data >10) & (self.data<12)]=True
        w = argwhere((self.data >10) & (self.data<12))
        w = [c for c in w if all([c[0]>200,c[0]<3800,c[1]>200,c[1]<3800])]
        for c in w:
            print(c[0],c[1])
        
        cube =([np.asarray(self.data[coord[0]-10:coord[0]+11,coord[1]-10:coord[1]+11]) for coord in w ])
        cube = [im for im in cube if all([len(im[im>2])<2,len(im[im>13])==0])]
        print("There is %d hot pixels in this ramp."%len(cube))
        cube = [ (im/np.nanmedian(im))-1 for im in cube]
        cube = [im/np.sum(im) for im in cube]
        return np.nanmedian(cube,axis=0)        
        #
        
        #hdu = fits.PrimaryHDU(data=np.asarray(cube))
        #hdu.writeto("/home/jonathan/c11.fits",overwrite=True)
        
        #hdu = fits.PrimaryHDU(data=np.asarray(np.nanmedian(cube,axis=0)))
        #hdu.writeto("/home/jonathan/c12.fits",overwrite=True)
    def plot(self,array,color='y',alpha=0.9):
        """
        **Description**:
            This function generate a 3D plot from a 5x5 array. The output of stackHotPixel should be used as input of plotIPC.
        **color**:
            color of the plot
        **alpha**:
            transparancy of the bar plot
            """
        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set_theme()
        A = np.array(array)*100
        A[A<.1]=0
        bkp = A[2,2]
        A[2,2]=6
        fig = plt.figure(figsize=(12,12))
        ax = fig.add_subplot(111,projection='3d')
        x,y = np.meshgrid(range(5),range(5))
        z = np.zeros(len(A.ravel()))
        dx = np.ones(len(x.ravel()))-0.1
        dy = np.ones(len(y.ravel()))-0.1
        ax.bar3d(x.ravel(),y.ravel(),z,dx,dy,A.ravel(),alpha=0.7)
        ax.set_zlim([0,8])
        for i in range(5):
            for ii in range(5):
                if np.any([ii!=2,i!=2]):
                    if A[ii,i]<0.001:
                        ax.text(i+0.4,ii+0.3,A[ii,i],"%d%%"%float(A[ii,i]),c='w',fontsize=8,alpha=0.9)
                    else:
                        ax.text(i+0.3,ii+0.3,A[ii,i],"%2.1f%%"%float(A[ii,i]),c='w',fontsize=8,alpha=0.9)
                else:
                    ax.text(i+0.35,ii+0.35,A[ii,i],"%d%%"%float(bkp),c='w',fontsize=8)
        ax.set_xlim([0,5])
        ax.set_ylim([0,5])
        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.grid(False)
        ax.zaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticklabels([])
        ax.xaxis.set_ticks([],[])
        ax.yaxis.set_ticks([],[])
        ax.zaxis.set_ticks([],[])
        for line in ax.zaxis.get_ticklines():
            line.set_visible(False)
        plt.tight_layout()
        plt.show()
        return fig
        #ax.get_zaxis().set_visible(False)
        
        #ax.set_axis_off()       
        
        
if '__main__' in __name__:
    darks = ["NIRPS_2022-05-05T15_31_01_590.fits",
            "NIRPS_2022-05-05T15_40_24_489.fits",
            "NIRPS_2022-05-05T15_49_47_388.fits",
            "NIRPS_2022-05-05T15_59_10_285.fits"]
    cmaps = [];
    from tqdm import tqdm
    for fname in darks:
        
        coupling = ccoupling()
        coupling.get_data(fname)
        cmaps.append(coupling.cap_coupling())
    capmap = np.nanmedian(cmaps,axis=0)
    hdu = fits.PrimaryHDU(data=capmap)
    hdu.header.add_comment("built using: ")
    from datetime import datetime
    fmt = "%Y-%m-%dT%H:%M:%S"
    fmtuid = "%Y%m%d%H%M%S"
    hdu.header['DATE'] = datetime.now().strftime(fmt)
    hdu.header['UID'] = datetime.now().strftime(fmtuid)
    
    for fname in darks:
        hdu.header.add_comment(fname)
    hdu.writeto(join(cpath,'capacitive_coupling_model.fits'),overwrite=True)
        
    _c = ((1.0-capmap[10,10]) +1.)**2 -1.
    print("Gain will be off by %.2f %%"%(_c*100))
    fig = coupling.plot(capmap[8:13,8:13])
    fig.savefig(join(cpath,'capacitive_coupling.png'))
