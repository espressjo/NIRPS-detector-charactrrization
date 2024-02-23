#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find NIRPS's ramp file based on unique identification numer.


Created on Thu Apr  8 10:15:54 2021

@author: jonathan st-antoine
"""
from os.path import isdir,join
from datetime import datetime,timedelta
from os import listdir
from astropy.io import fits
from getconfig import getString





class findramp():
    def __init__(self):
        self.datapath =  getString("DATAPATH",default='/nirps_raw/read') 
    def _kinit(self,uid):
        if isinstance(uid, int):
            uid = str(uid)
        y = int(uid[:4])
        m = int(uid[4:6])
        d = int(uid[6:8])
        H = int(uid[8:10])
        M = int(uid[10:12])
        S = int(uid[12:14])
        dt = datetime(y,m,d,H,M,S)
        delta = timedelta(days=1)
        dt1 = dt-delta
        dt2 = dt+delta
        ls = [];
        if isdir(join(self.datapath,dt1.strftime("r%Y%m%d"))):
            ls.append(join(self.datapath,dt1.strftime("r%Y%m%d")))
        if isdir(join(self.datapath,dt.strftime("r%Y%m%d"))):
            ls.append(join(self.datapath,dt.strftime("r%Y%m%d")))
        if isdir(join(self.datapath,dt2.strftime("r%Y%m%d"))):
            ls.append(join(self.datapath,dt2.strftime("r%Y%m%d")))
        return ls
        
    def findramp(self,uid):
        ls = self._kinit(uid)
        for d in ls:
            ls = [join(d,f) for f in listdir(d) if 'fits' in f]
            for fname in ls:
                try:
                    if uid==fits.getheader(fname)["UNIQUEID"]:
                        return fname
                except :
                    pass
if '__main__' in __name__:
    uid = 20210403171440
    fr = findramp()    
    print(fr.findramp(str(uid)))
