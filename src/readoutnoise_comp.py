#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 13:25:48 2022

@author: noboru
"""
from astropy.stats import sigma_clipped_stats as sc
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_theme()
fig,ax = plt.subplots()
gain_old = 1.27
gain_new = 1.27

files_old = ["/spirou_slow/NIRPS/characterization/.tmp/20210407053206.readout.to.npy",
             "/spirou_slow/NIRPS/characterization/.tmp/20210407041021.readout.to.npy",
             "/spirou_slow/NIRPS/characterization/.tmp/20210407042034.readout.to.npy",
             "/spirou_slow/NIRPS/characterization/.tmp/20210407043047.readout.to.npy",
             "/spirou_slow/NIRPS/characterization/.tmp/20210407044100.readout.to.npy",
             "/spirou_slow/NIRPS/characterization/.tmp/20210407045113.readout.to.npy",
             "/spirou_slow/NIRPS/characterization/.tmp/20210407050127.readout.to.npy"]
files_new = ["/nirps_raw/nirps/characterization/.tmp/NIRPS_2022-05-05T15_49_47_388.fits.readout.to.npy",
             "/nirps_raw/nirps/characterization/.tmp/NIRPS_2022-05-05T15_40_24_489.fits.readout.to.npy",
             "/nirps_raw/nirps/characterization/.tmp/NIRPS_2022-05-05T15_31_01_590.fits.readout.to.npy",
             "/nirps_raw/nirps/characterization/.tmp/NIRPS_2022-05-05T15_59_10_285.fits.readout.to.npy"]
#old data
ro_old = [];
err_old = [];
rel_err_old = 0
for f in files_old:#"/spirou_slow/NIRPS/characterization/.tmp/20210407053206.readout.to.npy"
    arr = np.load(f)
    stats = [sc(np.asarray(r)*gain_old) for r in arr]
    S,_,err = zip(*stats)
    _ro_old,_err_old = np.asarray(S),np.asarray(err)
    ro_old.append(_ro_old)
    err_old.append(_err_old)
    rel_err_old+=(_err_old/_ro_old)**2

ro_old = np.mean(ro_old,axis=0)
rel_err_old = [np.sqrt(t) for t in rel_err_old]
err_old = np.asarray(ro_old)*np.asarray(rel_err_old)
print("[COPL] CDS noise: %.3f+/-%.3f"%(ro_old[0],err_old[0]))
#new data
#arr = np.load("/nirps_raw/nirps/characterization/.tmp/NIRPS_2022-05-05T15_49_47_388.fits.readout.to.npy")
#new data
ro_new = [];
err_new = [];
rel_err_new = 0
for f in files_new:#"/spirou_slow/NIRPS/characterization/.tmp/20210407053206.readout.to.npy"
    arr = np.load(f)
    stats = [sc(np.asarray(r)*gain_new) for r in arr]
    S,_,err = zip(*stats)
    _ro_new,_err_new = np.asarray(S),np.asarray(err)
    ro_new.append(_ro_new)
    err_new.append(_err_new)
    rel_err_new+=(_err_new/_ro_new)**2
    
ro_new = np.mean(ro_new,axis=0)
rel_err_new = [np.sqrt(t) for t in rel_err_new]
err_new = np.asarray(ro_new)*np.asarray(rel_err_new)
print("[Chile] CDS noise: %.3f+/-%.3f"%(ro_new[0],err_new[0]))

fig,ax = plt.subplots()
#decay for old data
x_old = np.arange(2,len(ro_old)+2,1)
a = ro_old[1]/(3**(-1.0/2.0))
y2 = [a*i**(-1.0/2.0) for i in x_old]
ax.plot(x_old,y2,'--',markersize=2,label="N$^{-1/2}$ decay")


#old data
ax.errorbar(x_old,ro_old,yerr=err_old,c='b',fmt='')
ax.plot(x_old,ro_old,'d',c='b',markersize=2,label="NIRPS [COPL]")

#new data
x_new = np.arange(2,len(ro_new)+2,1)
ax.errorbar(x_new,ro_new,yerr=err_new,c='g',fmt='')
ax.plot(x_new,ro_new,'d',c='g',markersize=2,label="NIRPS [Chile]")


ax.set(title='Readout noise',xlabel='Read number',ylabel='Readout noise (electron)')
ax.set_xlim([2,len(ro_new)])
ax.legend()
fig.savefig("/nirps_raw/nirps/characterization/readoutnoise_comp_copl_chile.png")
plt.show()