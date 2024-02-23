#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 14:09:53 2022

@author: noboru
"""

from matplotlib import pyplot as plt
import seaborn as sns
from astropy.io import fits
from astropy.stats import sigma_clipped_stats as sc

sns.set_theme()
import numpy as np

old_nl = "/spirou_slow/NIRPS/characterization/nonlin.fits"
new_nl = "/nirps_raw/nirps/characterization/nonlin.fits"

hdu_old = fits.open(old_nl)
a_old = hdu_old[0].data
a_new = fits.getdata(new_nl)[0].ravel()

b_old = hdu_old[1].data
b_new = fits.getdata(new_nl)[1].ravel()


mn,_,std = sc(a_old)
a_old = a_old[(a_old<(mn+5*std)) & (a_old>(mn-5*std))]
mn,_,std = sc(a_new)
a_new = a_new[(a_new<(mn+5*std)) & (a_new>(mn-5*std))]

mn,_,std = sc(b_old)
b_old = b_old[(b_old<(mn+5*std)) & (b_old>(mn-5*std))]
mn,_,std = sc(b_new)
b_new = b_new[(b_new<(mn+5*std)) & (b_new>(mn-5*std))]
f,ax = plt.subplots(2,1)
plt.suptitle("Nonlinearity Coefficient C$_3$ and C$_2$")

ax[0].hist(a_new,bins=200,label="Chile")
ax[0].hist(a_old.ravel(),bins=200,label="COPL",alpha=0.3)
ax[0].legend() #.legend()
ax[0].set(title="C$_3$",ylabel="Counts")



ax[1].hist(b_new,bins=200,label="Chile")
ax[1].hist(b_old.ravel(),bins=200,label="COPL",alpha=0.3)
ax[1].legend() #.legend()
ax[1].set(title="C$_2$",ylabel="Counts",xlabel="bins")
plt.tight_layout()
f.savefig("/nirps_raw/nirps/characterization/nonlin_coeff_comp.png")
plt.show()