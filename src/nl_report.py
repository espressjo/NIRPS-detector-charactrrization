#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 10:19:02 2022

@author: noboru
"""
from astropy.io import fits
from nl_fit_artigau import nl
n = nl()
#with bias subtracted
f_meas = fits.getdata('/nirps_raw/nirps/characterization/20220505163017.bias.fits')[:,2000,2000]
n.fit_graph_rapport(f_meas,xpos=2000,ypos=2000)
#new nonlin.fits
n.make_report(nonlinf='nonlin.fits',cpath='/nirps_raw/nirps/characterization',sbias='superbias220517095751.fits')

#old nonlin
n = nl()
n.make_report(nonlinf='nonlin.fits',cpath='/spirou_slow/NIRPS/characterization',sbias='superbias_20210416.fits')
