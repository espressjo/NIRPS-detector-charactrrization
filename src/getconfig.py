#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 21:16:37 2024

@author: espressjo
"""


P = "./config"

def getString(key,default=""):
    value = default
    with open(P,'r') as f:
        for l in f.readlines():
            if key in l:
                value = l.replace(key,"").strip()
    return value