#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 12:05:51 2021

@author: nico
"""
import pandas as pd
import sys

if len(sys.argv) == 1:
    orig = "elst_results.csv"
else:
    orig = sys.argv[1]
    
df = pd.read_csv(orig, index_col=0)
n = df.copy()
convert = lambda x: x*627.503
for column in df.columns:
    try:
        n[column] = df[column].apply(convert)
    except: 
        continue
n.to_csv("{}_kcal.csv".format(orig[:-4]))