#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 12:05:51 2021

@author: nico
"""
import pandas as pd

df = pd.read_csv("results.csv", index_col=0)
n = df.copy()
convert = lambda x: x*627.503
for column in df.columns:
    try:
        n[column] = df[column].apply(convert)
    except: 
        continue
n.to_csv("results_kcal.csv")