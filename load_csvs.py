#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 08:31:39 2021

@author: nico
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
plt.rc('font', size=22)
plt.rc('axes', titlesize=22)
plt.rc('axes', labelsize=22)
plt.rc('xtick', labelsize=22)
plt.rc('ytick', labelsize=22)
plt.rc('figure', titlesize=22)
plt.rc('legend', fontsize=22)

enD = pd.read_csv("cc-pVDZ/results_en_kcal.csv", index_col=0)
densD = pd.read_csv("cc-pVDZ/results_densities.csv", index_col=0)

enDa = pd.read_csv("results_en_kcal.csv", index_col=0)
densDa = pd.read_csv("results_densities.csv", index_col=0)

enT = pd.read_csv("cc-pVTZ/results_en_kcal.csv", index_col=0)
densT = pd.read_csv("cc-pVTZ/results_densities.csv", index_col=0)

enD_SE = pd.read_csv("cc-pVDZ/results_SE_en_kcal.csv", index_col=0)
densD_SE = pd.read_csv("cc-pVDZ/results_SE_densities.csv", index_col=0)

enDa_SE = pd.read_csv("results_SE_en_kcal.csv", index_col=0)
densDa_SE = pd.read_csv("results_SE_densities.csv", index_col=0)

enT_SE = pd.read_csv("cc-pVTZ/results_SE_en_kcal.csv", index_col=0)
densT_SE = pd.read_csv("cc-pVTZ/results_SE_densities.csv", index_col=0)

names = ["cc-pVDZ", "aug-cc-pVDZ", "cc-pVTZ", "cc-pVDZ-SE", "aug-cc-pVDZ-SE", "cc-pVTZ-SE"]
for n,en in enumerate([enD, enDa, enT, enD_SE, enDa_SE, enT_SE]):
    en["E_MPk"] = en["E_FDET_MP"] - en["kernel_tot"]
    en["BSSE"] = en["E_ref_MP_CP"] - en["E_ref_MP"]
    en.to_csv("results_proc_{}.csv".format(names[n]))
