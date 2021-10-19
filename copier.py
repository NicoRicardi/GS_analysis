#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 21:52:17 2021

@author: nico
"""
import os
import subprocess as sp
import shutil as sh

sourcefol = "/home/nico/baobab/NEG_GS/XVI_2HCOOH/FT-SE"
destfol = "/home/nico/yggdrasil/NEG_GS/XVI_2HCOOH_FT-SE"
basename = "cy"

cyfols = [i for i in os.listdir(sourcefol) if i.startswith(basename)]
tocopy = ["Densmat_B.txt", "FDE_State0_tot_dens.txt"]
if not os.path.isdir(destfol):
    os.makedirs(destfol)
for cyfol in cyfols:
    dcyfol = os.path.join(destfol, os.path.basename(cyfol))
    if not os.path.isdir(dcyfol):
        os.makedirs(dcyfol)
    tmp = [os.path.join(sourcefol, cyfol, i) for i in os.listdir(os.path.join(sourcefol, cyfol)) if os.path.isdir(os.path.join(sourcefol, cyfol, i))]
    if len(tmp) == 0:
        continue
    ftmcfol = tmp[0]
    dftmcfol = os.path.join(dcyfol, os.path.basename(ftmcfol))
    if not os.path.isdir(dftmcfol):
        os.makedirs(dftmcfol)
    for fname in ["meta.json", "energies.txt"]:
        try:
            sh.copy(os.path.join(ftmcfol, fname), os.path.join(dftmcfol, fname))
        except:
            print("some issue in {}".format(ftmcfol))
    last_it = max([int(i[2:]) for i in os.listdir(ftmcfol) if i.startswith(basename)])
    lastfol = os.path.join(ftmcfol,"{}{}".format(basename, last_it))
    dlastfol = os.path.join(dftmcfol, os.path.basename(lastfol))
    if not os.path.isdir(dlastfol):
        os.makedirs(dlastfol)
    for tc in tocopy:
        sh.copy(os.path.join(lastfol, tc), os.path.join(dlastfol, tc))
    
    


