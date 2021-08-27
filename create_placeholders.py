#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 17:30:03 2021

@author: nico
"""
import os
import itertools as ittl

cwd = os.getcwd()
systems = ["7HQ_2MeOH", "7HQ_formate", "Uracil_5H2O", "XVI_2HCOOH"]
calcs = ["FT-ME", "MC-nopp", "MC-pp_Mulliken", "MC-pp_ChelPG"]
joblist = [i for i in ittl.product([cwd], systems, calcs)]

for job  in joblist:
    path = os.path.join(*job)
    with open(os.path.join(path,"placeholder.out"), "w") as f:
        f.write("")
    