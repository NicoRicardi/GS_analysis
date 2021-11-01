#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: nico
"""
import os
import itertools as ittl
import glob as gl

cwd = os.getcwd()
systems = ["7HQ_2MeOH", "7HQ_formate", "Uracil_5H2O", "XVI_2HCOOH"]
#calcs = ["FT-ME", "MC-nopp", "MC-pp_Mulliken", "MC-pp_ChelPG"]

joblist = []
for system in systems:
#    joblist.append([cwd, system, "FT-SE"])
    os.chdir(os.path.join(cwd, system))
    joblist.extend([[cwd, system, i] for i in gl.glob("FT*-MC-?E")])
    joblist.extend([[cwd, system, i] for i in gl.glob("FT?-?E")])
    joblist.extend([[cwd, system, i] for i in gl.glob("FT??-?E")])
os.chdir(cwd)
#joblist = [i for i in ittl.product([cwd], systems, calcs)]

for job  in joblist:
    path = os.path.join(*job)
#    print(path)
    with open(os.path.join(path,"placeholder.out"), "w") as f:
        f.write("")
