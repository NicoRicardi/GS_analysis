#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 13:58:16 2021

@author: nico
"""
import os
import glob as gl
import itertools as ittl
from CCDatabase.CCDatabase import correct_values

cwd = os.getcwd()
systems = ["7HQ_2MeOH", "7HQ_formate", "Uracil_5H2O", "XVI_2HCOOH"]
#dmfinders =[os.path.join(*i, "DMfinder.json") for i in ittl.product([cwd], systems, calcs)]

systems = ["7HQ_2MeOH", "7HQ_formate", "Uracil_5H2O", "XVI_2HCOOH"]
joblist = []
for system in systems:
    os.chdir(os.path.join(cwd, system))
    joblist.extend([[cwd, system, i] for i in gl.glob("FT*-MC-ME")])
    joblist.append([cwd, system, "FT-ME"])
    joblist.extend([[cwd, system, i] for i in gl.glob("FT*-MC-SE")])
    joblist.append([cwd, system, "FT-SE"])
os.chdir(cwd)
datafiles = [os.path.join(*i, "data.json") for i in joblist]
correct_values(datafiles, todel="E_FDET_HF")

joblist = []
for system in systems:
    os.chdir(os.path.join(cwd, system))
    joblist.extend([[cwd, system, i] for i in gl.glob("FT*-ME")])
    joblist.append([cwd, system, "FT-ME"])
    joblist.extend([[cwd, system, i] for i in gl.glob("FT*-SE")])
    joblist.append([cwd, system, "FT-SE"])
os.chdir(cwd)
datafiles = [os.path.join(*i, "data.json") for i in joblist]
correct_values(datafiles, todel="E_FDET_HF")

#correct_values(dmfinders, todel=["HF_ref", "HF_FDET_A"])