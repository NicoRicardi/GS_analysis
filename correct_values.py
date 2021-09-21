#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 13:58:16 2021

@author: nico
"""
import os
import itertools as ittl
from CCDatabase.CCDatabase import correct_values

cwd = os.getcwd()
systems = ["7HQ_2MeOH", "7HQ_formate", "Uracil_5H2O", "XVI_2HCOOH"]
calcs = ["FT-ME", "MC-nopp", "MC-pp_Mulliken", "MC-pp_ChelPG"]
datafiles = [os.path.join(*i, "data.json") for i in ittl.product([cwd], systems, calcs)]
#dmfinders =[os.path.join(*i, "DMfinder.json") for i in ittl.product([cwd], systems, calcs)]


correct_values(datafiles, todel="kernel_tot")
#correct_values(dmfinders, todel=["HF_ref", "HF_FDET_A"])