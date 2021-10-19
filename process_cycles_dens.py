#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 19:26:23 2021

@author: nico
"""
import CCDatabase.CCDatabase as ccd
from CCDatabase.DMfinder import get_all as pfunc
import pandas as pd
import glob as gl
import os
#import itertools as ittl

cqlist = ["densdiff_FDET_ref", "densdiff_iso_ref", "densdiff_iso_FDET",  "M_value"]
reqs = {}
reqs["densdiff_FDET_ref"] = ["DMfinder.json,{}".format(i) for  i in ["HF_ref", "HF_FDET_A", "HF_FDET_B"]]
reqs["densdiff_iso_ref"] = ["DMfinder.json,{}".format(i) for  i in ["HF_iso_A", "HF_iso_B", "HF_ref"]]
reqs["densdiff_iso_FDET"] = ["DMfinder.json,{}".format(i) for  i in ["HF_iso_A", "HF_iso_B", "HF_FDET_A", "HF_FDET_B"]]
reqs["M_value"] = ["DMfinder.json,{}".format(i) for  i in ["HF_iso_B", "HF_ref"]]
parser = "dmf"
parserfuncs = pfunc
      
cwd = os.getcwd()
systems = ["7HQ_2MeOH", "7HQ_formate", "Uracil_5H2O", "XVI_2HCOOH"]
calcs = ["FT-ME", "FT-SE"]
for calc in calcs:
    joblist = []
    for system in systems:
        base = os.path.join(cwd, system, calc)
        cyfols = [i.replace(os.path.join(base, ""), "") for i in gl.glob(os.path.join(base, "cy*", "FT*-MC-*"))]
        cyfols = [i for i in cyfols if not int([j for j in i if i.isnumeric()][-1]) % 2]  # only even
        print("cyfols", cyfols)
        joblist.extend([[cwd, system, calc, i] for i in cyfols])
    print("joblist", joblist)
#    df = ccd.collect_data(joblist, levelnames=["base","system","calc", "cy"], qlist=cqlist,
#                 reqs=reqs, ext="*.out", ignore="slurm*", parser=parser, parserfuncs=parserfuncs,
#                 parser_file="CCParser.json",parser_args=None, 
#                 parser_kwargs=None, check_input=True, funcdict = "ccp",
#                 to_console=True, to_log=False, printlevel=10)   
#    df.to_csv("FTMC_{}_dens.csv".format(calc[-2:]))
