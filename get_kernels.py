#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 19:26:23 2021

@author: nico
"""
import CCDatabase.CCDatabase as ccd
from CCDatabase.DMfinder import get_all as phf
from CCDatabase.DMfinder import find_fdet_dmfiles as pmp
import pandas as pd
import os
import itertools as ittl

cqlist = ["kernel_tot"]
reqs = {}
quantities_hf = ["DMfinder.json,{}".format(i) for  i in ["HF_FDET_A", "HF_FDET_B"]]
quantities_mp = ["DMfinder.json,{}".format(i) for  i in ["MP_FDET_A", "MP_FDET_B"]]
quantities_ccp = ["MP2_A,{}".format(i) for i in ["fde_Tfunc", "fde_Xfunc", "fde_Cfunc"]]
reqs["kernel_tot"] = quantities_hf + quantities_mp + quantities_ccp
parser = {k: "dmf_hf" for k in quantities_hf}
parser.update({k: "dmf_mp" for k in quantities_mp})
parser.update({k: "ccp" for k in quantities_ccp})
parserfuncs = {"dmf_hf": phf, "dmf_mp": pmp, "ccp": None}
parser_kwargs = {"dmf_hf": {}, "dmf_mp": {"filename": "Densmat_MP.txt", "prop_key": "MP_FDET"}, "ccp": {}}
      
cwd = os.getcwd()
systems = ["7HQ_2MeOH", "7HQ_formate", "Uracil_5H2O", "XVI_2HCOOH"]
calcs = ["MC-nopp", "MC-pp_Mulliken", "MC-pp_ChelPG", "FT-ME"]
joblist = [i for i in ittl.product([cwd], systems, calcs)]
df = ccd.collect_data(joblist, levelnames=["base","system","calc"], qlist=cqlist,
                 reqs=reqs, ext="*.out", ignore="slurm*", parser=parser, parserfuncs=parserfuncs,
                 parser_file="CCParser.json",parser_args=None, 
                 parser_kwargs=parser_kwargs, check_input=True, funcdict = "ccp",
                 to_console=True, to_log=False, printlevel=10)   
df.to_csv("results_kernels.csv")
