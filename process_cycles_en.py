#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 16:36:12 2021

@author: nico
"""
import CCDatabase.CCDatabase as ccd
from CCDatabase.DMfinder import get_all as phf
from CCDatabase.DMfinder import find_fdet_dmfiles as pmp
import glob as gl
import pandas as pd
import os
import itertools as ittl

cqlist = ["E_FDET_MP", "kernel_tot", "E_FDET_HF", "E_ref_HF", "E_ref_HF_CP", "E_ref_MP",
          "E_ref_MP_CP"]
reqs = {
      "E_ref_HF" : ["AB_MP2,scf_energy", "A_MP2,scf_energy", "B_MP2,scf_energy"], 
      "E_ref_HF_CP" : ["AB_MP2,scf_energy", "A_MP2_gh,scf_energy", "B_MP2_gh,scf_energy"], 
      "E_FDET_HF" : ["MP2_A,J_int", "MP2_A,V_AB", "MP2_A,AnucB", "MP2_A,BnucA",
                     "MP2_A,Exc_nad", "MP2_A,Ts_nad", "MP2_A,fde_delta_lin",
                     "MP2_A,scf_energy", "MP2_A,fde_expansion", "MP2_B,cycle_energies"]
      }
reqs["E_ref_MP"] = reqs["E_ref_HF"] + ["AB_MP2,mp_correction", "A_MP2,mp_correction", "B_MP2,mp_correction"]
reqs["E_ref_MP_CP"] = reqs["E_ref_HF"] + ["AB_MP2,mp_correction", "A_MP2_gh,mp_correction", "B_MP2_gh,mp_correction"]
reqs["E_FDET_MP"] = reqs["E_FDET_HF"] + ["MP2_A,mp_correction", "MP2_B,mp_correction"]

quantities_hf = ["DMfinder.json,{}".format(i) for  i in ["HF_FDET_A", "HF_FDET_B"]]
quantities_mp = ["DMfinder.json,{}".format(i) for  i in ["MP_FDET_A", "MP_FDET_B"]]
quantities_ccp_kernel = ["MP2_A,{}".format(i) for i in ["fde_Tfunc", "fde_Xfunc", "fde_Cfunc"]]
reqs["kernel_tot"] = quantities_hf + quantities_mp + quantities_ccp_kernel

quantities_ccp = list(set(reqs["E_ref_HF"] + reqs["E_ref_HF_CP"] + reqs["E_FDET_HF"]\
                          + reqs["E_ref_MP"] + reqs["E_ref_MP_CP"] + reqs["E_FDET_MP"]))
parser = {k: "dmf_hf" for k in quantities_hf}
parser.update({k: "dmf_mp" for k in quantities_mp})
parser.update({k: "ccp" for k in quantities_ccp_kernel + quantities_ccp})
parserfuncs = {"dmf_hf": phf, "dmf_mp": pmp, "ccp": None}
parser_kwargs = {"dmf_hf": {}, "dmf_mp": {"filename": "Densmat_MP.txt", "prop_key": "MP_FDET"}, "ccp": {}}

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
#                 parser_kwargs=parser_kwargs, check_input=True, funcdict = "ccp",
#                 to_console=True, to_log=False, printlevel=10)   
#    df.to_csv("FTMC_{}_en.csv".format(calc[-2:]))


