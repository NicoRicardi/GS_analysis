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

cqlist = ["E_FDET_HF", "E_linFDET_HF", "E_ref_HF", "E_ref_HF_CP", "E_ref_MP",
          "E_ref_MP_CP"]
reqs = {
      "E_ref_HF" : ["AB_MP2,scf_energy", "A_MP2,scf_energy", "B_MP2,scf_energy"], 
      "E_ref_HF_CP" : ["AB_MP2,scf_energy", "A_MP2_gh,scf_energy", "B_MP2_gh,scf_energy"], 
      "E_FDET_HF" : ["HF_A,J_int", "HF_A,V_AB", "HF_A,AnucB", "HF_A,BnucA",
                     "HF_A,Exc_nad", "HF_A,Ts_nad",
                     "HF_A,scf_energy", "HF_A,fde_expansion", "HF_B,cycle_energies"],
      "E_linFDET_HF" : ["HF_A,J_int", "HF_A,V_AB", "HF_A,AnucB", "HF_A,BnucA",
                     "HF_A,Exc_nad", "HF_A,Ts_nad", "HF_A,fde_delta_lin",
                     "HF_A,scf_energy", "HF_A,fde_expansion", "HF_B,cycle_energies"]
      }
reqs["E_ref_MP"] = reqs["E_ref_HF"] + ["AB_MP2,mp_correction", "A_MP2,mp_correction", "B_MP2,mp_correction"]
reqs["E_ref_MP_CP"] = reqs["E_ref_HF"] + ["AB_MP2,mp_correction", "A_MP2_gh,mp_correction", "B_MP2_gh,mp_correction"]

quantities_hf = ["DMfinder.json,{}".format(i) for  i in ["HF_FDET_A", "HF_FDET_B"]]

quantities_ccp = list(set(reqs["E_ref_HF"] + reqs["E_ref_HF_CP"] + reqs["E_FDET_HF"]\
                          + reqs["E_ref_MP"] + reqs["E_ref_MP_CP"]))

parser = {k: "dmf_hf" for k in quantities_hf}
parser.update({k: "ccp" for k in quantities_ccp})
parserfuncs = {"dmf_hf": phf, "ccp": None}
parser_kwargs = {"dmf_hf": {}, "ccp": {}}

cwd = os.getcwd()
#systems = ["7HQ_2MeOH", "7HQ_formate", "Uracil_5H2O", "XVI_2HCOOH"]
systems = ["7HQ_2MeOH"]
joblist = []
for system in systems:
    os.chdir(os.path.join(cwd, system))
    joblist.extend([[cwd, system, i] for i in gl.glob("FT*-MC-ME")])
    joblist.append([cwd, system, "FT-ME"])
    joblist.extend([[cwd, system, i] for i in gl.glob("FT*-MC-SE")])
    joblist.append([cwd, system, "FT-SE"])
#print("joblist", joblist)
os.chdir(cwd)
df = ccd.collect_data(joblist, levelnames=["base","system","calc"], qlist=cqlist,
             reqs=reqs, ext="*.out", ignore="slurm*", parser=parser, parserfuncs=parserfuncs,
             parser_file="CCParser.json",parser_args=None, 
             parser_kwargs=parser_kwargs, check_input=True, funcdict = "ccp",
             to_console=True, to_log=False, printlevel=10)   
df.to_csv("FTMC_en.csv")

joblist = []
for system in systems:
    os.chdir(os.path.join(cwd, system))
    joblist.extend([[cwd, system, i] for i in gl.glob("FT?-ME")])
    joblist.extend([[cwd, system, i] for i in gl.glob("FT??-ME")])
    joblist.append([cwd, system, "FT-ME"])
    joblist.extend([[cwd, system, i] for i in gl.glob("FT?-SE")])
    joblist.extend([[cwd, system, i] for i in gl.glob("FT??-SE")])
    joblist.append([cwd, system, "FT-SE"])
#print("joblist", joblist)
os.chdir(cwd)
df = ccd.collect_data(joblist, levelnames=["base","system","calc"], qlist=cqlist,
             reqs=reqs, ext="*.out", ignore="slurm*", parser=parser, parserfuncs=parserfuncs,
             parser_file="CCParser.json",parser_args=None, 
             parser_kwargs=parser_kwargs, check_input=True, funcdict = "ccp",
             to_console=True, to_log=False, printlevel=10)   
df.to_csv("FT_en.csv")


