#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 16:36:12 2021

@author: nico
"""
import CCDatabase.CCDatabase as ccd
import pandas as pd
import os
import itertools as ittl

cqlist = ["E_fdet_mp", "E_fdet_hf", "E_ref_HF", "E_ref_HF_CP", "E_ref_MP",
          "E_ref_MP_CP", "elst_change_ref", "elst_change_fdet"]
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

joblist = [("WORK","database","folder"), ("WORK","database","folder2")]
cwd = os.getcwd()
#systems = ["7HQ_2MeOH", "7HQ_formate", "Uracil_5H2O", "XVI_2HCOOH"]
systems = ["Uracil_5H2O"]
#calcs = ["FT-ME", "MC-nopp", "MC-pp_Mulliken", "MC_pp_ChelPG"]7
calcs = ["FT-ME"]
joblist = [i for i in ittl.product([cwd], systems, calcs)]
df = ccd.collect_data(joblist, levelnames=["base","system","calc"], qlist=cqlist,
                 reqs=reqs, ext="*.out", ignore="slurm*", parser=None,
                 parser_file="CCParser.json",parser_args=None, 
                 parser_kwargs=None, check_input=True, funcdict = "ccp",
                 to_console=True, to_log=False, printlevel=10)   

df.to_csv("results.csv")