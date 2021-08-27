#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 16:36:12 2021

@author: nico
"""
import CCDatabase.CCDatabase as ccd
from CCDatabase.get_elst_int_sum import elst_int_sum_iso as pfunc
import pandas as pd
import os
import itertools as ittl


cqlist = ["elst_change_ref", "elst_change_FDET"]
reqs = {}
reqs_sum_iso_pfunc = ["FT-ME/cy0,J_sum_iso", "FT-ME/cy0,AnucB_sum_iso","FT-ME/cy0,BnucA_sum_iso"]
reqs_sum_iso_ccp = ["A_MP2,tot_coulomb", "A_MP2,nuc_attr", "A_MP2,nuc_repu",
                    "B_MP2,tot_coulomb", "B_MP2,nuc_attr", "B_MP2,nuc_repu",
                    "FT-ME/cy0,V_AB"]
reqs_elst_AB = ["AB_MP2,tot_coulomb", "AB_MP2,nuc_attr", "AB_MP2,nuc_repu"]
reqs_elst_FDET = ["MP2_A,tot_coulomb", "MP2_A,nuc_attr", "MP2_A,nuc_repu",
    "MP2_A,J_int", "MP2_A,V_AB", "MP2_A,AnucB", "MP2_A,BnucA",
    "MP2_B,tot_coulomb", "MP2_B,nuc_attr", "MP2_B,nuc_repu"]
reqs["elst_change_ref"] = reqs_sum_iso_pfunc + reqs_sum_iso_ccp + reqs_elst_AB
reqs["elst_change_FDET"] = reqs_sum_iso_pfunc + reqs_sum_iso_ccp + reqs_elst_FDET
parser = {i: "ccp" for i in reqs_elst_AB + reqs_elst_FDET + reqs_sum_iso_ccp}
parser.update(**{i: "elst" for i in reqs_sum_iso_pfunc})
parserfuncs = {"ccp": None, "elst": pfunc}
      
cwd = os.getcwd()
systems = ["7HQ_2MeOH", "7HQ_formate", "Uracil_5H2O", "XVI_2HCOOH"]
calcs = ["FT-ME", "MC-nopp", "MC-pp_Mulliken", "MC-pp_ChelPG"]
joblist = [i for i in ittl.product([cwd], systems, calcs)]
df = ccd.collect_data(joblist, levelnames=["base","system","calc"], qlist=cqlist,
                 reqs=reqs, ext="*.out", ignore="slurm*", parser=parser, parserfuncs=parserfuncs,
                 parser_file="CCParser.json",parser_args=None, 
                 parser_kwargs=None, check_input=True, funcdict = "ccp",
                 to_console=True, to_log=False, printlevel=10)   
df.to_csv("results_elst.csv")

