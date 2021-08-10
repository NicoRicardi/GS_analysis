#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 17:12:11 2021

@author: nico
"""
#import json as js
#import subprocess as sp
#import CCParser as ccp
#import numpy as np
import glob as gl
import traceback
import subprocess as sp
import os
import CCJob as ccj
import CCJob.utils as ut
from CCJob.iterative import freeze_and_thaw, macrocycles, copy_density
import shutil as sh
from CCJob.Composable_templates import Tdefaults, Tinps, Tinp, Trem_kw, Tmolecule,\
 Tadc,  Tpc, Tbasis, chelpg_kw, Tfragments, Tfde
#import slurm_stuff_yggdrasil as slrm
import slurm_stuff_baobab as slrm
import sys
import logging

###############################################################################
###############################################################################
root = os.getcwd()
system = sys.argv[1]
systfol = os.path.join(root, system)
###############################################################################
logfile = os.path.join(systfol,"{}_ccj.log".format(system))
# set up logger. NB avoid homonimity with other module's loggers (e.g. ccp)
ccjlog = logging.getLogger("ccj")
#ut.setupLogger(to_console=True, to_log=True, printlevel=20)
ut.setupLogger(to_console=True, to_log=True, logname=logfile)
bs_kw = "gen"  
bs_string = ut.read_file("aug-cc-pVDZ.bas")
if bs_string[-1] == "\n":
    bs_string = bs_string[:-1]
rem_adc_basic = dict(method="adc(2)", basis=bs_kw, ee_specs=Tadc.substitute(Tdefaults["adc"]))
rem_hf_basic = dict(method="hf", basis=bs_kw)
extra_basic = Tbasis.substitute(**{"basis_specs": bs_string})
    
ut.logchdir(ccjlog,systfol)
# --- First things first : STRUCTURE ---
zr_file = ccj.find_file(systfol, extension="zr")
frags = ccj.zr_frag(zr_file)

# --- Electronic configuration ---
f_elconf = ccj.find_eleconfig(systfol)
elconf = ccj.read_eleconfig(fname=f_elconf)
found_elconf = False if f_elconf == 0 else True
# shortcuts for electronic configuration
if found_elconf:
    elconf_A   = dict(charge=elconf["charge_a"], 
                      multiplicity=elconf["multiplicity_a"])
    elconf_B   = dict(charge=elconf["charge_b"], 
                      multiplicity=elconf["multiplicity_b"])
    elconf_AB  = dict(charge=elconf["charge_tot"], 
                      multiplicity=elconf["multiplicity_tot"])
    elconf_rev = dict(**elconf)
    elconf_rev["charge_a"] = elconf["charge_b"]
    elconf_rev["charge_b"] = elconf["charge_a"]
    elconf_rev["multiplicity_a"] = elconf["multiplicity_b"]
    elconf_rev["multiplicity_b"] = elconf["multiplicity_a"]
        
#-----------------------------------------------------------------------------#
# 1 Freeze-and-Thaw, SE
#-----------------------------------------------------------------------------#
memory = 28000
rem_kw = dict(**rem_hf_basic, **{"memory": memory})
fde_kw = dict(Tdefaults["fde"], **{"method_a": "import_rhoA true", "method_b": "import_rhoB true", "expansion": "SE"})
specs_fnt = dict(rem_kw=rem_kw, fde_kw=fde_kw, extras=extra_basic, use_zr=False,
                fragments=frags, elconf=elconf, q_custom=slrm.slurm_add,
                maxiter=20, thresh=1e-9, en_file="energies.txt")  
queue_fnt = dict(**slrm.shabug_M)  
meta_fnt  = dict(method_A="HF", method_B="HF", opt="freeze-thaw", status=None)

# create calculation folder
meta_fnt["path"] = os.path.join(systfol, "FT-SE")
en_file = os.path.join(meta_fnt["path"], specs_fnt["en_file"])
already_done_fnt = ut.status_ok(path=meta_fnt["path"])

if already_done_fnt == False:
    # run calculation and update status ("checkpoint")
    try:
        ut.mkdif(meta_fnt["path"])
        ut.logchdir(ccjlog,meta_fnt["path"])
        if not os.path.exists(os.path.join(meta_fnt["path"],"Densmat_B.txt")):
            copy_density(os.path.join(systfol, "B_MP2_gh", "Densmat_SCF.txt"),
                         "Densmat_B.txt", header_src=False, alpha_only_src=False) 
        if not os.path.exists(os.path.join(meta_fnt["path"],"Densmat_A.txt")):
            copy_density(os.path.join(systfol, "A_MP2_gh", "Densmat_SCF.txt"),
                    "Densmat_A.txt", header_src=False, alpha_only_src=False)
        freeze_and_thaw(queue_fnt, **specs_fnt)  
        meta_fnt["status"] = "FIN"
        already_done_fnt = True
    except Exception as e:
        print(e)
        traceback.print_exc()
        meta_fnt["status"] = "FAIL"
    # serialize current meta information for later
    ut.save_status(meta_fnt)

#-----------------------------------------------------------------------------#    
# 2: FDE-MP2 using FT densities: (A-in-B) [SE]
#-----------------------------------------------------------------------------#
meta_ftmpa  = dict(method_A="MP2", method_B="import", opt=None, status=None,
              basename="emb")

# get name of calculation folder
meta_ftmpa["path"] = os.path.join(meta_fnt["path"], "MP2_A")
already_done_ftmpa = ut.status_ok(path=meta_ftmpa["path"])

# run calculation and update status ("checkpoint")
if already_done_ftmpa == False and already_done_fnt:
    memory = 87000
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory, "fde": "true"})
    frag_specs = dict(frag_a=frags["A"], frag_b=frags["B"])
    if found_elconf:
        frag_specs.update(elconf)
    frag_str = Tfragments.substitute(Tdefaults["molecule"], **frag_specs)
    fde_sect = Tfde.substitute(Tdefaults["fde"], **{"method_a": "import_rhoA true", "method_b": "import_rhoB true", "expansion": "SE"})
    extras = "\n".join([extra_basic]+[fde_sect])
    specs_ftmpa = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=frag_str, extras=extras)
    queue_ftmpa = dict(**slrm.weso_small1)  
    # calculation hasn't run yet, create folder in order to copy densmat
    ut.mkdif(meta_ftmpa["path"])
    
    # SPECIAL: copy density matrices
    iterDir = ut.get_last_iter_dir(active="A", path=meta_fnt["path"])
    sh.copy(os.path.join(iterDir, "Densmat_B.txt"), meta_ftmpa["path"],)
    copy_density(os.path.join(iterDir, "FDE_State0_tot_dens.txt"),
                 os.path.join(meta_ftmpa["path"],"Densmat_A.txt"),
                 header_src=False, alpha_only_src=False)
    
    # finally run 
    json_files = gl.glob(os.path.join(meta_ftmpa["path"],"*.json"))
    ut.run_job(specs_ftmpa, queue_ftmpa, meta_ftmpa, Tinp, q_custom=slrm.slurm_add,  
            batch_mode=False, create_folder=False)  # because we want to extract data  
    njsf = [i for i in gl.glob(os.path.join(meta_ftmpa["path"],"*.json")) if i not in json_files][0]  # whatever json was just added, i.e. default_ccpjson
    # serialize current meta information for later (we're still in the calc folder)
    ut.save_status(meta_ftmpa)
    data = ut.load_js(njsf)  # default ccp json_filename
    sp.call("echo {E_A} >> {en_file}".format(E_A=data["scf_energy"][-1][0], en_file=en_file), shell=True)

#-----------------------------------------------------------------------------#
# 3: FDE-MP2 using FT densities: (B-in-A) [SE]
#-----------------------------------------------------------------------------#

meta_ftmpb  = dict(method_A="MP2", method_B="import", opt=None, status=None,
              basename="emb")
# get name of calculation folder
meta_ftmpb["path"] = os.path.join(meta_fnt["path"], "MP2_B")
already_done_ftmpb = ut.status_ok(path=meta_ftmpb["path"])

# run calculation and update status ("checkpoint")
if already_done_ftmpb == False and already_done_fnt:
    memory = 87000
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory, "fde": "true"})
    frag_specs = dict(frag_a=frags["B"], frag_b=frags["A"])
    if found_elconf:
        frag_specs.update(elconf_rev)
    frag_str = Tfragments.substitute(Tdefaults["molecule"], **frag_specs)
    fde_sect = Tfde.substitute(Tdefaults["fde"], **{"method_a": "import_rhoA true", "method_b": "import_rhoB true", "expansion": "SE"})
    extras = "\n".join([extra_basic]+[fde_sect])
    specs_ftmpb = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=frag_str, extras=extras)
    queue_ftmpb = dict(**slrm.weso_small1)  
    # calculation hasn't run yet, create folder in order to copy densmat
    ut.mkdif(meta_ftmpb["path"])
    
    # SPECIAL: copy density matrices
    iterDir = ut.get_last_iter_dir(active="A", path=meta_fnt["path"])  # in case block for MP2_A did not run
    sh.copy(os.path.join(iterDir, "Densmat_B.txt"), os.path.join(meta_ftmpb["path"], "Densmat_A.txt"))
    # corresponds to density of B that is also used in the last cycle
    copy_density(os.path.join(iterDir, "FDE_State0_tot_dens.txt"),
                 os.path.join(meta_ftmpb["path"], "Densmat_B.txt"),
                 header_src=False, alpha_only_src=False)
    
    ut.run_job(specs_ftmpb, queue_ftmpb, meta_ftmpb, Tinp, q_custom=slrm.slurm_add,  
            batch_mode=False)   # because we want to extract data 
    
    # serialize current meta information for later (we're still in the calc folder)
    ut.save_status(meta_ftmpb)
