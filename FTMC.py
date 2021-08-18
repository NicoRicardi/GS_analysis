#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 18:17:22 2021

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

def symlinkif(src, dst, printout=False):
    """
    Note
    ----
    as os.symlink but never raises errors (checks if file already exists)
    
    Parameters
    ----------
    src: str
        source path
    dst: str
        destination path
    printout: bool
        whether it should print what happens
    """
    if os.path.exists(src) and not os.path.exists(dst):
        os.symlink(src, dst)
        if printout:
            print("created symlink!")
    if printout:
        print("Nothing done")

###############################################################################
###############################################################################
root = os.getcwd()
system = sys.argv[1]
expansion = sys.argv[2]
systfol = os.path.join(root, system)
###############################################################################
logfile = os.path.join(systfol,"{}_ccj.log".format(system))
# set up logger. NB avoid homonimity with other module's loggers (e.g. ccp)
ccjlog = logging.getLogger("ccj")
#ut.setupLogger(to_console=True, to_log=True, printlevel=20)
ut.setupLogger(to_console=True, to_log=True, logname=logfile)

default_ccpjson = "CCParser.json"
bs_kw = "gen"  
bs_string = ut.read_file("aug-cc-pVDZ.bas")
if bs_string[-1] == "\n":
    bs_string = bs_string[:-1]
rem_extras_basic = ["thresh = 14", "basis_lin_dep_thresh = 5"]
rem_adc_basic = dict(method="adc(2)", basis=bs_kw,
                     ee_specs=Tadc.substitute(Tdefaults["adc"]),
                     rem_extras="\n".join(rem_extras_basic))
rem_hf_basic = dict(method="hf", basis=bs_kw,
                    rem_extras="\n".join(rem_extras_basic))
extra_basic = Tbasis.substitute(**{"basis_specs": bs_string})
    
ut.logchdir(ccjlog,systfol)
# --- First things first : STRUCTURE ---
zr_file = ccj.find_file(systfol, extension="zr")
frags = ccj.zr_frag(zr_file)

frags_rev = {}
frags_rev["A"] = frags["B"]
frags_rev["B"] = frags["A"]
frags_rev["AB_ghost"] = frags["BA_ghost"]
frags_rev["AB"] = frags_rev["A"] + frags_rev["B"]

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
# (3.1) Macrocycles starting from each FnT iteration
#-----------------------------------------------------------------------------#  
memory = 57500 if expansion == "ME" else 87000
rem_kw = dict(**rem_hf_basic, **{"memory": memory})
fde_kw = Tfde.substitute(Tdefaults["fde"], **{"method_a": "import_rhoA true", "method_b": "import_rhoB true", "expansion": expansion})
specs_ftmc = dict(rem_kw=rem_kw, fde_kw=fde_kw, extras=extra_basic, use_zr=False,
                fragments=frags, elconf=elconf, q_custom=slrm.slurm_add,
                maxiter=20, thresh=1e-9, en_file="energies.txt")  
specs_ftmc_rev = dict(rem_kw=rem_kw, fde_kw=fde_kw, extras=extra_basic, use_zr=False,
                fragments=frags_rev, elconf=elconf_rev, q_custom=slrm.slurm_add,
                maxiter=20, thresh=1e-9, en_file="energies.txt")  
meta_ftmc  = dict(method_A="HF", method_B="HF", opt="macrocycles", status=None)
queue_ftmc = dict(**slrm.shabug_L if expansion == "SE" else slrm.weso_small1)  
        
ft_fol = os.path.join(systfol, "FT-{}".format(expansion))
dirbase = "cy"
iterdirs = [i for i in os.listdir(ft_fol) if dirbase in i]
# we run for both odd- and even-numbered cycles because we need MP2_B as well
# we stop before the last two because we already have them from the end of FT 
for n in range(0, len(iterdirs)-2):  
    cyfol = "{}{}".format(dirbase, n)
    ftmc_fol = os.path.join(ft_fol, cyfol, "FT{}-MC-{}".format(n, expansion))
    meta_ftmc["path"] = ftmc_fol
    if expansion == "ME" and n == 0:  # we already have this, just symlink
        symlinkif(os.path.join(systfol, "MC-nopp"), ftmc_fol)
        continue
    en_file = os.path.join(ftmc_fol, specs_ftmc["en_file"])
    already_done_ftmc = ut.status_ok(path=meta_ftmc["path"])
    if already_done_ftmc == False:
        ut.mkdif(ftmc_fol)
        if not os.path.isfile(os.path.join(ftmc_fol, "Densmat_B.txt")):
            sh.copy(os.path.join(ft_fol, cyfol, "Densmat_B.txt"),
                    os.path.join(ftmc_fol, "Densmat_B.txt"))
        if not os.path.isfile(os.path.join(ftmc_fol, "Densmat_A.txt")):
            if os.path.isfile(os.path.join(ft_fol, cyfol, "Densmat_A.txt")):
                sh.copy(os.path.join(ft_fol, cyfol, "Densmat_A.txt"),
                        os.path.join(ftmc_fol, "Densmat_A.txt"))
            else:
                copy_density(os.path.join(ft_fol, cyfol, "FDE_State0_tot_dens.txt"),
                             os.path.join(ftmc_fol, "Densmat_A.txt"),
                             header_src=False,
                             alpha_only_src=False)
        try:
            ut.logchdir(ccjlog,ftmc_fol)
            macrocycles(queue_ftmc, **specs_ftmc if n%2 == 0 else specs_ftmc_rev)
            meta_ftmc["status"] = "FIN"
            already_done_ftmc = True
        except Exception as e:  # Mainly NotConverged, but not only
            print(e)
            traceback.print_exc()
            meta_ftmc["status"] = "FAIL"
        ut.save_status(meta_ftmc)
    #-----------------------------------------------------------------------------#
    # FDE-MP2  A in B
    #-----------------------------------------------------------------------------#   
    meta_mpa = dict(method_A="MP2", method_B="import", opt=None, status=None,basename="emb")
    meta_mpa["path"] = os.path.join(meta_ftmc["path"], "MP2_A")
    already_done_mpa = ut.status_ok(path=meta_mpa["path"])
    if already_done_mpa == False and already_done_ftmc == True:
        memory = 87000 if expansion == "ME" else 712000
        rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory, "fde": "true"})
        frag_specs = dict(frag_a=frags["A"], frag_b=frags["B"]) if n%2 == 0 else dict(frag_a=frags_rev["A"], frag_b=frags_rev["B"])
        if found_elconf:
            frag_specs.update(elconf if n%2 == 0 else elconf_rev)
        frag_str = Tfragments.substitute(Tdefaults["molecule"], **frag_specs)
        fde_sect = Tfde.substitute(Tdefaults["fde"], **{"method_a": "import_rhoA true", "method_b": "import_rhoB true", "expansion": expansion})
        extras = "\n".join([extra_basic]+[fde_sect])
        specs_mpa = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=frag_str, extras=extras)
        queue_mpa = dict(**slrm.weso_small1 if expansion == "ME" else slrm.weso_big1)  
        ut.mkdif(meta_mpa["path"]) 
        iterDir = ut.get_last_iter_dir(active="A", path=meta_ftmc["path"],opt="macrocycles")
        sh.copy(os.path.join(iterDir, "Densmat_B.txt"), meta_mpa["path"])
        copy_density(os.path.join(iterDir, "FDE_State0_tot_dens.txt"),
                    os.path.join(meta_mpa["path"], "Densmat_A.txt"),
                    header_src=False, alpha_only_src=False)
        json_files = gl.glob(os.path.join(meta_mpa["path"],"*.json"))
        ut.run_job(specs_mpa, queue_mpa, meta_mpa, Tinp, q_custom=slrm.slurm_add, batch_mode=False)  # because we want to extract data
        try:
            njsf = [os.path.basename(i) for i in gl.glob(os.path.join(meta_mpa["path"],"*.json")) if i not in json_files][0]  # whatever json was just added, i.e. default_ccpjson
            assert njsf == default_ccpjson
        except IndexError:
            ccjlog.critical("The json file was already there. Will use default name: {}".format(default_ccpjson))
        except AssertionError:
            ccjlog.critical("CCParser's default seems to be \"{}\" change default_ccpjson in this script!!".format(njsf))
            default_ccpjson = njsf
        parser_jsfile = os.path.join(meta_mpa["path"], default_ccpjson)
        data = ut.load_js(parser_jsfile)  # default ccp json_filename
        ut.save_status(meta_mpa)
        sp.call("echo {E_A} >> {en_file}".format(E_A=data["scf_energy"][-1][0], en_file=en_file), shell=True)
    #-----------------------------------------------------------------------------#
    # Create Symbolic link  to MP2 B
    #-----------------------------------------------------------------------------#    
    if n % 2 == 0 and n > 1:
        dst = os.path.join(meta_ftmc["path"], "MP2_B")
        src = os.path.join(ft_fol, "{}{}".format(dirbase, n-1),
                               "FT{}-MC-{}".format(n-1, expansion),
                               "MP2_A")
        symlinkif(src, dst)