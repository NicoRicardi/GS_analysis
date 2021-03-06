#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 11:46:46 2021

@author: nico
"""
#import json as js
#import subprocess as sp
#import CCParser as ccp
import numpy as np
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
big = False if len(sys.argv) < 3 else True if sys.argv[2].lower() == "true" else False
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
# 1: HF+MP2 isolated (B)
#-----------------------------------------------------------------------------#

meta_B_MP2  = dict(method="MP2", status=None, basename="mp2")

# get name of calculation folder
meta_B_MP2["path"] = os.path.join(systfol, "B_MP2")
already_done = ut.status_ok(path=meta_B_MP2["path"])
batch_mode_bmp = False
# run calculation and update status ("checkpoint")
if already_done == False:
    memory = 87000
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory})
    mol_specs = dict(xyz=frags["B"])
    if found_elconf:
        mol_specs.update(elconf_B)
    molecule = Tmolecule.substitute(Tdefaults["molecule"], **mol_specs)
    specs_B_MP2 = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=molecule, extras=extra_basic)
    queue_B_MP2 = dict(** slrm.weso_small1)  
    ut.run_job(specs_B_MP2, queue_B_MP2, meta_B_MP2, Tinp,  
            batch_mode=batch_mode_bmp)  # because we want to extract data and copy matrices
    
    # save status
    ut.save_status(meta_B_MP2)
    #-----------------------------------------------------------------------------#
    # Saving HF and MP density matrices
    #-----------------------------------------------------------------------------#
    if not batch_mode_bmp:
        copy_density(os.path.join(meta_B_MP2["path"], "Densmat_SCF.txt"), 
                    os.path.join(systfol, "Densmat_B_nopp.txt"),
                    header_src=False, alpha_only_src=False)
        copy_density(os.path.join(meta_B_MP2["path"], "Densmat_MP.txt"), 
                os.path.join(systfol, "Densmat_B_MP.txt"),
                header_src=False, alpha_only_src=False)
    del specs_B_MP2, queue_B_MP2, batch_mode_bmp
del meta_B_MP2

#-----------------------------------------------------------------------------#
# 2: HF+MP2 (AB)
#-----------------------------------------------------------------------------#
meta_AB_MP2  = dict(method="MP2", status=None, basename="mp2")

# get name of calculation folder
meta_AB_MP2["path"] = os.path.join(systfol, "AB_MP2")
already_done = ut.status_ok(path=meta_AB_MP2["path"])
batch_mode_mpab = True

# run calculation and update status ("checkpoint")
if already_done == False:
    memory = 375000 if big else  240000
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory})
    mol_specs = dict(xyz=frags["AB"])
    if found_elconf:
        mol_specs.update(elconf_AB)
    molecule = Tmolecule.substitute(Tdefaults["molecule"], **mol_specs)
    specs_AB_MP2 = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=molecule, extras=extra_basic)
    
    queue_AB_MP2 = dict(**slrm.weso_big2 if big else slrm.weso_big3)  
    ut.run_job(specs_AB_MP2, queue_AB_MP2, meta_AB_MP2, Tinp,  
            batch_mode=batch_mode_mpab)  # because we want to extract data and copy matrices
    
    # save status
    ut.save_status(meta_AB_MP2)
    if not batch_mode_mpab:
        copy_density(os.path.join(meta_AB_MP2["path"], "Densmat_SCF.txt"),
                os.path.join(systfol, "Densmat_AB_HF.txt"),
                header_src=False, alpha_only_src=False)
        copy_density(os.path.join(meta_AB_MP2["path"], "Densmat_MP.txt"),
                os.path.join(systfol, "Densmat_AB_MP.txt"),
                header_src=False, alpha_only_src=False)
    del specs_AB_MP2, queue_AB_MP2, batch_mode_mpab
del meta_AB_MP2

#-----------------------------------------------------------------------------#
# 3: HF+MP2 ghost (A)
#-----------------------------------------------------------------------------#
meta_A_MP2g  = dict(method="MP2", status=None, basename="mp2_gh")

# get name of calculation folder
meta_A_MP2g["path"] = os.path.join(systfol, "A_MP2_gh")
already_done = ut.status_ok(path=meta_A_MP2g["path"])

# run calculation and update status ("checkpoint")
if already_done == False:
    memory = 87000
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory})
    mol_specs = dict(xyz=frags["AB_ghost"])
    if found_elconf:
        mol_specs.update(elconf_A)
    molecule = Tmolecule.substitute(Tdefaults["molecule"], **mol_specs)
    specs_A_MP2g = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=molecule, extras=extra_basic)
    queue_A_MP2g = dict(** slrm.weso_small1)  
    ut.run_job(specs_A_MP2g, queue_A_MP2g, meta_A_MP2g, Tinp,
            batch_mode=True)
    
    # save status
    ut.save_status(meta_A_MP2g)
    del specs_A_MP2g, queue_A_MP2g
del meta_A_MP2g

#-----------------------------------------------------------------------------#
# 4: HF+MP2 ghost (B)
#-----------------------------------------------------------------------------#
meta_B_MP2g  = dict(method="MP2", status=None, basename="mp2_gh")

# get name of calculation folder
meta_B_MP2g["path"] = os.path.join(systfol, "B_MP2_gh")
already_done = ut.status_ok(path=meta_B_MP2g["path"])

# run calculation and update status ("checkpoint")
if already_done == False:
    memory = 87000
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory})
    mol_specs = dict(xyz=frags["BA_ghost"])
    if found_elconf:
        mol_specs.update(elconf_B)
    molecule = Tmolecule.substitute(Tdefaults["molecule"], **mol_specs)
    specs_B_MP2g = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=molecule, extras=extra_basic)
    queue_B_MP2g = dict(** slrm.weso_small1)  
    ut.run_job(specs_B_MP2g, queue_B_MP2g, meta_B_MP2g, Tinp,  
            batch_mode=True)
    
    # save status
    ut.save_status(meta_B_MP2g)
    del specs_B_MP2g, queue_B_MP2g
del meta_B_MP2g

#-----------------------------------------------------------------------------#    
# 5: HF+MP2 isolated (A) getting Mulliken and ChelPG charges
#-----------------------------------------------------------------------------#
meta_A_MP2  = dict(method="MP2", status=None, basename="mp2")

# get name of calculation folder
meta_A_MP2["path"] = os.path.join(systfol, "A_MP2")
already_done = ut.status_ok(path=meta_A_MP2["path"])

# run calculation and update status ("checkpoint")
if already_done == False:
    memory = 87000
    rem_extras = rem_extras_basic + [chelpg_kw]
    rem_kw_tmp = ut.myupd(rem_adc_basic, {"rem_extras": "\n".join(rem_extras), "memory": memory})
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_kw_tmp)
    
    mol_specs = dict(xyz=frags["A"])
    if found_elconf:
        mol_specs.update(elconf_A)
    molecule = Tmolecule.substitute(Tdefaults["molecule"], **mol_specs)
    specs_A_MP2 = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=molecule, extras=extra_basic)
    queue_A_MP2 = dict(** slrm.weso_small1)  
    json_files = gl.glob(os.path.join(meta_A_MP2["path"],"*.json"))
    ut.run_job(specs_A_MP2, queue_A_MP2, meta_A_MP2, Tinp,
            batch_mode=False)  # because we want to extract data and copy matrices
    try:
        njsf = [os.path.basename(i) for i in gl.glob(os.path.join(meta_A_MP2["path"],"*.json")) if i not in json_files][0]  # whatever json was just added, i.e. default_ccpjson
        assert njsf == default_ccpjson
    except IndexError:
        ccjlog.critical("The json file was already there. Will use default name: {}".format(default_ccpjson))
    except AssertionError:
        ccjlog.critical("CCParser's default seems to be \"{}\" change default_ccpjson in this script!!".format(njsf))
        default_ccpjson = njsf
    # save status
    ut.save_status(meta_A_MP2)
    #-----------------------------------------------------------------------------#
    # Saving HF and MP density matrices
    #-----------------------------------------------------------------------------#
    copy_density(os.path.join(meta_A_MP2["path"], "Densmat_SCF.txt"),
                 os.path.join(systfol, "Densmat_A_nopp.txt"),
                 header_src=False, alpha_only_src=False)
    copy_density(os.path.join(meta_A_MP2["path"], "Densmat_MP.txt"), 
                 os.path.join(systfol, "Densmat_A_MP.txt"),
                 header_src=False, alpha_only_src=False)
    del specs_A_MP2, queue_A_MP2


#-----------------------------------------------------------------------------#
# 6: HF+MP prepol Mulliken(B)
#-----------------------------------------------------------------------------#
meta_B_MPp  = dict(method="HF", status=None, basename="mp_prepol")

# get name of calculation folder
meta_B_MPp["path"] = os.path.join(systfol, "B_MP_pp_Mulliken")
already_done = ut.status_ok(path=meta_B_MPp["path"])

    #-----------------------------------------------------------------------------#
    # Getting MP2 isolated (A) Mulliken and ChelPG charges
    #-----------------------------------------------------------------------------#
parser_jsfile = os.path.join(meta_A_MP2["path"], default_ccpjson)
results = ut.load_js(parser_jsfile)
if type(results["xyz"][0][0]) == str:
    npz = np.load(os.path.join(meta_A_MP2["path"], results["xyz"][0][0]), 
                  allow_pickle=True)
    mulliken_charges = np.hstack([npz["xyz"][0][:,1:],npz["mulliken"][1][:,1:]])  # the first one is MP2 charges
else:
    mulliken_charges = list(map(lambda x, y: x[1:]+[y[1]], results["xyz"][0][0],
                         results["mulliken"][1][0])) # the first one is MP2 charges
mulliken_charges_str = "\n".join(["    ".join(list(map(str, s))) for s in \
                               mulliken_charges])

# run calculation and update status ("checkpoint")
if already_done == False and len(mulliken_charges) != 0:
    memory = 87000
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory}) #asdf
    mol_specs = dict(xyz=frags["B"])
    if found_elconf:
        mol_specs.update(elconf_B)
    molecule = Tmolecule.substitute(Tdefaults["molecule"], **mol_specs)
    point_charges = Tpc.substitute(dict(point_charges=mulliken_charges_str))
    extras = [extra_basic] + [point_charges] 
    inp1 = Tinp.substitute(Tdefaults["inp"], rem_kw=rem_kw, molecule=molecule, extras="\n".join(extras))
    rem_extras = rem_extras_basic + ["max_scf_cycles = 0", "scf_guess = read"]
    rem_kw_tmp = ut.myupd(rem_hf_basic, {"memory": memory, "rem_extras": "\n".join(rem_extras)})
    rem_kw2 = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_kw_tmp)
    inp2 = Tinp.substitute(Tdefaults["inp"], rem_kw=rem_kw2, molecule="read", extras=extra_basic)
    
    specs_B_MPp = dict(inp1=inp1, inp2=inp2)
    queue_B_MPp = dict(** slrm.weso_small1)  
    # technically this calculation always fails. Nonetheless this template
    # works in our favour since CCParser does not know the second part fails
    ut.run_job(specs_B_MPp, queue_B_MPp, meta_B_MPp, Tinps,  
            batch_mode=False, q_custom=slrm.slurm_add)  # because we want to extract data and copy matrices
    
    # save status
    ut.save_status(meta_B_MPp)
    copy_density(os.path.join(meta_B_MPp["path"], "Densmat_SCF.txt"),
            "Densmat_B_pp_Mulliken.txt",
            header_src=False, alpha_only_src=False)
    del specs_B_MPp, queue_B_MPp
del meta_B_MPp

#-----------------------------------------------------------------------------#
#  7: HF+MP prepol ChelPG(B)
#-----------------------------------------------------------------------------#
if "results" not in globals():
    parser_jsfile = os.path.join(meta_A_MP2["path"], default_ccpjson)
    results = ut.load_js(parser_jsfile)
    #-----------------------------------------------------------------------------#
    # Getting MP2 isolated (A) Mulliken and ChelPG charges
    #-----------------------------------------------------------------------------#
if type(results["xyz"][0][0]) == str:
    npz = np.load(os.path.join(meta_A_MP2["path"], results["xyz"][0][0]), 
                  allow_pickle=True)
    chelpg_charges = np.hstack([npz["xyz"][0][:,1:],npz["chelpg"][0][:,1:]])
else:    
    chelpg_charges = list(map(lambda x, y: x[1:]+[y[1]], results["xyz"][0][0],
                         results["chelpg"][0][0]))
chelpg_charges_str = "\n".join(["    ".join(list(map(str, s))) for s in \
                               chelpg_charges])
meta_B_MPp  = dict(method="HF", status=None, basename="mp_prepol")

# get name of calculation folder
meta_B_MPp["path"] = os.path.join(systfol, "B_MP_pp_ChelPG")
already_done = ut.status_ok(path=meta_B_MPp["path"])

# run calculation and update status ("checkpoint")
if already_done == False and len(chelpg_charges) != 0:
    memory = 87000
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory})
    mol_specs = dict(xyz=frags["B"])
    if found_elconf:
        mol_specs.update(elconf_B)
    molecule = Tmolecule.substitute(Tdefaults["molecule"], **mol_specs)
    point_charges = Tpc.substitute(dict(point_charges=chelpg_charges_str))
    extras = [extra_basic] + [point_charges] 
    inp1 = Tinp.substitute(Tdefaults["inp"], rem_kw=rem_kw, molecule=molecule, extras="\n".join(extras))
    rem_extras = rem_extras_basic + ["max_scf_cycles = 0", "scf_guess = read"]
    rem_kw_tmp = ut.myupd(rem_hf_basic, {"memory": memory, "rem_extras": "\n".join(rem_extras)})
    rem_kw2 = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_kw_tmp)
    inp2 = Tinp.substitute(Tdefaults["inp"], rem_kw=rem_kw2, molecule="read", extras=extra_basic)

    specs_B_MPp = dict(inp1=inp1, inp2=inp2)
    queue_B_MPp = dict(** slrm.weso_small1)  
    # technically this calculation always fails. Nonetheless this template
    # works in our favour since CCParser does not know the second part fails
    ut.run_job(specs_B_MPp, queue_B_MPp, meta_B_MPp, Tinps,  
            batch_mode=False, q_custom=slrm.slurm_add)  # because we want to extract data and copy matrices
    
    # save status
    ut.save_status(meta_B_MPp)
    copy_density(os.path.join(meta_B_MPp["path"], "Densmat_SCF.txt"),
            "Densmat_B_pp_ChelPG.txt",
            header_src=False, alpha_only_src=False)
    del specs_B_MPp, queue_B_MPp
del meta_B_MPp
del meta_A_MP2  # only now because needed for charges
#-----------------------------------------------------------------------------#
#  8: Freeze-and-Thaw, ME
#-----------------------------------------------------------------------------#
memory = 87000 if big else 57500
rem_kw = dict(**rem_hf_basic, **{"memory": memory})
fde_kw = Tdefaults["fde"].copy()
specs_fnt = dict(rem_kw=rem_kw, fde_kw=fde_kw, extras=extra_basic, use_zr=False,
                fragments=frags, elconf=elconf, q_custom=slrm.slurm_add,
                maxiter=20, thresh=1e-9, en_file="energies.txt")  
queue_fnt = dict(**slrm.weso_small1 if big else slrm.shabug_L)  
meta_fnt  = dict(method_A="HF", method_B="HF", opt="freeze-thaw", status=None)

# create calculation folder
meta_fnt["path"] = os.path.join(systfol, "FT-ME")
en_file = os.path.join(meta_fnt["path"], specs_fnt["en_file"])
already_done_fnt = ut.status_ok(path=meta_fnt["path"])

if already_done_fnt == False:
    # run calculation and update status ("checkpoint")
    ut.mkdif(meta_fnt["path"])
    ut.logchdir(ccjlog,meta_fnt["path"])
    update = False
    if os.path.exists(os.path.join(systfol, "Densmat_B_nopp.txt")):
        fde_kw.update(method_a="import_rhoA true")
        update = True
        if not os.path.exists(os.path.join(meta_fnt["path"],"Densmat_B.txt")):
            sh.copy(os.path.join(systfol, "Densmat_B_nopp.txt"), "Densmat_B.txt")
    if os.path.exists(os.path.join(systfol, "Densmat_A_nopp.txt")):
        fde_kw.update(method_b="import_rhoB true")
        update = True
        if not os.path.exists(os.path.join(meta_fnt["path"],"Densmat_A.txt")):
            sh.copy(os.path.join(systfol, "Densmat_A_nopp.txt"), "Densmat_A.txt")
    if update:
        specs_fnt.update(fde_kw=fde_kw)
    try:
        freeze_and_thaw(queue_fnt, **specs_fnt)  
        meta_fnt["status"] = "FIN"
        already_done_fnt = True
    except Exception as e:
        ccjlog.error("{}\n{}".format(e,traceback.format_exc()))
        meta_fnt["status"] = "FAIL"
    # save status
    ut.save_status(meta_fnt)
#-----------------------------------------------------------------------------#
#  9: FDE-MP2 using FT densities: (B-in-A) [ME]
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
    fde_sect = Tfde.substitute(Tdefaults["fde"], **{"method_a": "import_rhoA true", "method_b": "import_rhoB true"})
    extras = "\n".join([extra_basic]+[fde_sect])
    specs_ftmpb = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=frag_str, extras=extras)
    queue_ftmpb = dict(** slrm.weso_small1)  
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
            batch_mode=True)   # because we want to extract data 
    
    # save status
    ut.save_status(meta_ftmpb)
    
#-----------------------------------------------------------------------------#    
# 10: FDE-MP2 using FT densities: (A-in-B) [ME]
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
    fde_sect = Tfde.substitute(Tdefaults["fde"], **{"method_a": "import_rhoA true", "method_b": "import_rhoB true"})
    extras = "\n".join([extra_basic]+[fde_sect])
    specs_ftmpa = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=frag_str, extras=extras)
    queue_ftmpa = dict(** slrm.weso_small1)  
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
            batch_mode=False)  # because we want to extract data  
    try:
        njsf = [os.path.basename(i) for i in gl.glob(os.path.join(meta_ftmpa["path"],"*.json")) if i not in json_files][0]  # whatever json was just added, i.e. default_ccpjson
        assert njsf == default_ccpjson
    except IndexError:
        ccjlog.critical("The json file was already there. Will use default name: {}".format(default_ccpjson))
    except AssertionError:
        ccjlog.critical("CCParser's default seems to be \"{}\" change default_ccpjson in this script!!".format(njsf))
        default_ccpjson = njsf
    # save status
    ut.save_status(meta_ftmpa)
    parser_jsfile = os.path.join(meta_ftmpa["path"], default_ccpjson)
    data = ut.load_js(parser_jsfile)  # default ccp json_filename
    sp.call("echo {E_A} >> {en_file}".format(E_A=data["scf_energy"][-1][0], en_file=en_file), shell=True)


#-----------------------------------------------------------------------------#
#  11+: Macrocycles for rhoA dependency with frozen rhoB for several rhoB
#-----------------------------------------------------------------------------#
densities = {i: "Densmat_B_{}.txt".format(i) for i in ["nopp", "pp_Mulliken", "pp_ChelPG"]}
mpb_fols = {"nopp": os.path.join(systfol, "B_MP2"),
            "pp_Mulliken": os.path.join(systfol, "B_MP_pp_Mulliken"),
            "pp_ChelPG": os.path.join(systfol, "B_MP_pp_ChelPG")}

memory = 87000 if big else 57500
rem_kw = dict(**rem_hf_basic, **{"memory": memory})
fde_kw = Tfde.substitute(Tdefaults["fde"], **{"method_a": "import_rhoA true", "method_b": "import_rhoB true"})
specs_mc = dict(rem_kw=rem_kw, fde_kw=fde_kw, extras=extra_basic, use_zr=False,
                fragments=frags, elconf=elconf, q_custom=slrm.slurm_add,
                maxiter=20, thresh=1e-9, en_file="energies.txt")  
meta_mc  = dict(method_A="HF", method_B="HF", opt="macrocycles", status=None)
queue_mc = dict(**slrm.weso_small1 if big else slrm.shabug_L)  
for ID, dmfile in densities.items():
    meta_mc["path"] = os.path.join(systfol, "MC-{}".format(ID))
    en_file = os.path.join(meta_mc["path"], specs_mc["en_file"])
    already_done_mc = ut.status_ok(path=meta_mc["path"])
    if already_done_mc == False:
        try:
            ut.mkdif(meta_mc["path"])
            ut.logchdir(ccjlog,meta_mc["path"])
            ##exceptional, not general
            if not os.path.exists(os.path.join(meta_mc["path"],"Densmat_B.txt")):
                sh.copy(os.path.join(systfol, dmfile), "Densmat_B.txt")
            if not os.path.exists(os.path.join(meta_mc["path"],"Densmat_A.txt")):
                sh.copy(os.path.join(systfol, "Densmat_A_nopp.txt"), "Densmat_A.txt")
                #-----------------------------------------------------------------------------#
                # HF Macrocycles 
                #-----------------------------------------------------------------------------#
            macrocycles(queue_mc, **specs_mc)
            meta_mc["status"] = "FIN"
            already_done_mc = True
        except Exception as e:  # Mainly NotConverged, but not only
            ccjlog.error("{}\n{}".format(e,traceback.format_exc()))
            meta_mc["status"] = "FAIL"
        ut.save_status(meta_mc)
    #-----------------------------------------------------------------------------#
    # FDE-MP2  A in B
    #-----------------------------------------------------------------------------#   
    meta_mpa = dict(method_A="MP2", method_B="import", opt=None, status=None,basename="emb")
    meta_mpa["path"] = os.path.join(meta_mc["path"], "MP2_A")
    already_done_mpa = ut.status_ok(path=meta_mpa["path"])
    if already_done_mpa == False and already_done_mc == True:
        memory = 87000
        rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory, "fde": "true"})
        frag_specs = dict(frag_a=frags["A"], frag_b=frags["B"])
        if found_elconf:
            frag_specs.update(elconf)
        frag_str = Tfragments.substitute(Tdefaults["molecule"], **frag_specs)
        fde_sect = Tfde.substitute(Tdefaults["fde"], **{"method_a": "import_rhoA true", "method_b": "import_rhoB true"})
        extras = "\n".join([extra_basic]+[fde_sect])
        specs_mpa = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=frag_str, extras=extras)
        queue_mpa = dict(**slrm.weso_small1)  
        ut.mkdif(meta_mpa["path"]) 
        iterDir = ut.get_last_iter_dir(active="A", path=meta_mc["path"],opt="macrocycles")
        sh.copy(os.path.join(iterDir, "Densmat_B.txt"), meta_mpa["path"])
        copy_density(os.path.join(iterDir, "FDE_State0_tot_dens.txt"),
                    os.path.join(meta_mpa["path"], "Densmat_A.txt"),
                    header_src=False, alpha_only_src=False)
        json_files = gl.glob(os.path.join(meta_mpa["path"],"*.json"))
        ut.run_job(specs_mpa, queue_mpa, meta_mpa, Tinp, q_custom=slrm.slurm_add, batch_mode=False)  # because we want to extract data
        
        try:
            njsf = [os.path.basename(i) for i in gl.glob(os.path.join(meta_mpa["path"],"*.json")) if i not in json_files][0]
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
    dst = os.path.join(meta_mc["path"], "MP2_B")
    symlinkif(mpb_fols[ID], dst)
