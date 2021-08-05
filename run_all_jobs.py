#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 11:46:46 2021

@author: nico
"""
#import json as js
#import subprocess as sp
#import CCParser as ccp
#import numpy as np
import glob as gl
import os
import CCJob as ccj
import CCJob.utils as ut
from CCJob.iterative import freeze_and_thaw, macrocycles, copy_density, NotConvergedError
import shutil as sh
from CCJob.Composable_templates import Tdefaults, Tinps, Tinp, Trem_kw, Tmolecule,\
 Tadc,  Tpc, Tbasis, chelpg_kw, Tfragments, Tfde
import slurm_stuff_yggdrasil as slrm
#import slurm_stuff_baobab as slrm


###############################################################################
##############################################################################
# --- First things first : STRUCTURE ---
systfol = os.getcwd()
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
    
bs_kw = "gen"  
bs_string = ut.read_file("aug-cc-pVDZ.bas")
rem_adc_basic = dict(method="adc(2)", basis=bs_kw, ee_specs=Tadc.substitute(Tdefaults["adc"]))
rem_hf_basic = dict(method="hf", basis=bs_kw)
extra_basic = [Tbasis.substitute(**{"basis_specs": bs_string})]
#-----------------------------------------------------------------------------#
# REF 1: HF+MP2 isolated (A) getting Mulliken and ChelPG charges
#-----------------------------------------------------------------------------#
meta_A_MP2  = dict(method="MP2", status=None, basename="mp2")

# get name of calculation folder
meta_A_MP2["path"] = os.path.join(systfol, "A_MP2")
already_done = ut.status_ok(path=meta_A_MP2["path"])

# run calculation and update status ("checkpoint")
if not already_done:
    memory = 14000
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"rem_extras": chelpg_kw, "memory": memory})
    
    mol_specs = dict(xyz=frags["A"])
    if found_elconf:
        mol_specs.update(elconf_A)
    molecule = Tmolecule.substitute(Tdefaults["molecule"], **mol_specs)
    specs_A_MP2 = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=molecule, extras=extra_basic)
    queue_A_MP2 = dict(**slrm.shabug_XS)  
    json_files = gl.glob("*.json")
    ut.run_job(specs_A_MP2, queue_A_MP2, meta_A_MP2, Tinp,
            batch_mode=True)
    njsf = [i for i in gl.glob("*.json") if i not in json_files][0]  # whatever json was just added, i.e. default_ccpjson
     
    # serialize current meta information for later (we're still in the calc folder)
    ut.save_status(meta_A_MP2)
    #-----------------------------------------------------------------------------#
    # Saving MP2 isolated (A) Mulliken and ChelPG charges
    #-----------------------------------------------------------------------------#
    parser_jsfile = os.path.join(meta_A_MP2["path"], njsf)
    results = ut.load_js(parser_jsfile)
    mulliken_charges = list(map(lambda x, y: x[1:]+[y[1]], results["xyz"][0][0],
                             results["mulliken"][0][0]))
    mulliken_charges_str = "\n".join(["    ".join(list(map(str, s))) for s in \
                                   mulliken_charges])
    chelpg_charges = list(map(lambda x, y: x[1:]+[y[1]], results["xyz"][0][0],
                             results["chelpg"][0][0]))
    chelpg_charges_str = "\n".join(["    ".join(list(map(str, s))) for s in \
                                   chelpg_charges])
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
del meta_A_MP2

#-----------------------------------------------------------------------------#
# REF 2: HF+MP2 ghost (A)
#-----------------------------------------------------------------------------#
meta_A_MP2g  = dict(method="MP2", status=None, basename="mp2_gh")

# get name of calculation folder
meta_A_MP2g["path"] = os.path.join(systfol, "A_MP2_gh")
already_done = ut.status_ok(path=meta_A_MP2g["path"])

# run calculation and update status ("checkpoint")
if not already_done:
    memory = 14000
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory})
    mol_specs = dict(xyz=frags["AB_ghost"])
    if found_elconf:
        mol_specs.update(elconf_A)
    molecule = Tmolecule.substitute(Tdefaults["molecule"], **mol_specs)
    specs_A_MP2g = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=molecule, extras=extra_basic)
    queue_A_MP2g = dict(**slrm.shabug_XS)  
    ut.run_job(specs_A_MP2g, queue_A_MP2g, meta_A_MP2g, Tinp,
            batch_mode=True)
    
    # serialize current meta information for later (we're still in the calc folder)
    ut.save_status(meta_A_MP2g)
    del specs_A_MP2g, queue_A_MP2g
del meta_A_MP2g

#-----------------------------------------------------------------------------#
# REF 3: HF+MP2 isolated (B)
#-----------------------------------------------------------------------------#

meta_B_MP2  = dict(method="MP2", status=None, basename="mp2")

# get name of calculation folder
meta_B_MP2["path"] = os.path.join(systfol, "B_MP2")
already_done = ut.status_ok(path=meta_B_MP2["path"])

# run calculation and update status ("checkpoint")
if not already_done:
    memory = 14000
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory})
    mol_specs = dict(xyz=frags["B"])
    if found_elconf:
        mol_specs.update(elconf_B)
    molecule = Tmolecule.substitute(Tdefaults["molecule"], **mol_specs)
    specs_B_MP2 = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=molecule, extras=extra_basic)
    queue_B_MP2 = dict(**slrm.shabug_XS)  
    ut.run_job(specs_B_MP2, queue_B_MP2, meta_B_MP2, Tinp,  
            batch_mode=True)
    
    # serialize current meta information for later (we're still in the calc folder)
    ut.save_status(meta_B_MP2)
    #-----------------------------------------------------------------------------#
    # Saving HF and MP density matrices
    #-----------------------------------------------------------------------------#
    copy_density(os.path.join(meta_B_MP2["path"], "Densmat_SCF.txt"), 
                os.path.join(systfol, "Densmat_B_nopp.txt"),
                header_src=False, alpha_only_src=False)
    copy_density(os.path.join(meta_B_MP2["path"], "Densmat_MP.txt"), 
            os.path.join(systfol, "Densmat_B_MP.txt"),
            header_src=False, alpha_only_src=False)
    del specs_B_MP2, queue_B_MP2
del meta_B_MP2
#-----------------------------------------------------------------------------#
# REF 4: HF+MP2 ghost (B)
#-----------------------------------------------------------------------------#
meta_B_MP2g  = dict(method="MP2", status=None, basename="mp2_gh")

# get name of calculation folder
meta_B_MP2g["path"] = os.path.join(systfol, "B_MP2_gh")
already_done = ut.status_ok(path=meta_B_MP2g["path"])

# run calculation and update status ("checkpoint")
if not already_done:
    memory = 14000
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory})
    mol_specs = dict(xyz=frags["BA_ghost"])
    if found_elconf:
        mol_specs.update(elconf_B)
    molecule = Tmolecule.substitute(Tdefaults["molecule"], **mol_specs)
    specs_B_MP2g = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=molecule, extras=extra_basic)
    queue_B_MP2g = dict(**slrm.shabug_XS)  
    ut.run_job(specs_B_MP2g, queue_B_MP2g, meta_B_MP2g, Tinp,  
            batch_mode=True)
    
    # serialize current meta information for later (we're still in the calc folder)
    ut.save_status(meta_B_MP2g)
    del specs_B_MP2g, queue_B_MP2g
del meta_B_MP2g
#-----------------------------------------------------------------------------#
# REF 5: HF+MP2 (AB)
#-----------------------------------------------------------------------------#
meta_AB_MP2  = dict(method="MP2", status=None, basename="mp2")

# get name of calculation folder
meta_AB_MP2["path"] = os.path.join(systfol, "AB_MP2")
already_done = ut.status_ok(path=meta_AB_MP2["path"])

# run calculation and update status ("checkpoint")
if not already_done:
    memory = 14000
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory})
    mol_specs = dict(xyz=frags["AB"])
    if found_elconf:
        mol_specs.update(elconf_AB)
    molecule = Tmolecule.substitute(Tdefaults["molecule"], **mol_specs)
    specs_AB_MP2 = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=molecule, extras=extra_basic)
    
    queue_AB_MP2 = dict(**slrm.shabug_XS)  
    ut.run_job(specs_AB_MP2, queue_AB_MP2, meta_AB_MP2, Tinp,  
            batch_mode=True)
    
    # serialize current meta information for later (we're still in the calc folder)
    ut.save_status(meta_AB_MP2)
    copy_density(os.path.join(meta_AB_MP2["path"], "Densmat_SCF.txt"),
            os.path.join(systfol, "Densmat_AB_HF.txt"),
            header_src=False, alpha_only_src=False)
    copy_density(os.path.join(meta_AB_MP2["path"], "Densmat_MP.txt"),
            os.path.join(systfol, "Densmat_AB_MP.txt"),
            header_src=False, alpha_only_src=False)
    del specs_AB_MP2, queue_AB_MP2
del meta_AB_MP2

#-----------------------------------------------------------------------------#
# REF 6: HF prepol Mulliken(B)
#-----------------------------------------------------------------------------#
meta_B_HFp  = dict(method="HF", status=None, basename="hf_prepol")

# get name of calculation folder
meta_B_HFp["path"] = os.path.join(systfol, "B_HF_ppM_Mulliken")
already_done = ut.status_ok(path=meta_B_HFp["path"])

# run calculation and update status ("checkpoint")
if not already_done and len(mulliken_charges) != 0:
    memory = 14000
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_hf_basic, **{"memory": memory})
    mol_specs = dict(xyz=frags["B"])
    if found_elconf:
        mol_specs.update(elconf_B)
    molecule = Tmolecule.substitute(Tdefaults["molecule"], **mol_specs)
    point_charges = Tpc.substitute(dict(point_charges=mulliken_charges_str))
    extras = [extra_basic] + [point_charges] 
    inp1 = Tinp.substitute(Tdefaults["inp"], rem_kw=rem_kw, molecule=molecule, extras="\n".join(extras))
    rem_extras = "\n".join(["max_scf_cycles = 0", "scf_guess = read"])
    rem_kw2 = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_hf_basic, **{"memory": memory, rem_extras: rem_extras})
    inp2 = Tinp.substitute(Tdefaults["inp"], rem_kw=rem_kw2, molecule="read", extras=extra_basic)
    
    specs_B_HFp = dict(inp1=inp1, inp2=inp2)
    queue_B_HFp = dict(**slrm.shabug_XS)  
    # technically this calculation always fails. Nonetheless this template
    # works in our favour since CCParser does not know the second part fails
    ut.run_job(specs_B_HFp, queue_B_HFp, meta_B_HFp, Tinps,  
            batch_mode=False, q_custom=slrm.slurm_add)  
    
    # serialize current meta information for later (we're still in the calc folder)
    ut.save_status(meta_B_HFp)
    copy_density(os.path.join(meta_B_HFp["path"], "Densmat_SCF.txt"),
            "Densmat_B_pp_Mulliken.txt",
            header_src=False, alpha_only_src=False)
    del specs_B_HFp, queue_B_HFp
del meta_B_HFp

#-----------------------------------------------------------------------------#
# REF 7: HF prepol ChelPG(B)
#-----------------------------------------------------------------------------#

meta_B_HFp  = dict(method="HF", status=None, basename="hf_prepol")

# get name of calculation folder
meta_B_HFp["path"] = os.path.join(systfol, "B_HF_ppM_ChelPG")
already_done = ut.status_ok(path=meta_B_HFp["path"])

# run calculation and update status ("checkpoint")
if not already_done and len(chelpg_charges) != 0:
    memory = 14000
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_hf_basic, **{"memory": memory})
    mol_specs = dict(xyz=frags["B"])
    if found_elconf:
        mol_specs.update(elconf_B)
    molecule = Tmolecule.substitute(Tdefaults["molecule"], **mol_specs)
    point_charges = Tpc.substitute(dict(point_charges=chelpg_charges_str))
    extras = [extra_basic] + [point_charges] 
    inp1 = Tinp.substitute(Tdefaults["inp"], rem_kw=rem_kw, molecule=molecule, extras="\n".join(extras))
    rem_extras = "\n".join(["max_scf_cycles = 0", "scf_guess = read"])
    rem_kw2 = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_hf_basic, **{"memory": memory, rem_extras: rem_extras})
    inp2 = Tinp.substitute(Tdefaults["inp"], rem_kw=rem_kw2, molecule="read", extras=extra_basic)

    specs_B_HFp = dict(inp1=inp1, inp2=inp2)
    queue_B_HFp = dict(**slrm.shabug_XS)  
    # technically this calculation always fails. Nonetheless this template
    # works in our favour since CCParser does not know the second part fails
    ut.run_job(specs_B_HFp, queue_B_HFp, meta_B_HFp, Tinps,  
            batch_mode=False, q_custom=slrm.slurm_add)  
    
    # serialize current meta information for later (we're still in the calc folder)
    ut.save_status(meta_B_HFp)
    copy_density(os.path.join(meta_B_HFp["path"], "Densmat_SCF.txt"),
            "Densmat_B_pp_ChelPG.txt",
            header_src=False, alpha_only_src=False)
    del specs_B_HFp, queue_B_HFp
del meta_B_HFp

#-----------------------------------------------------------------------------#
# REF 8: Freeze-and-Thaw, ME
#-----------------------------------------------------------------------------#
memory = 14000
rem_kw = dict(**rem_hf_basic, **{"memory": memory})
fde_kw = Tfde.substitute(Tdefaults["fde"], **{"method_a": "import_rhoA = true", "method_b": "import_rhoB = true"})
specs_fnt = dict(rem_kw=rem_kw, fde_kw=fde_kw, extras=extra_basic, 
                fragments=frags, elconf=elconf, q_custom=slrm.slurm_add,
                maxiter=10, thresh=1e-9)  
queue_fnt = dict(**slrm.shabug_XS)  
meta_fnt  = dict(method_A="HF", method_B="HF", opt="freeze-thaw", status=None)

# create calculation folder
meta_fnt["path"] = os.path.join(systfol, "FT-ME")
already_done = ut.status_ok(path=meta_fnt["path"])

if not already_done:
    # run calculation and update status ("checkpoint")
    try:
        if not os.path.exists(meta_fnt["path"]):
            os.makedirs(meta_fnt["path"])
        os.chdir(meta_fnt["path"])
        if not os.path.exists(os.path.join(meta_fnt["path"],"Densmat_B.txt")):
            sh.copy(os.path.join(systfol, "Densmat_B_nopp.txt"), "Densmat_B.txt")
        if not os.path.exists(os.path.join(meta_fnt["path"],"Densmat_A.txt")):
            sh.copy(os.path.join(systfol, "Densmat_A_nopp.txt"), "Densmat_A.txt")
        freeze_and_thaw(queue_fnt, **specs_fnt)  
        meta_fnt["status"] = "FIN"
    except NotConvergedError:
        meta_fnt["status"] = "FAIL"
    
    # serialize current meta information for later
    ut.save_status(meta_fnt)
    energies_file = os.path.join(meta_fnt["path"], "energies.json")
    energies = ut.load_js(energies_file)
#-----------------------------------------------------------------------------#    
# REF 9: FDE-MP2 using FT densities: (A-in-B) [ME]
#-----------------------------------------------------------------------------#
meta_ftmpa  = dict(method_A="MP2", method_B="import", opt=None, status=None,
              basename="emb")

# get name of calculation folder
meta_ftmpa["path"] = os.path.join(meta_fnt["path"], "MP2_A")
already_done = ut.status_ok(path=meta_ftmpa["path"])

# run calculation and update status ("checkpoint")
if not already_done and meta_fnt["status"] == "FIN":
    memory = 14000
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory, "fde": "true"})
    frag_specs = dict(frag_a=frags["A"], frag_b=frags["B"])
    if found_elconf:
        frag_specs.update(elconf)
    frag_str = Tfragments.substitute(Tdefaults["molecule"], **frag_specs)
    fde_sect = Tfde.substitute(Tdefaults["fde"], **{"method_a": "import_rhoA = true", "method_b": "import_rhoB = true"})
    extras = "\n".join([extra_basic]+[fde_sect])
    specs_ftmpa = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=frag_str, extras=extras)
    queue_ftmpa = dict(**slrm.shabug_XS)  
    # calculation hasn't run yet, create folder in order to copy densmat
    if not os.path.exists(meta_ftmpa["path"]):
        os.makedirs(meta_ftmpa["path"])
    
    # SPECIAL: copy density matrices
    iterDir = ut.get_last_iter_dir(active="A", path=meta_fnt["path"])
    sh.copy(os.path.join(iterDir, "Densmat_B.txt"), meta_ftmpa["path"],)
    copy_density(os.path.join(iterDir, "FDE_State0_tot_dens.txt"),
                 os.path.join(meta_ftmpa["path"],"Densmat_A.txt"),
                 header_src=False, alpha_only_src=False)
    
    # finally run in batch mode
    json_files = gl.glob("*.json")
    ut.run_job(specs_ftmpa, queue_ftmpa, meta_ftmpa, Tinp, q_custom=slrm.slurm_add,  
            batch_mode=True, create_folder=False)  
    njsf = [i for i in gl.glob("*.json") if i not in json_files][0]  # whatever json was just added, i.e. default_ccpjson
    # serialize current meta information for later (we're still in the calc folder)
    ut.save_status(meta_ftmpa)
    data = ut.load_js(njsf)  # default ccp json_filename
    energies.append(data["SCF"][-1]) # CHECK
    ut.dump_js(energies, energies_file)

#-----------------------------------------------------------------------------#
# REF 10: FDE-MP2 using FT densities: (B-in-A) [ME]
#-----------------------------------------------------------------------------#

meta_ftmpb  = dict(method_A="MP2", method_B="import", opt=None, status=None,
              basename="emb")
# get name of calculation folder
meta_ftmpb["path"] = os.path.join(meta_fnt["path"], "MP2_B")
already_done = ut.status_ok(path=meta_ftmpb["path"])

# run calculation and update status ("checkpoint")
if not already_done and meta_fnt["status"] == "FIN":
    memory = 14000
    rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory, "fde": "true"})
    frag_specs = dict(frag_a=frags["B"], frag_b=frags["A"])
    if found_elconf:
        frag_specs.update(elconf_rev)
    frag_str = Tfragments.substitute(Tdefaults["molecule"], **frag_specs)
    fde_sect = Tfde.substitute(Tdefaults["fde"], **{"method_a": "import_rhoA = true", "method_b": "import_rhoB = true"})
    extras = "\n".join([extra_basic]+[fde_sect])
    specs_ftmpb = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=frag_str, extras=extras)
    queue_ftmpb = dict(**slrm.shabug_XS)  
    # calculation hasn't run yet, create folder in order to copy densmat
    if not os.path.exists(meta_ftmpb["path"]):
        os.makedirs(meta_ftmpb["path"])
    
    # SPECIAL: copy density matrices
#    iterDir = ut.get_last_iter_dir(active="A", path=meta_fnt["path"])
    sh.copy(os.path.join(iterDir, "Densmat_B.txt"), os.path.join(meta_ftmpb["path"], "Densmat_A.txt"))
    # corresponds to density of B that is also used in the last cycle
    copy_density(os.path.join(iterDir, "FDE_State0_tot_dens.txt"),
                 os.path.join(meta_ftmpb["path"], "Densmat_B.txt"),
                 header_src=False, alpha_only_src=False)
    
    ut.run_job(specs_ftmpb, queue_ftmpb, meta_ftmpb, Tinp, q_custom=slrm.slurm_add,  
            batch_mode=True)
    
    # serialize current meta information for later (we're still in the calc folder)
    ut.save_status(meta_ftmpb)

#-----------------------------------------------------------------------------#
# REF 11+: Macrocycles for rhoA dependency with frozen rhoB for several rhoB
#-----------------------------------------------------------------------------#
densities = {i: "Densmat_B_{}.txt".format(i) for i in ["nopp", "pp_Mulliken", "pp_ChelPG"]}

memory = 14000
rem_kw = dict(**rem_hf_basic, **{"memory": memory})
fde_kw = Tfde.substitute(Tdefaults["fde"], **{"method_a": "import_rhoA = true", "method_b": "import_rhoB = true"})
specs_mc = dict(rem_kw=rem_kw, fde_kw=fde_kw, extras=extra_basic, 
                fragments=frags, elconf=elconf, q_custom=slrm.slurm_add,
                maxiter=10, thresh=1e-9)  
meta_mc  = dict(method_A="HF", method_B="HF", opt="macrocycles", status=None)
queue_mc = dict(**slrm.shabug_XS)  
for ID, dmfile in densities.items():
    meta_mc["path"] = os.path.join(systfol, "MC-{}".format(ID))
    already_done = ut.status_ok(path=meta_mc["path"])
    if not already_done:
        try:
            if not os.path.exists(meta_mc["path"]):
                os.makedirs(meta_mc["path"])
            os.chdir(meta_mc["path"])
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
            already_done=True
            energies_file = os.path.join(meta_mc["path"], "energies.json")
            energies = ut.load_js(energies_file)
        except NotConvergedError:
            meta_mc["status"] = "FAIL"
        ut.save_status(meta_mc)
    #-----------------------------------------------------------------------------#
    # FDE-MP2  A in B
    #-----------------------------------------------------------------------------#   
    meta_mpa = dict(method_A="MP2", method_B="import", opt=None, status=None,basename="emb")
    meta_mpa["path"] = os.path.join(meta_mc["path"], "MP2_A")
    already_done_mpa = ut.status_ok(path=meta_mpa["path"])
    if not already_done_mpa and already_done:
        memory = 14000
        rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory, "fde": "true"})
        frag_specs = dict(frag_a=frags["A"], frag_b=frags["B"])
        if found_elconf:
            frag_specs.update(elconf)
        frag_str = Tfragments.substitute(Tdefaults["molecule"], **frag_specs)
        fde_sect = Tfde.substitute(Tdefaults["fde"], **{"method_a": "import_rhoA = true", "method_b": "import_rhoB = true"})
        extras = "\n".join([extra_basic]+[fde_sect])
        specs_mpa = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=frag_str, extras=extras)
        queue_mpa = dict(**slrm.shabug_XS)  
        if not os.path.exists(meta_mpa["path"]):
           os.makedirs(meta_mpa["path"]) 
        iterDir = ut.get_last_iter_dir(active="A", path=meta_mc["path"],opt="macrocycles")
        sh.copy(os.path.join(iterDir, "Densmat_B.txt"), meta_mpa["path"])
        copy_density(os.path.join(iterDir, "FDE_State0_tot_dens.txt"),
                    os.path.join(meta_mpa["path"], "Densmat_A.txt"),
                    header_src=False, alpha_only_src=False)
        json_files = gl.glob("*.json")
        ut.run_job(specs_mpa, queue_mpa, meta_mpa, Tinp, q_custom=slrm.slurm_add, batch_mode=True)  #TODO slurm_stuff
        ut.save_status(meta_mpa)
        njsf = [i for i in gl.glob("*.json") if i not in json_files][0]
        data = ut.load_js(njsf)  # default ccp json_filename
        energies.append(data["SCF"][-1]) # CHECK
        ut.dump_js(energies, energies_file)
    #-----------------------------------------------------------------------------#
    # FDE-MP2  B in A
    #-----------------------------------------------------------------------------#
    meta_mpb = dict(method_A="MP2", method_B="import", opt=None, status=None,basename="emb")
    meta_mpb["path"] = os.path.join(meta_mc["path"], "MP2_B")
    already_done_mpb = ut.status_ok(path=meta_mpb["path"])
    if not already_done_mpb and already_done:
        memory = 14000
        rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **rem_adc_basic, **{"memory": memory, "fde": "true"})
        frag_specs = dict(frag_a=frags["B"], frag_b=frags["A"])
        if found_elconf:
            frag_specs.update(elconf_rev)
        frag_str = Tfragments.substitute(Tdefaults["molecule"], **frag_specs)
        fde_sect = Tfde.substitute(Tdefaults["fde"], **{"method_a": "import_rhoA = true", "method_b": "import_rhoB = true"})
        extras = "\n".join([extra_basic]+[fde_sect])
        specs_mpb = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=frag_str, extras=extras)
        queue_mpb = dict(**slrm.shabug_XS)  
        if not os.path.exists(meta_mpb["path"]):
           os.makedirs(meta_mpb["path"]) 
#        iterDir = ut.get_last_iter_dir(active="A", path=meta["path"],opt="macrocycles")
        sh.copy(os.path.join(iterDir, "Densmat_B.txt"), os.path.join(meta_mpb["path"], "Densmat_A.txt"))
        copy_density(os.path.join(iterDir, "FDE_State0_tot_dens.txt"),
                    os.path.join(meta_mpb["path"], "Densmat_B.txt"),
                    header_src=False, alpha_only_src=False)
        ut.run_job(specs_mpb, queue_mpb, meta_mpb, Tinp, q_custom=slrm.slurm_add, batch_mode=True)  
        ut.save_status(meta_mpb)
        
