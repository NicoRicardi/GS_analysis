#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 21:00:16 2021

@author: nico
"""
import os
import glob as gl
import shutil as sh

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
    elif printout:
        print("Nothing done")

def mkdif(folder):
    if not os.path.exists(folder):
        os.mkdir(folder)
        

root = os.getcwd()
systems = ["7HQ_2MeOH", "7HQ_formate", "Uracil_5H2O", "XVI_2HCOOH"]
for system in systems:
    for calc in ["FT-ME", "FT-SE"]:
        os.chdir(os.path.join(root, system, calc))
        cyfols = gl.glob("cy*")
        os.chdir(os.path.join(root, system))
        for cyfol in cyfols:
            it = int(cyfol[2:])
            if it % 2:
                continue
            src = os.path.join(calc, cyfol, "FT{}-MC-{}".format(it, calc[-2:]))
            dst_mc = "FT{}-MC-{}".format(it, calc[-2:])
            dst = "FT{}-{}".format(it, calc[-2:])
            if os.path.exists(dst_mc):
                sh.rmtree(dst_mc)
            if os.path.exists(dst):
                sh.rmtree(dst)   
            mkdif(dst)
            mkdif(dst_mc)
            MCcys = gl.glob(os.path.join(src,"cy*"))
            for cy in MCcys:
                n_cy = cy.split("cy")[-1]
                symlinkif(cy, os.path.join(dst_mc,"cy{}".format(n_cy)))
            bfol = os.path.join("..",  calc, "cy{}".format(it-1))
            symlinkif(bfol, os.path.join(dst_mc,"HF_B"))
            symlinkif(os.path.join("..", calc, "cy{}".format(it)), os.path.join(dst,"HF_A"))
            symlinkif(bfol, os.path.join(dst,"HF_B"))
#            symlinkif(src, dst)


