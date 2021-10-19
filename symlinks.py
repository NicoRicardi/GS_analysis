#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 21:00:16 2021

@author: nico
"""
import os
import glob as gl

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

systems = ["7HQ_2MeOH", "7HQ_formate", "Uracil_5H2O", "XVI_2HCOOH"]
for system in systems:
    for calc in ["FT-ME", "FT-SE"]:
        os.chdir(os.path.join(system, calc))
        cyfols = gl.glob("cy*")
        os.chdir("..")
        for cyfol in cyfols:
            it = int(cyfol[2:])
            if it % 2:
                continue
            src = os.path.join(calc, cyfol, "FT{}-MC-{}".format(it, calc[-2:]))
            dst = "FT{}-MC-{}".format(it, calc[-2:])
            symlinkif(src, dst)