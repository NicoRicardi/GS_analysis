#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 12:17:52 2021

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

root = os.getcwd()
systems = ["7HQ_2MeOH", "7HQ_formate", "Uracil_5H2O", "XVI_2HCOOH"]
for system in systems:
    os.chdir(system)
    symlinkif(os.path.join("FT-SE", "cy0", "FT0-MC-SE"), "MC-SE-nopp")