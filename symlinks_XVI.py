#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 22:23:24 2021

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
        
system = "XVI_2HCOOH"
calc = "FT-SE"
os.chdir(os.path.join(system, calc))
cyfols = gl.glob("cy*")
for cyfol in cyfols:
    it = int(cyfol[2:])
    if it % 2 or it == 0:
        continue
    dst = os.path.join("FT{}-MC-{}".format(it, calc[-2:]), "MP2_B")
    src = os.path.join("FT{}-MC-{}".format(it-1, calc[-2:]), "MP2_A")
    symlinkif(src, dst)

