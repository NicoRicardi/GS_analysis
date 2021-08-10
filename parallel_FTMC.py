#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 16:49:29 2021

@author: nico
"""

import subprocess as sp
import sys

systems = ["7HQ_2MeOH", "7HQ_formate", "Uracil_5H2O", "XVI_2HCOOH"]
expansion = sys.argv[1]

for system in systems:
    slrm = "--time=12:00 -p shared-cpu --mem=8000 -c 1 --nodes=1 "
    script = "FTMC.py"
    command = "sbatch {slrm} python3 {script} {system} {expansion}".format(slrm=slrm, script=script, system=system, expansion=expansion) 
    print(command)
    p = sp.run(command, stdout=sp.PIPE, stderr=sp.STDOUT, shell=True)
    print(p.stdout)