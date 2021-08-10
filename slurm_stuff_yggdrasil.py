#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 14:55:25 2019

@author: nico
"""
script = "QC5aut"

p_shared = "shared-cpu"
p_shabug = "shared-cpu,debug-cpu"  # do not add space!!
p_weso = "private-wesolowski-bigmem"
shared_XXS = dict(mem=10000, cpus=1, time=15, partition=p_shared,
               script=script)
shared_XS = dict(mem=16000, cpus=2, time=15, partition=p_shared,
               script=script)
shared_S = dict(mem=32000, cpus=3, time=15, partition=p_shared,
               script=script)
shared_M = dict(mem=64000, cpus=4, time=15, partition=p_shared,
               script=script)
shared_L = dict(mem=96000, cpus=8, time=15, partition=p_shared,
               script=script)
shared_XL = dict(mem=128000, cpus=16, time=15, partition=p_shared,
               script=script)
shared_XXL = dict(mem=188000, cpus=36, time=15, partition=p_shared,
               script=script)

shabug_XXS = dict(mem=10000, cpus=1, time=15, partition=p_shabug,
               script=script)
shabug_XS = dict(mem=16000, cpus=2, time=15, partition=p_shabug,
               script=script)
shabug_S = dict(mem=32000, cpus=3, time=15, partition=p_shabug,
               script=script)
shabug_M = dict(mem=64000, cpus=4, time=15, partition=p_shabug,
               script=script)
shabug_L = dict(mem=96000, cpus=8, time=15, partition=p_shabug,
               script=script)
shabug_XL = dict(mem=128000, cpus=16, time=15, partition=p_shabug,
               script=script)
shabug_XXL = dict(mem=188000, cpus=36, time=15, partition=p_shabug,
               script=script)

weso_big1 = dict(mem=1545830, cpus=36, time="0-12", partition=p_weso,
                 script=script)
weso_big2 = dict(mem=772915, cpus=18, time="0-12", partition=p_weso,
                 script=script)
weso_big3 = dict(mem=515276, cpus=12, time="0-12", partition=p_weso,
                 script=script)
weso_big4 = dict(mem=386457, cpus=9,  time="0-12", partition=p_weso,
                 script=script)
weso_big5 = dict(mem=309166, cpus=7,  time="0-2", partition=p_weso,
                 script=script)
weso_big6 = dict(mem=257638, cpus=6,  time="0-12", partition=p_weso,
                 script=script)
weso_big7 = dict(mem=220832, cpus=5,  time="0-12", partition=p_weso,
                 script=script)
weso_big8 = dict(mem=193228, cpus=4,  time="0-1", partition=p_weso,
                 script=script)

slurm_add = {"--nodes": 1,
            "--mail-type": "FAIL",
            " --mail-user": "niccolo.ricardi@unige.ch"}


