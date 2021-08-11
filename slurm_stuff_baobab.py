#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 14:55:25 2019

@author: nico
"""
script = "QC5aut"
p_shabug = "shared-cpu,debug-cpu" 

debug_XS = dict(mem=8000, cpus=1, time=15, partition="debug-cpu",
               script=script)
debug_S = dict(mem=16000, cpus=2, time=15, partition="debug-cpu",
               script=script)
debug_M = dict(mem=32000, cpus=4, time=15, partition="debug-cpu",
               script=script)
debug_L = dict(mem=64000, cpus=8, time=15, partition="debug-cpu",
               script=script)

shabug_XS = dict(mem=8000, cpus=1, time=15, partition=p_shabug,
               script=script)
shabug_S = dict(mem=16000, cpus=2, time=15, partition=p_shabug,
               script=script)
shabug_M = dict(mem=32000, cpus=4, time=15, partition=p_shabug,
               script=script)
shabug_L = dict(mem=64000, cpus=8, time=15, partition=p_shabug,
               script=script)

shaweso_XS = dict(mem=8000, cpus=1, time=360, partition="shared-cpu,private-wesolowski-cpu",
               script=script)
shaweso_S = dict(mem=16000, cpus=2, time=360, partition="shared-cpu,private-wesolowski-cpu",
               script=script)
shaweso_M = dict(mem=32000, cpus=4, time=360, partition="shared-cpu,private-wesolowski-cpu",
               script=script)
shaweso_L = dict(mem=64000, cpus=8, time=360, partition="shared-cpu,private-wesolowski-cpu",
               script=script)
# small nodes on wesolowski partition
# last number means 1/n share of node, i.e. '...small2' := half node
weso_small1 = dict(mem=94000, cpus=20, time=180, partition="private-wesolowski-cpu",
               script=script)
weso_small2 = dict(mem=48000, cpus=10, time=180, partition="private-wesolowski-cpu",
               script=script)
weso_small3 = dict(mem=32000, cpus=8,  time=180, partition="private-wesolowski-cpu",
               script=script)
# big nodes on wesolowski partition
weso_big1 = dict(mem=768000, cpus=32, time="0-12", partition="private-wesolowski-cpu",
                 script=script)
weso_big2 = dict(mem=384000, cpus=16, time="0-12", partition="private-wesolowski-cpu",
                 script=script)
weso_big3 = dict(mem=256000, cpus=12, time="0-12", partition="private-wesolowski-cpu",
                 script=script)
weso_big4 = dict(mem=192000, cpus=8,  time="0-12", partition="private-wesolowski-cpu",
                 script=script)


slurm_add = {"--nodes": 1,
             "--ntasks-per-node": 1,
             "--mail-type": "FAIL",
             " --mail-user": "niccolo.ricardi@unige.ch"}
