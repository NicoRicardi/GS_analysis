#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 16:52:01 2022

@author: nico
"""
import re
fname = "220604_NR.tex"
#fname = "todel.txt"
with open(fname,"r") as f:
    text = f.read()

p = re.compile(r"\\nr{")

it = p.finditer(text)

points = []
for match in it:
    cnt = 1
    for i, c in enumerate(text[match.start()+5:]):
        if c == "{":
            cnt += 1
        if c == "}":
            cnt -= 1
        if cnt == 0:
            break
    points.extend([match.start(), match.start() + 6 + i])

new = [text[points[2*n] + 4:points[2*n + 1] - 1]+text[points[2*n + 1]:points[2*n + 2]] for n in range((len(points) // 2) - 1)] 

newfile = text[:points[0]] + "".join(new) + text[points[-2] + 4:points[-1] - 1] + text[points[-1]:]

with open("processed.tex","w") as f:
    f.write(newfile)